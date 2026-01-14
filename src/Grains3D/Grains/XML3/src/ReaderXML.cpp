#include "ReaderXML.hh"

#include <cassert>
#include <cstdlib>   // atoi, atof
#include <iostream>
#include <mutex>
#include <string>

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMLSParser.hpp>
#include <xercesc/dom/DOMConfiguration.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XMLUniDefs.hpp>

XERCES_CPP_NAMESPACE_USE
using namespace std;

// =====================================================================================
// Shared Xerces init/term (so ReaderXML + WriterXML can both use Xerces safely)
//
// If your WriterXML.cpp already declares these as `extern` and calls them,
// keep the names/signatures exactly.
// =====================================================================================
namespace grains_xml_internal {

static std::mutex g_mtx;
static int        g_refcount = 0;

void xerces_acquire()
{
  std::lock_guard<std::mutex> lock(g_mtx);
  if (g_refcount++ == 0) {
    XMLPlatformUtils::Initialize();
  }
}

void xerces_release()
{
  std::lock_guard<std::mutex> lock(g_mtx);
  if (g_refcount <= 0) return;
  if (--g_refcount == 0) {
    XMLPlatformUtils::Terminate();
  }
}

} // namespace grains_xml_internal

// =====================================================================================
// TU-local state (replaces the old `static DOMBuilder* ReaderXML::m_parser`)
// =====================================================================================
namespace {

DOMLSParser* g_parser = nullptr;

// Request the DOM “LS” implementation (same as old gLS logic)
static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };

// RAII for XMLCh* from transcode (prevents leaks)
struct XStr {
  XMLCh* p = nullptr;
  explicit XStr(const std::string& s) : p(XMLString::transcode(s.c_str())) {}
  ~XStr() { XMLString::release(&p); }
  const XMLCh* get() const { return p; }
  XStr(const XStr&) = delete;
  XStr& operator=(const XStr&) = delete;
};

// Convert XMLCh* -> std::string (and release temp)
std::string toString(const XMLCh* x)
{
  if (!x) return {};
  char* c = XMLString::transcode(x);
  std::string s = c ? c : "";
  XMLString::release(&c);
  return s;
}

static DOMElement* asElement(DOMNode* n)
{
  if (!n) return nullptr;

  if (n->getNodeType() == DOMNode::ELEMENT_NODE)
    return static_cast<DOMElement*>(n);

  // Common mistake: passing a DOMDocument*
  if (n->getNodeType() == DOMNode::DOCUMENT_NODE) {
    auto* doc = static_cast<DOMDocument*>(n);
    return doc ? doc->getDocumentElement() : nullptr;
  }

  return nullptr; // TEXT_NODE / COMMENT_NODE / etc.
}


} // namespace

// ----------------------------------------------------------------------------
// Initializes the reader (legacy semantics)
// ----------------------------------------------------------------------------
void ReaderXML::initialize()
{
  grains_xml_internal::xerces_acquire();

  if (!g_parser) {
    DOMImplementation* impl =
      DOMImplementationRegistry::getDOMImplementation(gLS);
    DOMImplementationLS* implLS = dynamic_cast<DOMImplementationLS*>(impl);

    if (!implLS) {
      throw std::runtime_error("Xerces: DOMImplementationLS not available");
    }

    // Xerces 3.x: DOMBuilder is gone; use LSParser instead
    g_parser = implLS->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, nullptr);

    // Old code set features on DOMBuilder; for LSParser use DOMConfiguration
    DOMConfiguration* cfg = g_parser->getDomConfig();

    const bool doNamespaces = false;
    const bool doSchema = false;
    const bool schemaFullChecking = false;

    // Set exactly the same intent as the old code.
    if (cfg->canSetParameter(XMLUni::fgDOMNamespaces, doNamespaces))
      cfg->setParameter(XMLUni::fgDOMNamespaces, doNamespaces);

    if (cfg->canSetParameter(XMLUni::fgXercesSchema, doSchema))
      cfg->setParameter(XMLUni::fgXercesSchema, doSchema);

    if (cfg->canSetParameter(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking))
      cfg->setParameter(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);

    // Old: m_parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
    if (cfg->canSetParameter(XMLUni::fgDOMDatatypeNormalization, true))
      cfg->setParameter(XMLUni::fgDOMDatatypeNormalization, true);
  }
}

// ----------------------------------------------------------------------------
// Frees the reader (legacy semantics)
// ----------------------------------------------------------------------------
void ReaderXML::terminate()
{
  // Old code only called Terminate() and leaked the parser. We can safely
  // release the parser while keeping external behavior unchanged.
  if (g_parser) {
    g_parser->release();
    g_parser = nullptr;
  }

  grains_xml_internal::xerces_release();
}

// ----------------------------------------------------------------------------
// Returns a node from its name and the root node
// Legacy: first descendant with tag name (not just direct child)
// ----------------------------------------------------------------------------
DOMNode* ReaderXML::getNode(DOMElement* root, string const& name)
{
  if (!root) return nullptr;
  XStr tag(name);
  DOMNodeList* lst = root->getElementsByTagName(tag.get());
  return lst ? lst->item(0) : nullptr;
}

// ----------------------------------------------------------------------------
// Returns a node from its name and the root node
// Legacy: search among filtered element children returned by getNodes(DOMNode*)
// ----------------------------------------------------------------------------
DOMNode* ReaderXML::getNode(DOMNode* root, string const& name)
{
  DOMNode* node = nullptr;
  DOMNodeList* nodes = ReaderXML::getNodes(root);

  for (XMLSize_t i = 0; nodes && i < nodes->getLength() && node == nullptr; i++) {
    if (name == ReaderXML::getNodeName(nodes->item(i)))
      node = nodes->item(i);
  }
  return node;
}

// ----------------------------------------------------------------------------
// Returns the scalar value of an attribute
// Legacy: assumes attribute exists; conversion uses atof()
// ----------------------------------------------------------------------------
double ReaderXML::getNodeAttr_Double(DOMNode* root, string const& name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XStr attrName(name);
  DOMNode* value = nodeValues ? nodeValues->getNamedItem(attrName.get()) : nullptr;

  // Legacy would likely crash here; keep behavior close but avoid UB where possible.
  assert(value != nullptr);

  const XMLCh* v = value->getNodeValue();
  // mimic: atof(XMLString::transcode(...))
  char* tmp = XMLString::transcode(v);
  double out = std::atof(tmp ? tmp : "");
  XMLString::release(&tmp);
  return out;
}

// ----------------------------------------------------------------------------
// Returns the integer value of an attribute (legacy: atoi)
// ----------------------------------------------------------------------------
int ReaderXML::getNodeAttr_Int(DOMNode* root, const string& name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XStr attrName(name);
  DOMNode* value = nodeValues ? nodeValues->getNamedItem(attrName.get()) : nullptr;

  assert(value != nullptr);

  const XMLCh* v = value->getNodeValue();
  char* tmp = XMLString::transcode(v);
  int out = std::atoi(tmp ? tmp : "");
  XMLString::release(&tmp);
  return out;
}

// ----------------------------------------------------------------------------
// Returns the string value of an attribute (legacy: XMLString::transcode)
// ----------------------------------------------------------------------------
// string ReaderXML::getNodeAttr_String(DOMNode* root, string const& name)
// {
//   DOMNamedNodeMap* nodeValues = root->getAttributes();
//   XStr attrName(name);
//   DOMNode* node = nodeValues ? nodeValues->getNamedItem(attrName.get()) : nullptr;

//   assert(node != nullptr);

//   const XMLCh* value = node->getNodeValue();
//   return toString(value);
// }

string ReaderXML::getNodeAttr_String(DOMNode* root, string const& name)
{
  DOMElement* elem = asElement(root);
  if (!elem) {
    std::cerr << "getNodeAttr_String: node is not an element (type="
              << (root ? root->getNodeType() : -1)
              << ", attr=" << name << ")\n";
    return "";
  }

  DOMNamedNodeMap* attrs = elem->getAttributes();
  if (!attrs) return "";

  XStr attrName(name);
  DOMNode* node = attrs->getNamedItem(attrName.get());
  if (!node) return "";

  return toString(node->getNodeValue());
}



// ----------------------------------------------------------------------------
// Returns whether the node has an attribute named name
// ----------------------------------------------------------------------------
bool ReaderXML::hasNodeAttr(DOMNode* root, string const& name)
{
  DOMNamedNodeMap* nodeValues = root ? root->getAttributes() : nullptr;
  if (!nodeValues) return false;

  XStr attrName(name);
  DOMNode* node = nodeValues->getNamedItem(attrName.get());
  return node != nullptr;
}

// ----------------------------------------------------------------------------
// Returns the name of a node (legacy: transcode nodeName)
// ----------------------------------------------------------------------------
string ReaderXML::getNodeName(DOMNode const* root)
{
  if (!root) return "";
  return toString(root->getNodeName());
}

// ----------------------------------------------------------------------------
// Returns the next node wrt to a node
// Legacy: returns the SECOND child (item(1)) of root->getChildNodes()
// ----------------------------------------------------------------------------
DOMNode* ReaderXML::getNodeNext(DOMNode* root)
{
  DOMNodeList* nodes = root ? root->getChildNodes() : nullptr;
  return nodes ? nodes->item(1) : nullptr;
}

// ----------------------------------------------------------------------------
// Returns the scalar value of the node "<root>xxx</root>" (legacy: firstChild)
// ----------------------------------------------------------------------------
double ReaderXML::getNodeValue_Double(DOMNode* root)
{
  if (!root || !root->getFirstChild()) return 0.0;

  const XMLCh* v = root->getFirstChild()->getNodeValue();
  char* tmp = XMLString::transcode(v);
  double out = std::atof(tmp ? tmp : "");
  XMLString::release(&tmp);
  return out;
}

// ----------------------------------------------------------------------------
// Returns the string value of the node "<root>xxx</root>" (legacy: firstChild)
// ----------------------------------------------------------------------------
string ReaderXML::getNodeValue_String(DOMNode* root)
{
  if (!root || !root->getFirstChild()) return "";
  return toString(root->getFirstChild()->getNodeValue());
}

// ----------------------------------------------------------------------------
// Returns the list of nodes in a node (legacy: descendants by tag name)
// ----------------------------------------------------------------------------
DOMNodeList* ReaderXML::getNodes(DOMElement* root, string const& name)
{
  if (!root) return nullptr;
  XStr tag(name);
  return root->getElementsByTagName(tag.get());
}

// ----------------------------------------------------------------------------
// Returns the list of nodes in a node
// Legacy: clone node, remove non-element children, return clone->getChildNodes()
// NOTE: This intentionally keeps the clone alive (leaks) like the old version.
// ----------------------------------------------------------------------------
DOMNodeList* ReaderXML::getNodes(DOMNode* root)
{
  if (!root) return nullptr;

  DOMNode* allNodes  = root->cloneNode(true);
  DOMNodeList* nodes = allNodes->getChildNodes();

  for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
    DOMNode* node = nodes->item(i);
    // Old code used nodeType != 1; ELEMENT_NODE matches that intention
    if (node->getNodeType() != DOMNode::ELEMENT_NODE) {
      allNodes->removeChild(node);
      // NOTE: live NodeList + removal can skip items; this matches legacy behavior.
    }
  }

  return allNodes->getChildNodes();
}

// ----------------------------------------------------------------------------
// Root node of the document (legacy: parseURI + getDocumentElement)
// ----------------------------------------------------------------------------
DOMElement* ReaderXML::getRoot(string const& xmlFile)
{
  DOMDocument* doc = nullptr;
  DOMElement*  root = nullptr;

  if (!g_parser) ReaderXML::initialize();

  try {
    doc  = g_parser->parseURI(xmlFile.c_str());
    root = doc ? doc->getDocumentElement() : nullptr;
  }
  catch (const DOMException& e) {
    XERCES_STD_QUALIFIER cerr << "XML exception " << e.code
                              << XERCES_STD_QUALIFIER endl
                              << toString(e.msg)
                              << XERCES_STD_QUALIFIER endl;
  }
  catch (const XMLException& e) {
    // Xerces 3.x often throws XMLException on parse errors; old code didn’t catch it,
    // but catching it here prevents hard termination and is usually what you want.
    XERCES_STD_QUALIFIER cerr << "XML exception: "
                              << toString(e.getMessage())
                              << XERCES_STD_QUALIFIER endl;
  }

  return root;
}

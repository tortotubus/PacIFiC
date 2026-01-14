#include "WriterXML.hh"

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMError.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/dom/DOMLocator.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/util/IOException.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XMLUniDefs.hpp>

XERCES_CPP_NAMESPACE_USE
using std::string;

// Keep the public interface the same: only m_document is in the header.
DOMDocument* WriterXML::m_document = nullptr;

// If you already have shared init/term with ReaderXML, keep using it.
// Otherwise you can replace these with direct Initialize/Terminate.
namespace grains_xml_internal {
void xerces_acquire();
void xerces_release();
}

namespace {

// Request DOM “LS” implementation
static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };

DOMImplementation* getImpl()
{
  return DOMImplementationRegistry::getDOMImplementation(gLS);
}

DOMImplementationLS* getImplLS(DOMImplementation* impl)
{
  return dynamic_cast<DOMImplementationLS*>(impl);
}

// RAII transcode helper (prevents leaks)
struct XStr {
  XMLCh* p = nullptr;
  explicit XStr(const std::string& s) : p(XMLString::transcode(s.c_str())) {}
  explicit XStr(const char* s) : p(XMLString::transcode(s)) {}
  ~XStr() { XMLString::release(&p); }
  const XMLCh* get() const { return p; }
  XStr(const XStr&) = delete;
  XStr& operator=(const XStr&) = delete;
};

std::string toString(const XMLCh* x)
{
  if (!x) return {};
  char* c = XMLString::transcode(x);
  std::string s = c ? c : "";
  XMLString::release(&c);
  return s;
}

// Serializer and error handler are private to this TU
DOMLSSerializer* g_serializer = nullptr;

struct DomErrHandler final : DOMErrorHandler {
  bool had_error = false;

  bool handleError(const DOMError& err) override {
    had_error = true;
    const DOMLocator* loc = err.getLocation();

    std::string uri  = (loc && loc->getURI()) ? toString(loc->getURI()) : "<unknown>";
    auto line = loc ? loc->getLineNumber() : 0;
    auto col  = loc ? loc->getColumnNumber() : 0;

    std::cerr << "Xerces DOM write error"
              << " uri=" << uri
              << " line=" << line
              << " col=" << col
              << " : " << toString(err.getMessage())
              << "\n";
    return true; // keep going; Xerces may still abort write
  }

  void resetErrors() { had_error = false; }
};

DomErrHandler g_errHandler;

} // namespace

// ----------------------------------------------------------------------------
DOMElement* WriterXML::initialize(string const& root)
{
  assert(root != "");

  grains_xml_internal::xerces_acquire();

  DOMImplementation* impl = getImpl();
  DOMImplementationLS* implLS = getImplLS(impl);
  if (!impl || !implLS) {
    grains_xml_internal::xerces_release();
    throw std::runtime_error("Xerces: DOMImplementation(LS) not available");
  }

  XStr rootName(root);
  m_document = impl->createDocument(nullptr, rootName.get(), nullptr);

  g_serializer = implLS->createLSSerializer();

  // Configure serializer
  if (DOMConfiguration* cfg = g_serializer->getDomConfig()) {
    if (cfg->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
      cfg->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);

    // Attach error handler so we see real diagnostics
    if (cfg->canSetParameter(XMLUni::fgDOMErrorHandler, &g_errHandler))
      cfg->setParameter(XMLUni::fgDOMErrorHandler, &g_errHandler);
  }

  return m_document->getDocumentElement();
}
void WriterXML::terminate(std::string const& file)
{
  assert(file != "");

  // IMPORTANT: don't terminate Xerces inside the try.
  // Let locals destruct while Xerces is still alive.
  try {
    DOMImplementation* impl = getImpl();
    DOMImplementationLS* implLS = getImplLS(impl);
    if (!impl || !implLS) {
      throw std::runtime_error("Xerces: DOMImplementation(LS) not available");
    }
    if (!m_document || !g_serializer) {
      throw std::runtime_error("WriterXML::terminate called before initialize()");
    }

    LocalFileFormatTarget target(file.c_str());

    DOMLSOutput* output = implLS->createLSOutput();
    output->setByteStream(&target);

    // Avoid XMLString::transcode + release entirely:
    static const XMLCh UTF8[] = { chLatin_U, chLatin_T, chLatin_F, chDigit_8, chNull };
    output->setEncoding(UTF8);

    g_serializer->write(m_document, output);

    output->release();

    g_serializer->release();
    g_serializer = nullptr;

    m_document->release();
    m_document = nullptr;
  }
  catch (const IOException& e) {
    std::cerr << "Xerces IOException writing " << file << ":\n"
              << "  " << toString(e.getMessage()) << "\n";
  }
  catch (const XMLException& e) {
    std::cerr << "Xerces XMLException writing " << file << ":\n"
              << "  " << toString(e.getMessage()) << "\n";
  }
  catch (const DOMException& e) {
    std::cerr << "Xerces DOMException writing " << file << " (code=" << e.code << "):\n"
              << "  " << toString(e.msg) << "\n";
  }
  catch (const std::exception& e) {
    std::cerr << "Exception writing " << file << ":\n"
              << "  " << e.what() << "\n";
  }

  // Now it's safe: all locals from try are destroyed already.
  grains_xml_internal::xerces_release();
}


// ----------------------------------------------------------------------------
DOMElement* WriterXML::createNode(string const& name)
{
  XStr tag(name);
  return m_document->createElement(tag.get());
}

DOMElement* WriterXML::createNode(DOMElement* root, string const& name)
{
  DOMElement* node = createNode(name);
  root->appendChild(node);
  return node;
}

void WriterXML::createNodeAttr(DOMElement* root, string const& attr, string const& value)
{
  XStr a(attr);
  XStr v(value);
  root->setAttribute(a.get(), v.get());
}

void WriterXML::createNodeAttr(DOMElement* root, string const& attr, double value)
{
  std::ostringstream os;
  os << value;
  createNodeAttr(root, attr, os.str());
}

void WriterXML::createNodeValue(DOMElement* root, string const& value)
{
  XStr txt(value);
  DOMText* node = m_document->createTextNode(txt.get());
  root->appendChild(node);
}

void WriterXML::createNodeValue(DOMElement* root, double const& value)
{
  std::ostringstream os;
  os << value;
  createNodeValue(root, os.str());
}

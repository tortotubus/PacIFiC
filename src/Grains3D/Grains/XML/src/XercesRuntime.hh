#ifndef GRAINS_XML_XERCES_RUNTIME_HH
#define GRAINS_XML_XERCES_RUNTIME_HH

// Shared with the legacy XML3 wrappers. Every native object must be destroyed
// before releasing its runtime reference. External Xerces users must retain
// their own Initialize/Terminate pairing while their objects are alive.
namespace grains_xml_internal {
void xerces_acquire();
void xerces_release();
struct XercesRuntime {
    XercesRuntime() { xerces_acquire(); }
    ~XercesRuntime() { xerces_release(); }
    XercesRuntime(const XercesRuntime&) = delete;
    XercesRuntime& operator=(const XercesRuntime&) = delete;
};
}
#endif

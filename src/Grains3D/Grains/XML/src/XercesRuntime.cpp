#include "XercesRuntime.hh"
#include <xercesc/util/PlatformUtils.hpp>
#include <mutex>
#include <cstddef>

namespace grains_xml_internal {
namespace {
std::mutex runtimeMutex;
std::size_t references = 0;
}
void xerces_acquire() {
    std::lock_guard<std::mutex> lock(runtimeMutex);
    // Increment only after successful initialization.
    if (references == 0) xercesc::XMLPlatformUtils::Initialize();
    ++references;
}
void xerces_release() {
    std::lock_guard<std::mutex> lock(runtimeMutex);
    if (references && --references == 0) xercesc::XMLPlatformUtils::Terminate();
}
}

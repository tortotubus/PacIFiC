#include <dlfcn.h>
#include <filesystem>
#include <vector>

std::filesystem::path get_library_dir() {
  Dl_info info;
  dladdr((void *)&get_library_dir, &info); // find our own symbol
  std::filesystem::path lib_path =
      std::filesystem::path(info.dli_fname).parent_path();
  return lib_path;
}

std::filesystem::path get_dtd_dir() {
  std::filesystem::path libdir = get_library_dir();
  std::vector<std::filesystem::path> candidates = {
      libdir / "../share/grains/Dtd", // install prefix
      libdir / "./Dtd"                  // in-tree build
  };

  for (auto &p : candidates) {
    std::error_code ec;
    if (std::filesystem::exists(p, ec))
      return std::filesystem::canonical(p, ec);
  }

  throw std::runtime_error("Cannot find DTD directory near " + libdir.string());
}

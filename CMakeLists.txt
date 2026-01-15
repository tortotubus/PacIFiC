# ================================================================================
#
# Top level project settings
#
# ================================================================================

cmake_minimum_required(VERSION 3.12.0...4.0.0)

project(PacIFiC 
  VERSION 1.0.0
  LANGUAGES C CXX
)

# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

add_subdirectory(src)


# ================================================================================
#
# Packaging 
#
# ================================================================================


set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")

# Keep junk out of the tarball
set(CPACK_SOURCE_IGNORE_FILES
  "/build/;/.git/;/.github/;/.idea/;/.vscode/;~$;[.]swp$;/install/;/_CPack_Packages/;CMakeCache.txt;"
)

include(CPack) 
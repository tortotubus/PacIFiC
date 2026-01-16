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

# ==========================================================
#
# Custom Build Configurations
#
# ==========================================================

# ==========================================================
# Custom build type: RelWithWarnings
# ==========================================================
set(CMAKE_CXX_FLAGS_RELWITHWARNINGS "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHWARNINGS "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "" FORCE)

# Common warning flags
set(_WARN_FLAGS
"-Werror=return-type \
-Werror=format-security \
-Werror=format=2 \
-Werror=nonnull \
-Werror=array-bounds \
-Werror=stringop-overflow \
-Werror=stringop-truncation \
-Werror=alloc-size-larger-than=9223372036854775807 \
-Werror=address \
-Werror=maybe-uninitialized"
)

# Optional runtime checks
set(_SANITIZER_FLAGS
    "-fsanitize=address,undefined -fno-omit-frame-pointer"
)

# Combine everything
set(CMAKE_CXX_FLAGS_RELWITHWARNINGS
    "${CMAKE_CXX_FLAGS_RELWITHWARNINGS} ${_WARN_FLAGS} ${_SANITIZER_FLAGS}"
    CACHE STRING "Flags used by the RelWithWarnings configuration" FORCE
)
set(CMAKE_C_FLAGS_RELWITHWARNINGS
    "${CMAKE_C_FLAGS_RELWITHWARNINGS} ${_WARN_FLAGS} ${_SANITIZER_FLAGS}"
    CACHE STRING "Flags used by the RelWithWarnings configuration" FORCE
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


# --- Utility target ---
add_custom_target(source_tar
  COMMAND "${CMAKE_CPACK_COMMAND}"
          --config "${CMAKE_BINARY_DIR}/CPackSourceConfig.cmake"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  COMMENT "Generating source tarball with CPack (TGZ)"
  VERBATIM
)
# ================================================================================
#
# Top level project settings
#
# ================================================================================

cmake_minimum_required(VERSION 3.12.0...4.0.0)

find_package(Git QUIET)

# ---- Version from git (tag -> SemVer, else fallback) ----
find_package(Git QUIET)

set(PACIFIC_FALLBACK_VERSION "0.0.1")
set(PACIFIC_VERSION_SEMVER   "${PACIFIC_FALLBACK_VERSION}")
set(PACIFIC_VERSION_FULL     "${PACIFIC_FALLBACK_VERSION}")

if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  # Prefer tags like v1.2.3 (or 1.2.3 if you don't use the v prefix)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --tags --match "v[0-9]*" --abbrev=12 --dirty --always
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  # Examples:
  #   v1.2.3
  #   v1.2.3-4-gabc1234d5678
  #   v1.2.3-4-gabc1234d5678-dirty
  #   abc1234d5678 (no tags reachable)
  
  set(PACIFIC_VERSION_FULL "${GIT_DESCRIBE}")

  # Extract SemVer from tag if present
  string(REGEX MATCH "^v([0-9]+\\.[0-9]+\\.[0-9]+)" _m "${GIT_DESCRIBE}")
  if(CMAKE_MATCH_1)
    set(PACIFIC_VERSION_SEMVER "${CMAKE_MATCH_1}")
  else()
    # No reachable SemVer tag: keep fallback numeric version
    # but still include commit id in FULL below.
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )

    if(GIT_SHA)
      set(PACIFIC_GIT_SHA ${GIT_SHA})
      set(PACIFIC_VERSION_FULL "${PACIFIC_VERSION_SEMVER}-g${GIT_SHA}")
    endif()


    string(TIMESTAMP PACIFIC_BUILD_TIME "%Y-%m-%d %H:%M:%S UTC" UTC)

    execute_process(
      COMMAND "${GIT_EXECUTABLE}" diff --quiet
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE GIT_DIRTY
      ERROR_QUIET
    )
    if(NOT GIT_DIRTY EQUAL 0)
      set(PACIFIC_VERSION_FULL "${PACIFIC_VERSION_FULL}-dirty")
    endif()
  endif()

  # If it *was* a tag-based describe, drop the leading "v" in FULL for readability
  string(REGEX REPLACE "^v" "" PACIFIC_VERSION_FULL "${PACIFIC_VERSION_FULL}")
endif()

message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "PacIFiC version ${PACIFIC_VERSION_FULL}")

project(PacIFiC 
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES C CXX
)
 

option(PACIFIC_ENABLE_MPI "Enable MPI-dependent targets" ON)

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
 
set(_SANITIZER_FLAGS
    "-fsanitize=address,undefined -fno-omit-frame-pointer"
)
 
set(CMAKE_CXX_FLAGS_RELWITHWARNINGS
    "${CMAKE_CXX_FLAGS_RELWITHWARNINGS} ${_WARN_FLAGS} ${_SANITIZER_FLAGS}"
    CACHE STRING "Flags used by the RelWithWarnings configuration" FORCE
)
set(CMAKE_C_FLAGS_RELWITHWARNINGS
    "${CMAKE_C_FLAGS_RELWITHWARNINGS} ${_WARN_FLAGS} ${_SANITIZER_FLAGS}"
    CACHE STRING "Flags used by the RelWithWarnings configuration" FORCE
)


# ------------------------------------------------------------------------------
# Handle optional MPI install prefix
# ------------------------------------------------------------------------------

include(GNUInstallDirs)

if(DEFINED CMAKE_MPI_INSTALL_LIBDIR AND CMAKE_MPI_INSTALL_LIBDIR)
  set(PACIFIC_MPI_INSTALL_LIBDIR "${CMAKE_MPI_INSTALL_LIBDIR}")
else()
  set(PACIFIC_MPI_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()


if(DEFINED CMAKE_MPI_INSTALL_BINDIR AND CMAKE_MPI_INSTALL_BINDIR)
  set(PACIFIC_MPI_INSTALL_BINDIR "${CMAKE_MPI_INSTALL_BINDIR}")
else()
  set(PACIFIC_MPI_INSTALL_BINDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
endif()

if(DEFINED CMAKE_MPI_INSTALL_INCLUDEDIR AND CMAKE_MPI_INSTALL_INCLUDEDIR)
  set(PACIFIC_MPI_INSTALL_INCLUDEDIR "${CMAKE_MPI_INSTALL_INCLUDEDIR}")
else()
  set(PACIFIC_MPI_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
endif()


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
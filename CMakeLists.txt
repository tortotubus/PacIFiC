# ================================================================================
#
# Top level project settings
#
# ================================================================================

cmake_minimum_required(VERSION 3.12.0...4.0.0)

# ---- Version selection policy (user -> git -> fallback) ----
find_package(Git QUIET)

set(PACIFIC_FALLBACK_VERSION "0.0.1")

# User-facing overrides (set via -D...=... or cmake-gui)
set(PACIFIC_VERSION_FULL_OVERRIDE   "" CACHE STRING "Override full version string (highest priority).")
set(PACIFIC_VERSION_SEMVER_OVERRIDE "" CACHE STRING "Override SemVer (e.g. 1.2.3). Used if FULL override is empty.")
option(PACIFIC_VERSION_USE_GIT "If no user override, compute version from git when available" ON)

# Defaults
set(PACIFIC_VERSION_SEMVER "${PACIFIC_FALLBACK_VERSION}")
set(PACIFIC_VERSION_FULL   "${PACIFIC_FALLBACK_VERSION}")

# 1) User-provided FULL override (highest priority)
if(PACIFIC_VERSION_FULL_OVERRIDE)
  set(PACIFIC_VERSION_FULL "${PACIFIC_VERSION_FULL_OVERRIDE}")

  # SemVer: prefer explicit semver override, else try to parse from FULL, else fallback
  if(PACIFIC_VERSION_SEMVER_OVERRIDE)
    set(PACIFIC_VERSION_SEMVER "${PACIFIC_VERSION_SEMVER_OVERRIDE}")
  else()
    string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _m "${PACIFIC_VERSION_FULL}")
    if(CMAKE_MATCH_1)
      set(PACIFIC_VERSION_SEMVER "${CMAKE_MATCH_1}")
    endif()
  endif()

# 2) User-provided SemVer override (next priority)
elseif(PACIFIC_VERSION_SEMVER_OVERRIDE)
  set(PACIFIC_VERSION_SEMVER "${PACIFIC_VERSION_SEMVER_OVERRIDE}")
  set(PACIFIC_VERSION_FULL   "${PACIFIC_VERSION_SEMVER_OVERRIDE}")

# 3) Git-provided (if enabled and available)
elseif(PACIFIC_VERSION_USE_GIT AND GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --tags --match "v[0-9]*" --abbrev=12 --dirty --always
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  if(GIT_DESCRIBE)
    string(REGEX REPLACE "^v" "" PACIFIC_VERSION_FULL "${GIT_DESCRIBE}")
    string(REGEX MATCH "^([0-9]+\\.[0-9]+\\.[0-9]+)" _m "${PACIFIC_VERSION_FULL}")
    if(CMAKE_MATCH_1)
      set(PACIFIC_VERSION_SEMVER "${CMAKE_MATCH_1}")
    else()
      # No reachable SemVer tag; keep fallback semver but still set a meaningful FULL
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      if(GIT_SHA)
        set(PACIFIC_VERSION_FULL "${PACIFIC_VERSION_SEMVER}-g${GIT_SHA}")
      endif()
    endif()
  endif()
# Fallback already set above
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

set(CMAKE_CXX_FLAGS_RELWITHWARNINGS "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELWITHWARNINGS  "${CMAKE_C_FLAGS_RELEASE}")

set(_WARN_FLAGS
  "-Werror=return-type -Werror=format-security -Werror=format=2 \
   -Werror=nonnull -Werror=array-bounds -Werror=stringop-overflow \
   -Werror=stringop-truncation \
   -Werror=alloc-size-larger-than=9223372036854775807 \
   -Werror=address -Werror=maybe-uninitialized"
)

set(_SANITIZER_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer")

string(APPEND CMAKE_CXX_FLAGS_RELWITHWARNINGS
       " ${_WARN_FLAGS} ${_SANITIZER_FLAGS}")
string(APPEND CMAKE_C_FLAGS_RELWITHWARNINGS
       " ${_WARN_FLAGS} ${_SANITIZER_FLAGS}")

mark_as_advanced(CMAKE_CXX_FLAGS_RELWITHWARNINGS CMAKE_C_FLAGS_RELWITHWARNINGS)


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
# Build and install third-party dependencies as part of CMake build
#
# ================================================================================

add_subdirectory(third_party)  

# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

add_subdirectory(src)
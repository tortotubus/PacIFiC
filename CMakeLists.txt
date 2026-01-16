# ================================================================================
#
# Top level project settings
#
# ================================================================================

cmake_minimum_required(VERSION 3.12.0...4.0.0)

project(PacIFiC 
  VERSION 0.0.1
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



if(DEFINED CMAKE_MPI_INSTALL_LIBDIR AND CMAKE_MPI_INSTALL_LIBDIR)
  set(PACIFIC_MPI_INSTALL_LIBDIR "${CMAKE_MPI_INSTALL_LIBDIR}")
else()
  set(PACIFIC_MPI_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()


if(DEFINED CMAKE_MPI_INSTALL_BINDIR AND CMAKE_MPI_INSTALL_BINDIR)
  set(PACIFIC_MPI_INSTALL_BINDIR "${CMAKE_MPI_INSTALL_BINDIR}")
else()
  set(PACIFIC_MPI_INSTALL_BINDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
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
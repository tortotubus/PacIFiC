if (PACIFIC_BUILD_THIRD_PARTY_HDF5)
    set(HDF5_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "HDF5 install prefix")
    add_subdirectory(hdf5)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_ZLIB)
    set(ZLIBNG_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for zlib-ng")
    add_subdirectory(zlib)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_XERCESC)
    set(XERECESC_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for XercesC")
    add_subdirectory(xercesc)
endif()


if (PACIFIC_BUILD_THIRD_PARTY_GTEST)
    set(GTEST_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for GTest")
    add_subdirectory(gtest)
endif()

if(PACIFIC_BUILD_THIRD_PARTY_PETSC)
    set(PACIFIC_PETSC_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "PETSc install prefix")
    add_subdirectory(petsc)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_LIBXML2)
    set(PACIFIC_LIBXML2_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "libxml2 install prefix")
    add_subdirectory(libxml2)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_EIGEN3)
    set(PACIFIC_EIGEN3_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Eigen3 install prefix")
    add_subdirectory(eigen)
endif()

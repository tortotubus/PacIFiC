

option(BUILD_THIRD_PARTY_HDF5 "Build HDF5" OFF)
set(HDF5_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "HDF5 install prefix")
if (BUILD_THIRD_PARTY_HDF5)
    add_subdirectory(hdf5)
endif()

option(BUILD_THIRD_PARTY_ZLIB "Build zlib" OFF)
set(ZLIBNG_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for zlib-ng")
if (BUILD_THIRD_PARTY_ZLIB)
    add_subdirectory(zlib)
endif()

option(BUILD_THIRD_PARTY_XERCESC "Build XercesC" OFF)
set(XERECESC_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for XercesC")
if (BUILD_THIRD_PARTY_XERCESC)
    add_subdirectory(xercesc)
endif()


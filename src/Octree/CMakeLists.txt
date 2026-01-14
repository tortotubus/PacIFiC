include(ExternalProject)

ExternalProject_Add(basilisk
	GIT_REPOSITORY        https://github.com/tortotubus/beam-basilisk.git
    UPDATE_DISCONNECTED   1
	SOURCE_DIR            ${CMAKE_BINARY_DIR}/_ext/basilisk-src/
	BINARY_DIR            ${CMAKE_BINARY_DIR}/_ext/basilisk-bin/
	INSTALL_DIR           ${CMAKE_BINARY_DIR}/_ext/basilisk-install/

	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
	-DCMAKE_BUILD_TYPE=Debug

	BUILD_COMMAND 
	${CMAKE_COMMAND}
	--build <BINARY_DIR>
	--config Debug
	--target qcc

	INSTALL_COMMAND
	${CMAKE_COMMAND}
	--install <BINARY_DIR>
	--prefix <INSTALL_DIR>
	--config Debug
	--component QCC
)

ExternalProject_Get_Property(basilisk install_dir)
ExternalProject_Get_Property(basilisk binary_dir)

set(QCC_EXECUTABLE "${install_dir}/bin/qcc" CACHE FILEPATH "Path to built qcc")


#
#
#


set(OCTREE_DLMFD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DLMFD" CACHE PATH "Octree DLMFD directory" FORCE)
set(OCTREE_VTKHYPERTREE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/VTKHyperTree" CACHE PATH "Octree VTKHyperTree directory" FORCE)

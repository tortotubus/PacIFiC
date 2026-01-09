# Paths to the real DTDs (choose whatever is correct for your tree)
set(CONVEX_DTD "${CMAKE_CURRENT_BINARY_DIR}/Convex.dtd")
set(CONTACT_DTD "${CMAKE_CURRENT_BINARY_DIR}/Contact.dtd")

# Generate the configured DTD into the build dir (or wherever you want)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Contact.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Contact.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Convex.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Convex.dtd"
  @ONLY
)

# Make them absolute + use forward slashes (important esp. on Windows) 
file(REAL_PATH "${CONVEX_DTD}" CONVEX_DTD)
file(REAL_PATH "${CONTACT_DTD}" CONTACT_DTD)
string(REPLACE "\\" "/" CONVEX_DTD "${CONVEX_DTD}")
string(REPLACE "\\" "/" CONTACT_DTD "${CONTACT_DTD}")


configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains2D_InFluid.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains2D_InFluid.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains2D_InFluid.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains2D_InFluid.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains2D.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains2D.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains3D_InFluid.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains3D_InFluid.dtd"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Grains3D.dtd.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Grains3D.dtd"
  @ONLY
)

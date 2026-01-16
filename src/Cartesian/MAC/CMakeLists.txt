cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(MAC
  VERSION 0.0.1
  LANGUAGES CXX
)

# ================================================================================
#
# USER PROVIDED DEPENDENCIES
#
# ================================================================================

find_package(MPI REQUIRED COMPONENTS C)

# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================
 
set(MAC_SOURCES "")

list(APPEND MAC_SOURCES
  ArrayVector/src/boolArray2D.cc
  ArrayVector/src/boolArray3D.cc
  ArrayVector/src/boolVector.cc
  ArrayVector/src/doubleArray2D.cc
  ArrayVector/src/doubleArray3D.cc
  ArrayVector/src/doubleArray4D.cc
  ArrayVector/src/doubleArray5D.cc
  ArrayVector/src/doubleArray6D.cc
  ArrayVector/src/doubleArray7D.cc
  ArrayVector/src/doubleVector.cc
  ArrayVector/src/intArray2D.cc
  ArrayVector/src/intArray3D.cc
  ArrayVector/src/intVector.cc
  ArrayVector/src/longLongIntArray3D.cc
  ArrayVector/src/longLongIntVector.cc
  ArrayVector/src/size_t_array2D.cc
  ArrayVector/src/size_t_array3D.cc
  ArrayVector/src/size_t_vector.cc
  ArrayVector/src/stringArray2D.cc
  ArrayVector/src/stringVector.cc
)
 
list(APPEND MAC_SOURCES
  Base/src/erf.cc
  Base/src/MAC_Application.cc
  Base/src/MAC_ApplicationRestorer.cc
  Base/src/MAC_ArithmeticExp.cc
  Base/src/MAC_ArrayExp.cc
  Base/src/MAC_assertions.cc
  Base/src/MAC_BalancedBinaryTree.cc
  Base/src/MAC_BalancedBinaryTreeIterator.cc
  Base/src/MAC_BalancedBinaryTreeNode.cc
  Base/src/MAC_BinStored.cc
  Base/src/MAC_BoolArray2D.cc
  Base/src/MAC_Bool.cc
  Base/src/MAC_BooleanExp.cc
  Base/src/MAC_BoolVector.cc
  Base/src/MAC.cc
  Base/src/MAC_Collection.cc
  Base/src/MAC_Communicator.cc
  Base/src/MAC_CommunicatorExp.cc
  Base/src/MAC_ComparisonExp.cc
  Base/src/MAC_ConditionalExp.cc
  Base/src/MAC_ConstantExp.cc
  Base/src/MAC_Container.cc
  Base/src/MAC_Context.cc
  Base/src/MAC_ContextPair.cc
  Base/src/MAC_ContextSimple.cc
  Base/src/MAC_ConvertTypeExp.cc
  Base/src/MAC_CutPointsExp.cc
  Base/src/MAC_Data.cc
  Base/src/MAC_DataWithContext.cc
  Base/src/MAC_DataWithContextExp.cc
  Base/src/MAC_DerivativeExp.cc
  Base/src/MAC_DistributedPartition.cc
  Base/src/MAC_DoubleArray2D.cc
  Base/src/MAC_DoubleArray3D.cc
  Base/src/MAC_Double.cc
  Base/src/MAC_DoubleComparator.cc
  Base/src/MAC_DoubleComparatorExact.cc
  Base/src/MAC_DoubleComparatorFloat.cc
  Base/src/MAC_DoubleVector.cc
  Base/src/MAC_Error.cc
  Base/src/MAC_Exec.cc
  Base/src/MAC_Expression.cc
  Base/src/MAC_ExternalAPI.cc
  Base/src/MAC_ExtractionExp.cc
  Base/src/MAC_FileToModule.cc
  Base/src/MAC_GroupExp.cc
  Base/src/MAC_HashTableSet.cc
  Base/src/MAC_HashTableSetIterator.cc
  Base/src/MAC_IndexSet.cc
  Base/src/MAC_IntArray2D.cc
  Base/src/MAC_IntArray3D.cc
  Base/src/MAC_Int.cc
  Base/src/MAC_InterpolExp.cc
  Base/src/MAC_IntVector.cc
  Base/src/MAC_Iterator.cc
  Base/src/MAC_KeyItemPair.cc
  Base/src/MAC_KeywordDataIterator.cc
  Base/src/MAC_KeywordDataPair.cc
  Base/src/MAC_Lex.cc
  Base/src/MAC_Lexical.cc
  Base/src/MAC_List.cc
  Base/src/MAC_ListIdentity.cc
  Base/src/MAC_ListItem.cc
  Base/src/MAC_ListIterator.cc
  Base/src/MAC_Map.cc
  Base/src/MAC_MapIterator.cc
  Base/src/MAC_MathFunctionExp.cc
  Base/src/MAC_MembershipExp.cc
  Base/src/MAC_MemoryTracer.cc
  Base/src/MAC_Module.cc
  Base/src/MAC_ModuleComparator.cc
  Base/src/MAC_ModuleExplorer.cc
  Base/src/MAC_ModuleIterator.cc
  Base/src/MAC_ModulePattern.cc
  Base/src/MAC_NumberedDoubleVectors.cc
  Base/src/MAC_Object.cc
  Base/src/MAC_ObjectReader.cc
  Base/src/MAC_ObjectRegister.cc
  Base/src/MAC_ObjectTest.cc
  Base/src/MAC_ObjectWriter.cc
  Base/src/MAC_Randomizer.cc
  Base/src/MAC_RelationalExp.cc
  Base/src/MAC_Root.cc
  Base/src/MAC_Sequence.cc
  Base/src/MAC_SequentialCommunicator.cc
  Base/src/MAC_Set.cc
  Base/src/MAC_SigalExp.cc
  Base/src/MAC_SortExp.cc
  Base/src/MAC_StringArray2D.cc
  Base/src/MAC_String.cc
  Base/src/MAC_StringExp.cc
  Base/src/MAC_StringVector.cc
  Base/src/MAC_System.cc
  Base/src/MAC_SystemExp.cc
  Base/src/MAC_Timer.cc
  Base/src/MAC_TransferExp.cc
  Base/src/MAC_UnaryArithmeticExp.cc
  Base/src/MAC_Variable.cc
  Base/src/MAC_VariableExp.cc
  Base/src/MAC_Vector.cc
  Base/src/MAC_VectorExp.cc
  Base/src/MAC_VectorIterator.cc
  Base/src/MAC_Yacc.cc
  Base/src/main.cc
)
 
list(APPEND MAC_SOURCES
  FiniteVolume/src/FV_BoundaryCondition.cc
  FiniteVolume/src/FV.cc
  FiniteVolume/src/FV_DiscreteField.cc
  FiniteVolume/src/FV_DiscreteField_Centered.cc
  FiniteVolume/src/FV_DiscreteField_Centered_SpatialSchemes.cc
  FiniteVolume/src/FV_DiscreteField_SpatialSchemes.cc
  FiniteVolume/src/FV_DiscreteField_Staggered.cc
  FiniteVolume/src/FV_DiscreteField_Staggered_SpatialSchemes.cc
  FiniteVolume/src/FV_DiscreteField_Tensor.cc
  FiniteVolume/src/FV_DiscreteField_Vertex.cc
  FiniteVolume/src/FV_DiscreteField_Vorticity.cc
  FiniteVolume/src/FV_DomainAndFields.cc
  FiniteVolume/src/FV_DomainBuilder.cc
  FiniteVolume/src/FV_MatlabPostProcessingWriter.cc
  FiniteVolume/src/FV_Mesh.cc
  FiniteVolume/src/FV_MorePostProcessing.cc
  FiniteVolume/src/FV_OneStepIteration.cc
  FiniteVolume/src/FV_ParaviewPostProcessingWriter.cc
  FiniteVolume/src/FV_PostProcessingWriter.cc
  FiniteVolume/src/FV_SplitSystem.cc
  FiniteVolume/src/FV_StepByStepProgression.cc
  FiniteVolume/src/FV_StreamFunction.cc
  FiniteVolume/src/FV_SystemNumbering.cc
  FiniteVolume/src/FV_TimeIteratorAdapter.cc
  FiniteVolume/src/FV_TimeIterator.cc
  FiniteVolume/src/FV_Vorticity.cc
  FiniteVolume/src/geomVector.cc
)
  
list(APPEND MAC_SOURCES
  LinearAlgebra/src/LA_BlockSeqMatrix.cc
  LinearAlgebra/src/LA_BlockSeqMatrixIterator.cc
  LinearAlgebra/src/LA.cc
  LinearAlgebra/src/LA_CRSmatrix.cc
  LinearAlgebra/src/LA_CRSmatrixIterator.cc
  LinearAlgebra/src/LA_DistImplementation.cc
  LinearAlgebra/src/LA_DistMatrix.cc
  LinearAlgebra/src/LA_DistScatter.cc
  LinearAlgebra/src/LA_DistVector.cc
  LinearAlgebra/src/LA_Implementation.cc
  LinearAlgebra/src/LA_Matrix.cc
  LinearAlgebra/src/LA_MatrixIterator.cc
  LinearAlgebra/src/LA_PairOfMatrixIterator.cc
  LinearAlgebra/src/LA_PelMatrix.cc
  LinearAlgebra/src/LA_PelMatrixIterator.cc
  LinearAlgebra/src/LA_Scatter.cc
  LinearAlgebra/src/LA_SeqImplementation.cc
  LinearAlgebra/src/LA_SeqMatrix.cc
  LinearAlgebra/src/LA_SeqScatter.cc
  LinearAlgebra/src/LA_SeqVector.cc
  LinearAlgebra/src/LA_ShiftedIndexMatrixIterator.cc
  LinearAlgebra/src/LA_Solver.cc
  LinearAlgebra/src/LA_StorableVectors.cc
  LinearAlgebra/src/LA_Vector.cc
)
 
# list(APPEND MAC_SOURCES
#   ExternalAPI/PETSc_3.23.4/src/EXT_PETScAPI.cc
#   ExternalAPI/PETSc_3.23.4/src/EXT_PETScImplementation.cc
#   ExternalAPI/PETSc_3.23.4/src/EXT_PETScMatrix.cc
#   ExternalAPI/PETSc_3.23.4/src/EXT_PETScScatter.cc
#   ExternalAPI/PETSc_3.23.4/src/EXT_PETScSolver.cc
#   ExternalAPI/PETSc_3.23.4/src/EXT_PETScVector.cc
# )

add_library(MAC SHARED ${MAC_SOURCES})
add_library(PacIFiC::MAC ALIAS MAC)

set_target_properties(MAC
  PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)


target_link_libraries(MAC 
  PUBLIC
    MPI::MPI_C
)

include(GNUInstallDirs)

target_include_directories(MAC
  PUBLIC
    # Build tree include dirs (absolute paths OK here)
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ArrayVector/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Base/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolume/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/include>

    # Install tree include dir (must NOT be in the source tree)
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/MAC>
)


# ================================================================================
#
# Install 
#
# ================================================================================

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS MAC
  EXPORT MACTargets
  LIBRARY DESTINATION ${PACIFIC_MPI_INSTALL_LIBDIR}     
  COMPONENT MAC_Runtime
  RUNTIME DESTINATION ${PACIFIC_MPI_INSTALL_BINDIR}      
  COMPONENT MAC_Runtime
  INCLUDES DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}  
)

# ---- 2. Install all public headers -------------------------------------------------
# If your project has multiple include roots, add them all.
install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/ArrayVector/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Base/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolume/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/include/
  DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}/MAC
  COMPONENT MAC_Devel
  FILES_MATCHING PATTERN "*.h*"      # match .h, .hh, .hpp
)

# ---- 4. Generate CMake package configuration --------------------------------------
set(MAC_INSTALL_CMAKEDIR ${PACIFIC_MPI_INSTALL_LIBDIR}/cmake/MAC)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MACConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfig.cmake
  INSTALL_DESTINATION ${MAC_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# ---- 5. Install the package + targets ---------------------------------------------
install(EXPORT MACTargets
  FILE MACTargets.cmake
  NAMESPACE PacIFiC::
  DESTINATION ${MAC_INSTALL_CMAKEDIR}
  COMPONENT MAC_Devel
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfigVersion.cmake
  DESTINATION ${MAC_INSTALL_CMAKEDIR}
  COMPONENT MAC_Devel
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT MACTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/MACTargets.cmake
  NAMESPACE PacIFiC::
)

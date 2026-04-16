#include <gtest/gtest.h>

#include <FV_PostProcessingWriter.hh>
#include <FV_ParaviewPostProcessingWriterVTKHDF.hh>

// Forward declarations used by the plugin factory API.
class MAC_ModuleExplorer;
class MAC_Communicator;
class FV_Mesh;
class FV_DiscreteField;
class FV_TimeIterator;

#include <list>
#include <string>

namespace {

class FVParaviewPostProcessingWriterVTKHDFTest : public ::testing::Test {
protected:
  void SetUp() override {
    // TODO(colive):
    // 1) Build a minimal MAC/FV runtime context.
    // 2) Create MAC_ModuleExplorer with:
    //      postprocessing_writer = "paraview_vtkhdf"
    //      results_directory / files_rootname
    // 3) Provide a valid communicator, mesh, and field list.
  }

  void TearDown() override {
    // TODO(colive): Destroy runtime objects created in SetUp().
  }

  // TODO(colive): Populate these when the integration fixture is wired.
  MAC_ModuleExplorer const *exp_ = 0;
  MAC_Communicator const *com_ = 0;
  FV_Mesh const *mesh_ = 0;
  std::list<FV_DiscreteField const *> fields_;
};

TEST_F(FVParaviewPostProcessingWriterVTKHDFTest,
       PluginFactoryCreatesVTKHDFWriter) {
  GTEST_SKIP() << "Skeleton test: wire exp_/com_/mesh_/fields_ first.";

  // Example target flow (enable once fixture is ready):
  // FV_PostProcessingWriter* writer = FV_PostProcessingWriter::make(
  //     /*owner=*/0, "paraview_vtkhdf", exp_, com_, fields_, mesh_,
  //     /*binary=*/true);
  // ASSERT_NE(writer, nullptr);
  // delete writer;
}

TEST_F(FVParaviewPostProcessingWriterVTKHDFTest, WriteCycleSmoke) {
  GTEST_SKIP() << "Skeleton test: requires valid FV_TimeIterator + mesh/fields.";

  // Example target flow (enable once fixture is ready):
  // FV_PostProcessingWriter* writer = FV_PostProcessingWriter::make(
  //     /*owner=*/0, "paraview_vtkhdf", exp_, com_, fields_, mesh_,
  //     /*binary=*/true);
  // ASSERT_NE(writer, nullptr);
  //
  // FV_TimeIterator const* t_it = ...;
  // writer->write_cycle(t_it, /*cycle_number=*/0);
  //
  // TODO(colive): Assert output file exists and contains /VTKHDF group.
  // delete writer;
}

} // namespace


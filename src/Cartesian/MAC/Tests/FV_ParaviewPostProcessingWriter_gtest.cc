#include <gtest/gtest.h>

#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <FV_PostProcessingWriter.hh>
#include <FV_TimeIterator.hh>

#include <MAC_Double.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_IntVector.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>

#include <doubleVector.hh>
#include <intVector.hh>

#include <filesystem>
#include <system_error>
#include <list>
#include <sstream>
#include <string>

namespace {

class FVParaviewPostProcessingWriterSmokeTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Ensure the built-in communicator is available.
    com_ = MAC_Exec::communicator();
    ASSERT_NE(com_, nullptr);

    // Create a dedicated output directory for this test run.
    out_dir_ =
        (std::filesystem::current_path() / "mac_paraview_writer_smoke")
            .string();
    if (com_->rank() == 0) {
      std::error_code ec;
      std::filesystem::remove_all(out_dir_, ec);
      ec.clear();
      std::filesystem::create_directories(out_dir_, ec);
    }
    com_->barrier();
  }

  void TearDown() override {
    if (writer_ != nullptr)
      writer_->destroy();
    for (FV_DiscreteField *f : fields_owned_) {
      if (f != nullptr)
        f->destroy();
    }
    fields_owned_.clear();
    if (time_it_ != nullptr)
      time_it_->destroy();
    if (mesh_ != nullptr)
      mesh_->destroy();
    if (writer_exp_ != nullptr)
      writer_exp_->destroy();
    if (time_exp_ != nullptr)
      time_exp_->destroy();
    if (mesh_exp_ != nullptr)
      mesh_exp_->destroy();
    if (writer_mod_ != nullptr)
      writer_mod_->destroy();
    if (time_mod_ != nullptr)
      time_mod_->destroy();
    if (mesh_mod_ != nullptr)
      mesh_mod_->destroy();
    
    // std::filesystem::remove_all(out_dir_);
  }

  void build_minimal_mesh() {
    mesh_mod_ = MAC_Module::create(MAC_Root::object(), "MeshRoot");
    MAC_Module *fv_mesh = MAC_Module::create(mesh_mod_, "FV_Mesh");
    MAC_Module *split = MAC_Module::create(mesh_mod_, "splitting_strategy");

    // 5 points per direction -> 4x4 cells
    doubleVector x(5);
    x(0) = 0.0;
    x(1) = 0.25;
    x(2) = 0.50;
    x(3) = 0.75;
    x(4) = 1.0;
    doubleVector y(5);
    y(0) = 0.0;
    y(1) = 0.25;
    y(2) = 0.50;
    y(3) = 0.75;
    y(4) = 1.0;

    fv_mesh->add_entry("vertices_coordinate_0",
                       MAC_DoubleVector::create(fv_mesh, x));
    fv_mesh->add_entry("vertices_coordinate_1",
                       MAC_DoubleVector::create(fv_mesh, y));
    mesh_mod_->add_module(fv_mesh);

    intVector dd(2);
    dd(0) = static_cast<int>(com_->nb_ranks());
    dd(1) = 1;
    split->add_entry("domain_decomposition", MAC_IntVector::create(split, dd));
    split->add_entry("security_bandwidth", MAC_Int::create(split, 1));
    mesh_mod_->add_module(split);

    mesh_exp_ = MAC_ModuleExplorer::create(MAC_Root::object(), mesh_mod_);
    mesh_ = FV_Mesh::create(MAC_Root::object(), mesh_exp_, 2);
  }

  void build_fields() {
    // Scalar centered field.
    FV_DiscreteField *phi = FV_DiscreteField::create(
        MAC_Root::object(), mesh_, "phi", "centered", /*nb_components=*/1,
        /*depth=*/1);
    phi->set_postprocessing_options("at_cell_centers", "phi");
    phi->build_field_numbering();
    const size_t nx_phi = phi->get_local_nb_dof(0, 0);
    const size_t ny_phi = phi->get_local_nb_dof(0, 1);
    for (size_t i = 0; i < nx_phi; ++i) {
      for (size_t j = 0; j < ny_phi; ++j) {
        if (phi->DOF_is_unknown(i, j, 0, 0)) {
          phi->set_DOF_value(i, j, 0, 0, 0, static_cast<double>(i + 10 * j));
        }
      }
    }
    fields_owned_.push_back(phi);

    // Vector centered field (2 comps, writer will pad to 3 for VTK vector).
    FV_DiscreteField *u = FV_DiscreteField::create(
        MAC_Root::object(), mesh_, "u", "centered", /*nb_components=*/2,
        /*depth=*/1);
    u->set_postprocessing_options("at_cell_centers", "u");
    u->build_field_numbering();
    const size_t nx_u = u->get_local_nb_dof(0, 0);
    const size_t ny_u = u->get_local_nb_dof(0, 1);
    for (size_t i = 0; i < nx_u; ++i) {
      for (size_t j = 0; j < ny_u; ++j) {
        if (u->DOF_is_unknown(i, j, 0, 0)) {
          u->set_DOF_value(i, j, 0, 0, 0, static_cast<double>(i));
        }
        if (u->DOF_is_unknown(i, j, 0, 1)) {
          u->set_DOF_value(i, j, 0, 1, 0, static_cast<double>(-j));
        }
      }
    }
    fields_owned_.push_back(u);
  }

  void build_writer_config() {
    writer_mod_ = MAC_Module::create(MAC_Root::object(), "FV_ResultSaver");
    writer_mod_->add_entry("results_directory",
                           MAC_String::create(writer_mod_, out_dir_));
    writer_mod_->add_entry("files_rootname",
                           MAC_String::create(writer_mod_, "paraview_smoke"));
    writer_exp_ = MAC_ModuleExplorer::create(MAC_Root::object(), writer_mod_);
  }

  void build_time_iterator() {
    time_mod_ = MAC_Module::create(MAC_Root::object(), "FV_TimeIterator");
    time_mod_->add_entry("time_start", MAC_Double::create(time_mod_, 0.0));
    time_mod_->add_entry("time_end", MAC_Double::create(time_mod_, 1.0));
    time_mod_->add_entry("time_step", MAC_Double::create(time_mod_, 1.0));
    time_exp_ = MAC_ModuleExplorer::create(MAC_Root::object(), time_mod_);
    time_it_ = FV_TimeIterator::create(MAC_Root::object(), time_exp_, 0.0);
  }

  std::string out_dir_;
  MAC_Communicator const *com_ = nullptr;

  MAC_Module *mesh_mod_ = nullptr;
  MAC_Module *time_mod_ = nullptr;
  MAC_Module *writer_mod_ = nullptr;

  MAC_ModuleExplorer *mesh_exp_ = nullptr;
  MAC_ModuleExplorer *time_exp_ = nullptr;
  MAC_ModuleExplorer *writer_exp_ = nullptr;

  FV_Mesh *mesh_ = nullptr;
  FV_TimeIterator *time_it_ = nullptr;
  FV_PostProcessingWriter *writer_ = nullptr;
  std::vector<FV_DiscreteField *> fields_owned_;
};

TEST_F(FVParaviewPostProcessingWriterSmokeTest, WritesOneCycleToVtrAndPvd) {
  build_minimal_mesh();
  build_fields();
  build_writer_config();
  build_time_iterator();

  std::list<FV_DiscreteField const *> fields;
  for (FV_DiscreteField *f : fields_owned_) {
    fields.push_back(f);
  }
  writer_ = FV_PostProcessingWriter::make(MAC_Root::object(), "paraview",
                                          writer_exp_, com_, fields, mesh_,
                                          /*a_binary=*/false);
  ASSERT_NE(writer_, nullptr);

  writer_->write_cycle(time_it_, /*cycle_number=*/0);

  // const std::filesystem::path pvd =
  //     std::filesystem::path(out_dir_) / "paraview_smoke.pvd";
  // const std::filesystem::path vtr =
  //     std::filesystem::path(out_dir_) / "paraview_smokeT0.vtr";

  // EXPECT_TRUE(std::filesystem::exists(pvd)) << pvd.string();
  // EXPECT_TRUE(std::filesystem::exists(vtr)) << vtr.string();
}

} // namespace

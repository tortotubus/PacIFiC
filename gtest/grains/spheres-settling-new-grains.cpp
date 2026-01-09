#include <gtest/gtest.h>

#include "Grains.hh"
#include "GrainsBuilderFactory.hh"
#include "ReaderXML.hh"
#include <cstdlib>
#include <mpi.h>
#include <string>

#include <filesystem>
#include <iostream>

TEST(GrainsMPI, SpheresSettling) {

  std::string test_name = "spheres-settling-new-grains";

  int initialized = 0;
  MPI_Initialized(&initialized);
  ASSERT_TRUE(initialized) << "MPI must be initialized in test main()";

  int rankproc = 0, nprocs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::filesystem::path filename_path =
      std::filesystem::path(GRAINS_TEST_DATA_DIR) / test_name /
      "Grains/Init/insert.xml";
  std::string filename = filename_path.string();
  std::string filename_exe;

  // Create a temporary input file with the proper XML header
  filename_exe = GrainsBuilderFactory::init(filename, rankproc, nprocs);

  // Creates the Grains application
  ReaderXML::initialize();
  DOMElement *rootNode = ReaderXML::getRoot(filename_exe);
  string option = ReaderXML::getNodeAttr_String(rootNode, "Type");

  Grains *grains = NULL;
  grains = GrainsBuilderFactory::create(rootNode);

  // Initial output message
  grains->initialOutputMessage();

  // Tasks to perform before time-stepping
  grains->do_before_time_stepping(rootNode);
  ReaderXML::terminate();

  // Delete the temporary input file
  // if (rankproc == 0) {
    // string cmd = "/bin/rm " + filename_exe;
    // GrainsExec::m_return_syscmd = system(cmd.c_str());
  // }

  // Run the simulation
  grains->Simulation();

  // Tasks to perform after time-stepping
  grains->do_after_time_stepping();

  // Delete the Grains application
  delete grains;

  if (Component::getNbCreatedComponents())
    cout << "Warning: " << Component::getNbCreatedComponents()
         << " component(s) was/were not properly destroyed on proc " << rankproc
         << endl;
}

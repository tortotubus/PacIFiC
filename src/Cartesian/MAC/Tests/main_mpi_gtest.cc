#include <gtest/gtest.h>
#include <mpi.h>

#include <MAC_Bool.hh>
#include <MAC_Context.hh>
#include <MAC_Exec.hh>
#include <MAC_Variable.hh>

int main(int argc, char **argv) {
  int initialized = 0;
  MPI_Initialized(&initialized);

  const bool owns_mpi = (initialized == 0);
  if (owns_mpi) {
    MPI_Init(&argc, &argv);
  }

  // Match normal MAC execution mode: mark MPI as available in MAC context.
  MAC_Variable const *with_mpi = MAC_Variable::object("BS_with_MPI");
  MAC_Context const *ctx = MAC_Exec::execution_context();
  if (!ctx->has_variable(with_mpi)) {
    MAC_Exec::add_variable_to_execution_context(with_mpi,
                                                MAC_Bool::create(0, true));
  }

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  int finalized = 0;
  MPI_Finalized(&finalized);
  if (owns_mpi && !finalized) {
    MPI_Finalize();
  }

  return result;
}

#include <mpi.h>
#include <gtest/gtest.h>

int main(int argc, char** argv)
{    
    setenv("GRAINS_HOME", GRAINS_HOME, 1);  // overwrite=1

    // Initialize MPI before GoogleTest
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}

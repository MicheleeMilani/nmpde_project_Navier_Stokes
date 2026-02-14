#include "NavierStokes3D.hpp"
#include <filesystem>
namespace fs = std::filesystem;

// Main function.
int main(int argc, char *argv[])
{

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string param_file = argc > 1 ? argv[1] : "../src/case_3D.prm";

  dealii::Timer timer;
  // Start the timer
  timer.restart();

  NavierStokes problem;

  problem.parse_parameters(param_file);
  problem.setup();
  problem.solve();

  // Stop the timer
  timer.stop();

  // Output the elapsed time
  if(rank == 0)
    std::cout << "Time taken to solve Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;

  return 0;
}

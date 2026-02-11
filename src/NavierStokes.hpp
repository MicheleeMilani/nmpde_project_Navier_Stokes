#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include "Preconditioners.hpp"
#include "IncludesFile.hpp"

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class NavierStokes
{
public:
  // Physical dimension (2D)
  static constexpr unsigned int dim = 2;


  NavierStokes() 
                  : 
                  mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), 
                  mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), 
                  pcout(std::cout, mpi_rank == 0), 
                  inlet_velocity(), 
                  mesh(MPI_COMM_WORLD)
                  
  {}

    // Parse parameters from a file.
  void parse_parameters(const std::string &parameter_file);

  // Setup system.
  void
  setup();

  // Solve system.
  void
  solve();


  std::vector<double> vec_drag;
  std::vector<double> vec_lift;
  std::vector<double> vec_drag_coeff;
  std::vector<double> vec_lift_coeff;

  std::vector<double> time_prec;
  std::vector<double> time_solve;


  // Function for inlet velocity based on 3 tests
  class InletVelocity : public Function<dim>
  {
  public:

  InletVelocity()
    : Function<dim>(dim + 1)
  {
  }
    
    virtual void
    vector_value(const Point<dim> & p, Vector<double> &values) const override
    {
      values[0] = 4.0 * u_m * p[1] * (H - p[1])/ (H*H);
      
      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }
    
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component == 0) {
        return 4.0 * u_m * p[1] * (H - p[1]) / (H*H);
      }
      else
        return 0;
    }
    
    double getMeanVelocity() const
    {
      return 2.0 * u_m / 3.0;
    }
    
  protected:
    double H = 0.41;
    double u_m = 1.5;
  };

protected:

  // Declare parameters.
  void declare_parameters(ParameterHandler &prm);

  // Assemble system the first time to create mass-stiffness matrixes 
  void assemble(const double &time);
  // Assemble at each time step to only compute the convection matrix that changes overtime
  void assemble_time_step(const double &time);

  // Solve the problem for one time step.
  void solve_time_step(double time);

  // Output results.
  void output(const unsigned int &time_step, std::vector<double>) const;

  //Compute drag and lift
  std::vector<double> compute_lift_drag();

  void compute_pressure_difference();

  // MPI parallel. /////////////////////////////////////////////////////////////
  
  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Reynolds number.
  double Re;

  // Kinematic viscosity [m2/s].
  double nu;

  // Density
  double rho;

  // Forcing term.
  Functions::ZeroFunction<dim> forcing_term;
 
  // Final time.
  double T;

  double drag;
  double lift;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  std::string mesh_file_name;

  // Polynomial degree used for velocity.
  unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  unsigned int degree_pressure;

  // TIme step.
  double deltat;

  // g(x).
  Functions::ZeroFunction<dim> function_g;

  // h(x).
  Functions::ZeroFunction<dim> function_h;

  // Initial condition.
  Functions::ZeroFunction<dim> u_0;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;
  // Stiffness Matrix
  TrilinosWrappers::BlockSparseMatrix stiffness_matrix;
  // Mass Matrix
  TrilinosWrappers::BlockSparseMatrix mass_matrix;
  // Convection matrix at step (u_k)
  TrilinosWrappers::BlockSparseMatrix convection_matrix;

  // Pressure mass matrix, needed for preconditioning, also referred as Mp
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  TrilinosWrappers::MPI::BlockVector previous_solution;

};

#endif

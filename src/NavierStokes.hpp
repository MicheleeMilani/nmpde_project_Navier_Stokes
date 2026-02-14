#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include "Preconditioners.hpp"

#include <fstream>
#include <filesystem>
#include <iostream>
#include <mpi.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class NavierStokes
{
public:
  // Physical dimension (2D)
  static constexpr unsigned int dim = 2;

  NavierStokes() : 
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
    InletVelocity() : Function<dim>(dim + 1), u_m(0.3), test_case(1) {} // deafult values

    void initialize(unsigned int test_id, double u_m_val) {
      test_case = test_id;
      u_m = u_m_val;
    }

    double getMeanVelocity() const {
      return 2.0 * u_m / 3.0;
    }

    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override {
      double factor = 1.0;
      if (test_case == 3) {
        factor = std::sin(M_PI * this->get_time() / 8.0);
      }
      values[0] = 4.0 * u_m * p[1] * (H - p[1]) / (H * H) * factor;
      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double value(const Point<dim> &p, const unsigned int component = 0) const override {
      if (component == 0) {
        double factor = (test_case == 3) ? std::sin(M_PI * this->get_time() / 8.0) : 1.0;
        return 4.0 * u_m * p[1] * (H - p[1]) / (H * H) * factor;
      }
      return 0.0;
    }

    unsigned int get_test_case() const { return test_case; }

  protected:
    double H = 0.41;
    double u_m;
    unsigned int test_case;
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

  // Drag and lift.
  double drag;
  double lift;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  std::string mesh_file_name;

  // Polynomial degree used for velocity.
  unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  unsigned int degree_pressure;

  // Preconditioner type.
  unsigned int preconditioner_type;

  // TIme step.
  double deltat;
  unsigned int save_frequency;

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

#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

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


using namespace dealii;

  class PreconditionBlockIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
    }

  protected:
  };

// SIMPLE preconditioner.
    //
    // The application of the P_SIMPLE can be divided into two steps:
    // 1. Solve the block lower triangular system:
    //      [ F    0 ] [ sol1_u ] = [ src_u ]
    //      [ B   -S ] [ sol1_p ]   [ src_p ]
    // 2. Solve the subsequent system to correct the intermediate solution:
    //      [ I    D^-1*B^T ] [ dst_u ] = [ sol1_u ]
    //      [ 0      alpha  ] [ dst_p ]   [ sol1_p ]
    //
  class PreconditionSIMPLE 
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::MPI::BlockVector &sol_owned
               )
    {
      F = &F_;
      B = &B_; 
      B_T = &B_t;
    
      neg_diag_D_inv.reinit(sol_owned.block(0));
      diag_D_inv.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D_inv.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D_inv[i] = 1.0 / temp;
        neg_diag_D_inv[i] = - 1.0 / temp;     
      }

      // Create S_tilde =B * (D^-1) * B^T,
      // note: Using negative (D^-1) to create - S_tilde 
      B->mmult(negative_S_tilde, *B_T, neg_diag_D_inv); 

      // Initialize the preconditioners
      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(negative_S_tilde);
      
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const 
    {
      const unsigned int maxiter = 10000;
      const double tol = 1e-2;
      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());

      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      // 1. Solve the block lower triangular system:
    //      [ F    0 ]      [ sol1_u ] =       [ src_u ]
    //      [ B   -S_tilde ] [ sol1_p ]   [ src_p ]

      //Step 1.1 Solve Fsol1_u = src_u

      // Store in temporaries the results
      TrilinosWrappers::MPI::Vector sol1_u = src.block(0);  
      TrilinosWrappers::MPI::Vector sol1_p = src.block(1);
      
      TrilinosWrappers::MPI::Vector temp_1 = src.block(1);

      solver_gmres.solve(*F, sol1_u, src.block(0), preconditioner_F);

      B->vmult(temp_1, sol1_u); //temp_1 = B * sol1_u
      temp_1 -= src.block(1); //temp1 = src_p - B * sol1_u

      //Step 1.2 Solve -S_tilde * sol1_p = src_p - temp1 (RHS)
      SolverControl solver_S(maxiter, tol * temp_1.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      //Note we have already constructed S-tilde as - S_tilde 
      solver_cg.solve(negative_S_tilde, sol1_p, temp_1, preconditioner_S);

      // temp_1.reinit(dst.block(0));

      //Step 2
      // Step 2: Solve the correction system
        //         [ I   D^-1*B^T ]
        //         [ 0     alphaI  ] dst = sol1
        // =====================================================

        //2.1 scaling alpha

        dst.block(1) = sol1_p;
        dst.block(1) *= 1. / alpha; //scaling 1/alpha * sol1_p
      
        //2.2 dst(0) = sol1_u - inv(D)B.T dst(1)

        dst.block(0) = sol1_u; // Start with sol1_u.
        TrilinosWrappers::MPI::Vector tmp = src.block(0); //to have same dim as sol1_u
        B_T->vmult(tmp, dst.block(1)); ////tmp = BT*dst.block(1)
        tmp.scale(diag_D_inv); //tmp = inv(D)*tmp
        dst.block(0) -= tmp; //sol1_u - tmp
        
    }
  protected:
    const double alpha = 0.5; // parameter (0,1]

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix negative_S_tilde;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    TrilinosWrappers::MPI::Vector neg_diag_D_inv;
    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;
  };

  // Yosida preconditioner -- the inverse of Mu is replaced by the inverse of it's diagonal's elements
  class PreconditionYosida 
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::SparseMatrix &M_,
               const TrilinosWrappers::MPI::BlockVector &sol_owned)
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;
      M = &M_;
      
      diag_D_inv.reinit(sol_owned.block(0));
      neg_diag_D_inv.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D_inv.locally_owned_elements())
      {
        //Note : we have assembled M as M/deltat
        diag_D_inv[i] = ( 1.0 / M->diag_element(i));  //  dt * (Mii)^-1
        neg_diag_D_inv[i] = ( -1.0 / M->diag_element(i));  //  dt * (Mii)^-1
      }

      // Create negative_S_tilde
      B->mmult(negative_S_tilde, *B_T, neg_diag_D_inv);
    
      // Initialize the preconditioners
      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(negative_S_tilde);
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const 
    {
      const unsigned int maxiter = 100000;
      const double tol = 1e-2;

      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      // Store in temporaries the results
      TrilinosWrappers::MPI::Vector yu = src.block(0);
      TrilinosWrappers::MPI::Vector yp = src.block(1);
      TrilinosWrappers::MPI::Vector tmp = src.block(1);
      TrilinosWrappers::MPI::Vector tmp2 = src.block(0);

      //Step 1
      // Step 1.1) yu = F^-1 * src.0
      solver_gmres.solve(*F, yu, src.block(0), preconditioner_F);
      
      //Step 1.2) yp = negative_S_tilde^-1(src1-B*yu)
      B->vmult(tmp, yu); //tmp = B*yu
      tmp.add(-1.0, src.block(1)); // tmp = src.block(1) - tmp
      // neg_S*yp = (src(1) - Byu)==tmp(RHS)
      SolverControl solver_S(maxiter, tol * tmp.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      solver_cg.solve(negative_S_tilde, yp, tmp, preconditioner_S);

      //Step 2) 
      // Step 2.1) dst1 = yp
      dst.block(1) = yp; 

      // Step 2.2) dst0 = yu - F^-1*B_T*yp
        //Step 2.2.1) 
      B_T->vmult(tmp2, dst.block(1)); //tmp2 = B_T*yp (rhs)

      //Solve the linear system 
      res.reinit(src.block(0)); //to store the result of the  lin sys F res = tmp2
      dst.block(0) = yu; //init final velocity dest 
      SolverControl solver_F2(maxiter, tol * tmp2.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres2(solver_F2);
      solver_gmres2.solve(*F, res, tmp2, preconditioner_F); // res = F^-1 * tmp2 
      dst.block(0).sadd(-1,res); //update final velocity dest dstu = yu - res

    }

  protected:

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    const TrilinosWrappers::SparseMatrix *M;
    TrilinosWrappers::SparseMatrix negative_S_tilde;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    TrilinosWrappers::MPI::Vector neg_diag_D_inv;
    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;

    mutable TrilinosWrappers::MPI::Vector res;
  };

#endif
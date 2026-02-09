#include "NavierStokes.hpp"

void 
NavierStokes::setup()
{
// Create the mesh.
{
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridGenerator::channel_with_cylinder(mesh_serial,
                                        0.03, 
                                        2,    
                                        2.0,  
                                        true);
    const Point<dim> center(0.2, 0.2);
    static const SphericalManifold<dim> manifold_description(center);

    mesh_serial.set_manifold(0, manifold_description); 
    mesh_serial.refine_global(2);

    GridTools::partition_triangulation(mpi_size, mesh_serial);

    const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
        << std::endl;
    }

pcout << "-----------------------------------------------" << std::endl;

// Initialize the finite element space.
{
    pcout << "Initializing the finite element space" << std::endl;

    const FE_Q<dim> fe_scalar_velocity(degree_velocity);
    const FE_Q<dim> fe_scalar_pressure(degree_pressure);

    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                        dim,
                                        fe_scalar_pressure,
                                        1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
        << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
        << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
        << std::endl;

    quadrature = std::make_unique<QGauss<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
        << std::endl;

    quadrature_face = std::make_unique<QGauss<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per face = " << quadrature_face->size()
        << std::endl;
}

pcout << "-----------------------------------------------" << std::endl;

// Initialize the DoF handler.
{
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0]    = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1]    = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
}

pcout << "-----------------------------------------------" << std::endl;

// Initialize the linear system.
{
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
        {
        for (unsigned int d = 0; d < dim + 1; ++d)
            {
            if (c == dim && d == dim) // pressure-pressure term
                coupling[c][d] = DoFTools::none;
            else // other combinations
                coupling[c][d] = DoFTools::always;
            }
        }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
        {
        for (unsigned int d = 0; d < dim + 1; ++d)
            {
            if (c == dim && d == dim) // pressure-pressure term
                coupling[c][d] = DoFTools::always;
            else // other combinations
                coupling[c][d] = DoFTools::none;
            }
        }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
        block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);

    // Update old solution
    solution_old.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);

    solution_owned = 0.0;
    solution = 0.0;
    solution_old = 0.0;
}
}

void NavierStokes::assemble()
{

}

void
NavierStokes::solve()
{

}

void
NavierStokes::output(const unsigned int &time_step)
{

}

void NavierStokes::run()
{

}

void NavierStokes::compute_lift_drag(const double &U_ref, const double &L_ref)
{

}


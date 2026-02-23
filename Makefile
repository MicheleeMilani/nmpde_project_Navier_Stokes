# Directory where .geo files are located
MESH_DIR = mesh
GMSH = gmsh

# Default rule (prints help)
help:
    @echo "Usage:"
    @echo "  make mesh NAME=turek DIM=2          -> Generates turek.msh (2D)"
    @echo "  make mesh NAME=Cylinder3D DIM=3     -> Generates Cylinder3D.msh (3D)"
    @echo "  make mesh-all                       -> Generates ALL project meshes"
    @echo "  make clean-mesh                     -> Removes all .msh files"

# Generates a specific mesh (e.g., make mesh NAME=turek DIM=2)
mesh:
    $(GMSH) $(MESH_DIR)/$(NAME).geo -$(DIM) -format msh2 -o $(MESH_DIR)/$(NAME).msh

# Generates all meshes (useful for initial setup)
mesh-all:
    $(GMSH) $(MESH_DIR)/Cylinder2D.geo -2 -o $(MESH_DIR)/Cylinder2D.msh
    $(GMSH) $(MESH_DIR)/turek.geo -2 -o $(MESH_DIR)/turek.msh
    $(GMSH) $(MESH_DIR)/turek_lite.geo -2 -o $(MESH_DIR)/turek_lite.msh
    $(GMSH) $(MESH_DIR)/turek_super_lite.geo -2 -o $(MESH_DIR)/turek_super_lite.msh
    $(GMSH) $(MESH_DIR)/Cylinder3D.geo -3 -o $(MESH_DIR)/Cylinder3D.msh
    $(GMSH) $(MESH_DIR)/mesh-cube.geo -3 -o $(MESH_DIR)/mesh-cube.msh

# Removes binary mesh files to clean the repository
clean-mesh:
    rm -f $(MESH_DIR)/*.msh

.PHONY: help mesh mesh-all clean-mesh
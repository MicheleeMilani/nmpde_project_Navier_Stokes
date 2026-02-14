# Directory dove si trovano i file .geo
MESH_DIR = mesh
GMSH = gmsh

# Regola di default (stampa un aiuto)
help:
	@echo "Utilizzo:"
	@echo "  make mesh NAME=turek DIM=2          -> Genera turek.msh (2D)"
	@echo "  make mesh NAME=Cylinder3D DIM=3     -> Genera Cylinder3D.msh (3D)"
	@echo "  make mesh-all                       -> Genera TUTTE le mesh del progetto"
	@echo "  make clean-mesh                     -> Rimuove tutti i file .msh"

# Genera una mesh specifica (es. make mesh NAME=turek DIM=2)
mesh:
	$(GMSH) $(MESH_DIR)/$(NAME).geo -$(DIM) -format msh2 -o $(MESH_DIR)/$(NAME).msh

# Genera tutte le mesh (comodo per il primo setup)
mesh-all:
	$(GMSH) $(MESH_DIR)/Cylinder2D.geo -2 -o $(MESH_DIR)/Cylinder2D.msh
	$(GMSH) $(MESH_DIR)/turek.geo -2 -o $(MESH_DIR)/turek.msh
	$(GMSH) $(MESH_DIR)/turek_lite.geo -2 -o $(MESH_DIR)/turek_lite.msh
	$(GMSH) $(MESH_DIR)/turek_super_lite.geo -2 -o $(MESH_DIR)/turek_super_lite.msh
	$(GMSH) $(MESH_DIR)/Cylinder3D.geo -3 -o $(MESH_DIR)/Cylinder3D.msh
	$(GMSH) $(MESH_DIR)/mesh-cube.geo -3 -o $(MESH_DIR)/mesh-cube.msh

# Rimuove i file binari della mesh per pulire il repo
clean-mesh:
	rm -f $(MESH_DIR)/*.msh

.PHONY: help mesh mesh-all clean-mesh
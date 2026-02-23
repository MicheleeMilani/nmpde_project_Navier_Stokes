### Organizing the source code
Please place all your sources into the `src` folder.

Binary files must not be uploaded to the repository (including executables).

Mesh files should not be uploaded to the repository. If applicable, upload `gmsh` scripts with suitable instructions to generate the meshes (and ideally a Makefile that runs those instructions). If not applicable, consider uploading the meshes to a different file sharing service, and providing a download link as part of the building and running instructions.


### Mesh Generation
Before compiling the solver, you must generate the required mesh files. A `Makefile` is provided in the root directory to automate this process using Gmsh.

* **To generate a specific mesh** (e.g., the 3D Cylinder):
    ```bash
    $ make mesh NAME=Cylinder3D DIM=3
    ```
* **To generate all project meshes at once:**
    ```bash
    $ make mesh-all
    ```
* **To clean up generated mesh files:**
    ```bash
    $ make clean-mesh
    ```

---


### Compiling
To build the executable, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```
The executable will be created into `build`, and can be executed through
```bash
$ ./navier_stokes
$ ./navier_stokes3D
```

## Git
We have used `git` as version control and have several branches for different version of code we wrote, listed here

  - `master` : serial naive imeplementation
  - `omp` : serial naive imeplementation with openmp
  - `tree` : serial Octree implementation
  - `cu` : CUDA native implementation
  - `cu_tree` : CUDA imeplementation with flattened octree
  - `cl` : incomplete OpenCL implementation
  - `cross_prod` : old implementation with cross product

As the entirety of the project is too large to fit in a moodle upload, please check our github repo [a_bit_advanced_rt](https://github.com/hpsc-iitb/a_bit_advanced_rt) for everything.

The project uses the `cmake` build system and were run on a Linux based system with `gcc-7` and `CUDA-10.1`. 

### To build
  - make a new directory inside the root dir (out of source builds are preferred), for example `build`
  - `cd build`
  - `cmake ..`
  - `make`

Note: Branches `tree` and `cu_tree` require [SFML](https://www.sfml-dev.org/) to build and visualize. Download cmake for your system, extract, rename the folder to SFML and put the folder in the root of project dir, like

```
a_bit_advanced_rt
|- build
|- include
|- models
|- SFML
|- src
|- CMakeLists.txt
```
Alternatively, you can edit the `CMakeLists.txt` to point to the `SFML` dir. Look for the line 
```
set(SFML_DIR "${CMAKE_SOURCE_DIR}/SFML/lib/cmake/SFML")
```

### To run
  - copy one of the meshes to the build directory, and make sure the `main.cpp / main.cu` has the correct file name. Look for the line
  ```c++
    std::string file = "shadow";
  ```
  - run using `./raytrace`, the default executable name, it will either produce a `a.ppm` image (can be opened with GIMP) or display it in SFML, depending on the branch built.
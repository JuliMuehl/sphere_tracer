Basic Pathtracer Written in C supports 
diffuse, microfacet and dielectric materials, 
homogeneous media, 
NEE (Next Event Estimation) based pathspace importance sampling (for non participating media).

![out.png](md_assets/out.png)
# Dependencies
The pathtracer itself only uses the standard library and pthreads for parallelization. However for the purpose of displaying the current state of the image in a window GLFW, OpenGL and GLEW are used.
# Compiling on Linux
To build the executable run
```bash
$ make 
```
If you don't have gcc installed you need to change the first line in the Makefile to e.g. in case of clang
``` Makefile
CC=clang
```
Then you can run the Program with
```bash
./sphere_tracer
```
At the moment the scene has to be customized by changing the src/sphere_tracer.c file this will be changed in the future. 

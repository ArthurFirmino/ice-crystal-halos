# ice-crystal-halos

Repository for the code from my MSc project titled [**Physically Based Rendering of Ice Crystal Halos**](https://arthurfirmino.github.io/msc_project.pdf).

Consists of a phase function plugin for the [Mitsuba](https://www.mitsuba-renderer.org/) renderer written in C++, 
and a MATLAB script for pre-computing the tabulated phase function.

## Sample Renders
![Alt Text](https://github.com/ArthurFirmino/ice-crystal-halos/raw/master/figures/cza.png)
![Alt Text](https://github.com/ArthurFirmino/ice-crystal-halos/raw/master/figures/22halo.png)

## Usuage

### Generating the Phase Function file

  1. Compile (with optimizations) *ics_function* and *write_out* using MATLAB Coder
  2. Run *script.m* producing **.data** file
  
### Compiling the plugin
  
  1. Read Mitsuba documentation for details on how to compile Mitsuba
  2. Copy *icecrystal.cpp* to the *mitsuba/src/phase* directory
  3. Add ```plugins += env.SharedLibrary('icecrystal', ['icecrystal.cpp'])``` to the SConscript file
  4. Compile Mitsuba
  
### Rendering

  1. Follow Mitsuba documentation for creating a scene with a heterogeneous participating media
  2. Set the medium's phase function to **icecrystal**
  3. Set the string parameter **filename** under **icecrystal** to the path of the generated file
  4. Render away!

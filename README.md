# SPH High Speed Collision Simulation

This Smoothed Particle Hydrodynamics C program simulates a high speed colision of a metal cilinder in a wall. The code and all parameters are based in the Fortran program presented in the book: Liu G. R. and Liu M. B., Smoothed particle hydrodynamics: a meshfree particle method, World Scientific, 3rd printing, 2003.

# Usage

Download and extract or clone the project, then run make in the terminal to compile.

```sh
$ make
```

After the compile finishes, run the sph executable directing the output to a file data.dat for instance.

```sh
$ ./sph > data.dat
```

The simulation will take a while and when it finishes, run the load 'plot.gnu' command into GNUPLOT to generate the "frames" of the simulation in the pict folder (which has to be created first).

```sh
gnuplot> load 'plot.gnu'
```

A video of the resulting simulation can be found in Youtube:
Plastic collision: https://www.youtube.com/watch?v=kVukwV27IQQ
Elastic collision: https://www.youtube.com/watch?v=1aEDRHqGNiU

Other simulation results based in the same code
Two circular bodies central collision: https://www.youtube.com/watch?v=0Ibab7-Jmrk
Two circular bodies deslocated collision: https://www.youtube.com/watch?v=T0jQkL9EL4Q
Two circular bodies lateral collision: https://www.youtube.com/watch?v=GW9220yNQmU
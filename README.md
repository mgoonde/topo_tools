# topo_tools

set of small tools to see how they work

Installation:

 First untar and compile nauty:
 
    tar -xf nauty26r11.tar.gz

    cd nauty26r11

    ./configure

    make

    cd ..

 then compile the files in this directory

    make all


Programs currently included:

- coords_hash.x: 
  from input coordinates generate connectivity, canonical ordering, and hash.
  Run:

     ./coords_hash.x < input_coords

- site_hash.x:
  from input configuration, go through every site, generate local configuration
  with a spherical cutoff Rcut, generate connectivity and hash of this small local
  configuration. Run:

     ./site_hash.x < input_si_tetra

- distance_tool.x:
  from two configurations in input, calculate the smallest distance from configuration_1
  to configuration_2, permutations included. Run:

     ./distance_tool.x < input_distance


Configurable files: 
 
- color_cutoff.dat: 
   file containing information about cutoff for different atomic types when generating the 
   connectivity matrix.

- rcut.dat:
   file containing Rcut parameter for constructing smaller local environment

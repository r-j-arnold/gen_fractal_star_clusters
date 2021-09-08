# gen_fractal_star_clusters
Generate star forming regions with kinematic and spatial substructure.

[![DOI](https://zenodo.org/badge/399783453.svg)](https://zenodo.org/badge/latestdoi/399783453)

Becky Arnold 25/08/21 (r.j.arnold@keele.ac.uk)

Code that generates simulated star forming regions via an adapted box 
fractal method.

It outputs three files, r.txt which holds the x, y and z coordinates of the
stars in columns, v.txt which holds the vx, vy and vz stellar velocities 
in columns, and m.txt which holds the stellar masses.

To run the code:

- Must be run with python, preferably python3.
- On the command line specify the number of stars, required probability 
  of survival and desired half mass radius in parsecs via the command line.
- The probability of survival must be between 0 and 1, and is how likely
  a star is to survive each generation of the fractal. For example 0.2
  would result in a highly substructured cluster, and 1 will give a smooth
  spherical cluster.

Example of how to run the code: python3 make_box_fractal.py 1000 0.5 5

This would generate a cluster of 1000 stars that is moderatly spatially 
substructured and has a half mass radius of 5 parsecs.

This code was written using with the following package versions, though it
may work with earlier/later versions.

Python version: 3.6.8
Numpy version: 1.19.5
Matplotlib version: 3.3.4


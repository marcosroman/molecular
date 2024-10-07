molecular
=========

![simulations](./simtilevids.gif)

simple molecular dynamics simulator
(particles in a 3D box with periodic boundary conditions)
using molecular dynamics

the particles interact in a pairwise fashion via the 
**Lennard-Jones potential**
```math
V_{LJ} (r) = 4*\epsilon*((\sigma/r)^12-(\sigma/r)^6)
```
(where $r$ is the pairwise distance)

this simple potential represents an interaction
that is
	repulsive at short distances
 	and slightly attractive past a fixed treshold

```
+--------------------------------------------------------------------+
|  +      *       +      +      +       +      +      +       +      |
|         *                                                          |
|         *                                                          |
|         *   V_{LJ}(r)/epsilon                                      |
|         *                                                          |
|         *                                                          |
|         *                                                          |
|         *             *********************************************|
|         *         *****                                            |
|         *      ***                                                 |
|         *     **                                                   |
|          *   *                                                     |
|          *  *                                                      |
|          * *                                                       |
|  +      + **    +      +      +       +      +      +       +      |   r
+--------------------------------------------------------------------+ -----
0.5     1      1.5     2     2.5      3     3.5     4      4.5     5   sigma
```


[spanish version](./README_es.md)


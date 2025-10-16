In OpenMX PDOS calculations, the labels d1 to d5 represent specific real d-orbitals that the projected density of states is calculated for, following a standard convention. Specifically, these labels correspond to the following orbitals:

- d1: $$d_{3z^2 - r^2}$$ orbital
- d2: $$d_{x^2 - y^2}$$ orbital
- d3: $$d_{xy}$$ orbital
- d4: $$d_{xz}$$ orbital
- d5: $$d_{yz}$$ orbital

This mapping follows the convention detailed in OpenMX support forums and documentation, which clarifies that the ordering implemented in PDOS output files is consistent across s-, p-, and d-orbitals (for example, s1=s, p1=px, p2=py, p3=pz, and similarly d1 through d5 for the various d-type orbitals).[1][6][8]

These orbitals use real spherical harmonics, not complex ones, which is typical for most DFT codes that work with real basis sets. The assignment is based on global Cartesian coordinates, meaning the defined $$x$$, $$y$$, and $$z$$ axes are those of the simulation cell.[6]

### Orbital Order and Usage
- This information is highly relevant for analyzing the PDOS contributions from individual d components, which is critical when interpreting the electronic structure, magnetism, or chemical bonding properties of transition metal compounds and other systems with d electrons.[8][1][6]
- The same file naming applies for s and p orbitals (s1, p1, p2, p3), so if you see outputs like *.PDOS.Tetrahedron.atomX.d1 through d5, you can directly map each to the corresponding real d orbital as listed above.[1][8]

[1](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=1388)
[2](https://www.openmx-square.org/openmx_man3.9/node70.html)
[3](https://pranabdas.github.io/espresso/hands-on/pdos/)
[4](http://abacus.deepmodeling.com/en/v3.3.0/advanced/elec_properties/dos.html)
[5](https://www.openmx-square.org/openmx_man3.9/node25.html)
[6](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3440)
[7](https://www.openmx-square.org/openmx_man3.5/node45.html)
[8](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=1835)
[9](https://wiki.fysik.dtu.dk/ase/ase/calculators/openmx.html)
[10](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=1936)
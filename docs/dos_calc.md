# DOS and pDOS Calculations

This tutorial explains how to perform Density of States (DOS) and projected Density of States (pDOS) calculations using OpenMX, following the conventional scheme.

## 1. Self-Consistent Field (SCF) Calculation

First, perform an SCF calculation for your system. As an example, consider carbon diamond. Prepare an input file `Cdia.dat` in your `work` directory with the following keywords for DOS calculation:

```
Dos.fileout      on
Dos.Erange       -25.0  20.0
Dos.Kgrid        12 12 12
```

- `Dos.fileout on` enables DOS output.
- `Dos.Erange -25.0 20.0` sets the energy range for DOS (in eV), where 0.0 is the chemical potential.
- `Dos.Kgrid 12 12 12` sets the k-point grid for DOS calculation.

Run the SCF calculation with:

```
./openmx Cdia.dat
```

After completion, you will find `cdia.Dos.val` (text) and `cdia.Dos.vec` (binary) in the `work` directory. These files contain eigenvalues and eigenvectors needed for DOS analysis.

## 2. DOS and pDOS Calculation

Next, compile the DOS analysis tool. Move to the `source` directory and compile:

```
make DosMain
```

Copy the resulting `DosMain` executable to your `work` directory. Then, run:

```
./DosMain cdia.Dos.val cdia.Dos.vec
```

The program will interactively prompt you for options:

- **Method:** Choose between Tetrahedron (1) or Gaussian Broadening (2).
- **DOS or PDOS:** Select DOS (1) or PDOS (2).
- **Atoms for PDOS:** If PDOS is selected, specify atom indices (e.g., `1 2`).

For example, to calculate PDOS for atom 1 using the Tetrahedron method:

```
Which method do you use?, Tetrahedron(1), Gaussian Broadeninig(2)
1
Do you want Dos(1) or PDos(2)?
2
Number of atoms=2
Which atoms for PDOS : (1,...,2), ex 1 2
1
```

The program will generate files such as:

- `cdia.PDOS.Tetrahedron.atom1.s1`
- `cdia.PDOS.Tetrahedron.atom1.p1`
- `cdia.PDOS.Tetrahedron.atom1.p2`
- `cdia.PDOS.Tetrahedron.atom1.p3`

Each file contains:
- Column 1: Energy (eV)
- Column 2: DOS or PDOS (eV⁻¹)
- Column 3: Integrated DOS or PDOS

For spin-polarized calculations, columns 2 and 3 are for up and down spins, and columns 4 and 5 are the corresponding integrated values.

If you select the Gaussian broadening method, you will be asked to set the broadening parameter `a` (in eV), which controls the width of the Gaussian function.

## References
1. https://openmx-square.org/openmx_man3.9/node70.html
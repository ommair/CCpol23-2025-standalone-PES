# CCpol23+(2025)-standalone-PES
The model improves the description of short- and long-range three-body interactions and extends the framework introduced in the seminal work:

CCpol23+ (2025): Three-Body Interaction Energy Evaluation for Water Trimers

This repository contains the reference implementation of the CCpol23+ (2025) three-body potential energy surface (PES) for water trimers.
The model provides improved short- and long-range three-body interaction energies, extending the methodology introduced in:

Góra et al., “Predictions for water clusters from a first-principles two- and three-body force field,” J. Chem. Phys. 140, 194101 (2014)
DOI: 10.1063/1.4875097

Current version written by: Ommair Ishaque (Dec 12, 2025)

## Overview

CCpol23+ (2025) computes three-body interaction energies of water trimers using an improved FS3 potential energy surface.
The trimer configuration space is separated using the average monomer separation:

PES Component	Region Definition	Purpose
Short-range	Ravg < 5.5 Å	Exchange-repulsion and short-range induction
Long-range	Ravg > 4.5 Å	Smooth decay of long-range 3-body terms
Overlap	4.5 Å – 5.5 Å	Blended SR/LR interaction via s(x)

The final three-body energy is:

F_S3 = s * F_S3_short + (1 - s) * F_S3_long

🔧 Smoothing Function
Switching variable
x = (Ravg - Rmin) / (Rmax - Rmin)

Rmin = 4.5 Å
Rmax = 5.5 Å

Cubic switching function s(x)
s = 1                      for x < 0
s = 1 + x^2 (2x - 3)       for 0 < x < 1
s = 0                      for x > 1


This ensures continuity of the PES and smooth transitions between the short- and long-range fits.

## Project Structure
File / Folder	Description
input_geo	Input geometries of water trimers (atomic O and H sites only). Includes 13 example trimers.
fs3_modules/	Modules implementing the short-range PES, long-range PES, damping functions, and switching function.
Makefile	Build script for compiling all modules and linking the executable.
execute	Main executable (generated after running make).
README.md	Documentation file.
* **Input Requirements:** input_geo

The program reads a file named input_geo containing trimer geometries in Cartesian coordinates (Å).
Only atomic O and H positions are required.

Example format:

O x y z
H x y z
H x y z
O x y z
...


Each trimer is processed independently.

The program outputs:

U0 (two-body reference)

NB(ind) (non-bonded induction)

3B(FS3) (three-body short/long-range potential)

Total energy = U0 + NB(ind) + 3B(FS3)

All energies are given in kcal/mol.

## Compilation & Execution
Build the program
make

Run the executable
./execute

Remove build files
make clean

## Example Output
**** CCpol23+ (2025) ****
**** kcal/mol ****

    Ravg                 U0              NB(ind)            3B(FS3)          U0+NB+3B
2.7825691668     -2.9796101980     -4.2734095396     -0.5809315073     -7.8339512449
3.2355450187     -8.0284188892     -1.4547084661     -0.1192495859     -9.6023769412
3.5590989501     -6.7298255593     -0.7317844679     -0.0386210144     -7.5002310416
...

## Scientific Background

The CCpol23+ (2025) model refines the three-body representation originally introduced in:

Góra et al., Journal of Chemical Physics 140, 194101 (2014)

Enhancements include:

Updated short-range fit using Ravg < 5.5 Å trimers

Improved long-range analytic damping behavior

Smooth blending across the overlap region (4.5–5.5 Å)

More physically consistent decay of three-body interactions

Increased accuracy for large Ravg configurations relevant to liquid water simulations

## License

Specify your chosen license here, e.g.:

MIT License (recommended for scientific software)

## Author

Ommair Ishaque
Ph.D. Candidate – Computational Physics / Molecular Simulations
University of Delaware
GitHub: ommair

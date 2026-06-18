# Sebastians_GRIT
GRIT-Code ([Gr]een-Functions for [I]nteracting [T]hermodynamics)

===============

General Remarks

===============

This code is at the very end of a relatively long chain of software packages (Quantum Espresso, Wannier90, Wannier-Tools, Electron-Phonon-Wannier and Python for the post-processing of the individual steps).

One needs to compute

- Fermi-energy with a nscf calculation (QE)
- Electronic energies on a fine grid (either WT or EPW)
- Phononic energies on a fine grid (either QE or EPW)
- Electron-phonon-couplings on a fine grid (EPW)

The original EPW-code doesn't support the functionality of printing the couplings to file. Therefore, use 

https://github.com/Skalhoef/Sebastians_Personalized_QE_n_EPW_Codes.git 

with input flags "print_fine_Fermi = .true." and "prtgkk_sebbe = .true.".

There are several example-values for various materials (metals and insulators) in the code.


===================

Commands for compiling and running the code

==================


In order to clean the directory of any executables

    make -f Makefile_PC clean


In order to compile the code (after installing / loading the neccessary packages)
 
    make -f Makefile_PC all



In order to run the code

    mpirun -np <number_of_processes> ./grit.x

    e.g. 

    mpirun -np 1 ./grit.x
    mpirun -np 2 ./grit.x
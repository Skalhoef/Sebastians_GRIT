# Sebastians_GRIT
GRIT-Code ([Gr]een-Functions for [I]nteracting [T]hermodynamics)

===============

General Remarks

===============

This code is at the very end of a relatively long chain of software packages (Quantum Espresso, Wannier90, Wannier-Tools, Electron-Phonon-Wannier and Python for the post-processing of the individual steps).

It is therefore essential to read and understand the loaders and the input-files to understand how the code expects certain numbers (e.g. couplings-constants or band-numbers) of individual materials to be stored. Example-files might be added in the future (especially if requested).



===================

Useful Commands

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
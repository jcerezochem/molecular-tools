molecular-tools
===============

Basic tools for molecular modeling manipulations

Molecular Tools is a buch of Fortran90 modules that undertake a number of basic molecular modelling / computational chemistry analysis and manipulations, such as structure file reading and writing (including pdb, gro, g96, xyz), and parsing output files from programs as Gaussian (fchk and log) and Molpro retrieving informations such as geometry, frequencies, normal modes, Hessian, Gradient, molecular orbital coefficients... A set of tools developed using these modules is distributed altogether.

This set of subroutines has been developed during my PhD program, so they might be too rudimentary for advanced purposes, but they can be of some help in order to develop new small projects written in Fortran90 handling chemical information. Such tasks can also be undertaken by means of more advanced parsers such as cclib ([cclib.sourceforge.net](http://cclib.sourceforge.net)) and, in this sense, our contribution is more limited but, because of that, it might be simpler to use. Furthermore, note that this project also adds support for fchk Gaussian files (actually, it is mainly focus on analysing output from such [software](http://www.gaussian.com)). Also, another remarkable point is that all tools are developed with command-line workflow in mind, so that it is easy to use in within a script.

Enjoy!


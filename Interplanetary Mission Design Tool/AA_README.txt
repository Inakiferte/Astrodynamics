<^> INSTRUCTION FOR THE CODE <^>
By: Iñaki Fernandez, Jorge Simón and Javier Sanchez.
Master In Space And Aeronautical Engineering.
Astrodynamics Course: Interplanetary Mission Analysis
Design Tool.
<^>==========================<^>

========================== FAST EXECUTION =========================
For fast try out of the code, you can find the main_PCP.m
script which is self explanatory. This will compute the PCP
for the chosen departure days and time of flight from Earth
to another plante. The code is ready to be execute.
===================================================================

========================== DIRECTORIES ============================
In this directory you can find 5 subdirectories:

- assigment1 : Here you will find all the code from assigmnet 1,
which allows to compute the state vector for a given day. This functions
are used in state_vector_JC.m , state_vector_JC_P.m & state_vector_JC_PC.m which
are found in the main directory.

- compute_orbits: Here you will find the main_orbits.m code. This allows to compute
the orbits for a chosen departure day and time of flight. To avoid errors, run the
main_orbits.m inside this folder. It requires the functions inside assignment1. Here
THE HYBRID SOLVER IS IMPLEMENTED!!!!!

- deltaV : Here the final deltaV matrix will be stored.

- orbital_elements_planets: Here the orbital elements of all the planet are stored, which
are afterwards readed by the main_PCP.m

- PCP_Plots: Here the PCP plots will be stored.

- verification: Here you can find the code used for the verification part. We mantained 
the hole folder since all functions are used.   

=====================================================================

===================== INTERESTING FUNCTIONS & SCRIPTS ===============
We want to remark that all the functions used in the main_PCP are 
explained inside each function, remarking the inputs and the outputs.

- main_PCP.m : main script to create PCPs.

- Mesh_test.m : this script is used to compute the time study of the 
report.

- z_solver_v2.m : transfer function solver using bisection method.

=====================================================================

If there is any problem just write to: inakiphy@gmail.com or jorge.simon@estudiantat.upc.edu or javier.sanchez.rodriguez@estudiantat.upc.edu 


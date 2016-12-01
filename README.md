# Fisher_waves

Splitting step algorithm to simulate the 2d sFKPP equation


Run make to compile on UNIX systems. 
Run [./fisher_waves ly ] in terminal, where ly is an integer ( lateral system size )  



Files: 

fisher_waves.cpp -> main file ( run simulations and write down data )

random_no_generators.h -> random number routines ( gamma, gauss and poisson distribution )

simulation.h -> contains the splitting step algorithm 

splitting_routines.h -> routines for the splitting step algorithm

initial_conditions.h -> Initial conditions (only step function at the moment)

read_write_msg.h -> input and output routines
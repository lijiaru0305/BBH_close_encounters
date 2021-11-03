# BBH_close_encounters
We use the open-source N-body code Rebound to simulation the close encounter dynamics of the multiple (~2) BHs in AGN disks. 

The code in this repository can output the rate of properties of the close encounters.
BBH_close_encounters.py: the main program that runs the experiment
drag_ecc_inc_mt.so: storing the C functions for external forces
src/disk_force.c: the source code for the damping force 

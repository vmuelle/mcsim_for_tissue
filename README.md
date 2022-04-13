#Introduction

This projct helps to simulate the properties of skin tissue regarding interference with light.
For the simulation MCML.exe is used.
This executable takes a file with input parameters. 
These Parameters are derived from Moco, New insight into the origin of remote PPG signals in visible light and infrared(2018).


#Execution

To execute a simulation with mcml, run the file launch_simulation_mcml.m (This is a matlab program).
It uses the file Tissue.m to calculate the parameters given by Moco for wavelength between 400 and 700 nm and generates input files for mcml  
There will be 4 files for different skin settings. 
n/c means normal or compressed skin.
d/s means skin in diastolic or systolic state.
These files are then simulated one by one with mcml.exe and the reults will be stored in data_files/outputs.
The output files will have the same n/c d/s signature and a number corresponding to the wavelength of light the simulation was run with.

To further look into the results, run lookmcml.m
Here you also have to manually edit the code and give the name of the simulated output (.mco) to look into a specific result. 

To visualize the penetration depth and the depth origin of the simulated tissue over different wavelenght, run visPD_DO.m

To execute a simulation with mcxyz, run the file launch_simulation_mcxyzn.m (This is a matlab program).
It uses the file Tissue.m to calculate the parameters given by Moco. Instead of a cylindrical coordinate system, a cartesian coordinate system is used. 
This allows to simulate more complex structures than just layers. In this program a circular object is inserted into the layer structure. 
The simulation is only executed for normal tissue in diastolic state at wavelenght 400nm.
The results of the simulation are stored in several binary files starting with moco_params and can be used for further visualization tasks. 
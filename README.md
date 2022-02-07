#Introduction

This projct helps to simulate the properties of skin tissue regarding interference with light.
For the simulation MCML.exe is used.
This executable takes a file with input parameters. 
These Parameters are derived from Moco, New insight into the origin of remote PPG signals in visible light and infrared(2018).


#Execution

To make a file of input parameters, execute writecoefficientstofile.m (This is a matlab program).
To change parameters, edit the code of this file. The coefficients are at the top of this file. 
The files will be stored in data_files/inputs. 
There will be 4 files for different skin settings. 
n/c means normal or compressed skin.
d/s means skin in diastolic or systolic state.
These files can now be simulated one by one with mcml.exe.
When running mcml.exe the reults will be stored in data_files/outputs.
The output files will have the same n/c d/s signature and a number corresponding to the wavelength of light the simulation was run with.
To further look into the results, run lookmcml.m
Here you also have to manually edit the code and give the name of the simulated output (.mco) to look into a specific result. 
To visualize the penetration depth and the depth origin of the simulated tissue, run visPD_DO.m
## About
This short script can help you visualize ESP points from Gaussian16 logfiles and pyRESP logfiles using PyMol. 

![image](https://github.com/Olecon/Visualization-ESP-from-Gaussian-pyRESP/assets/107347740/bd6a98b4-f596-412e-89e6-ff1ac3dafb75)

## Examples
load Gaussian16 file:
```
run main.py
load_gaussian_esp $PATH_to_Gaussian_file
```
* load_gaussian_esp ./test/Gaussian.log

load PyRESP file:
```
run main.py
load_RESP_files $PATH_to_outfile $PATH_to_datfile Molecule
```
* load_RESP_files ./test/RESP_charges.dat ./test/RESP.dat 1:3

# Prerequisites
* python (3.6)
* colour 0.1.5
* numpy (1.26.4)
* cclib (1.8)

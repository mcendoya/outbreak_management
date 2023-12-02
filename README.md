# Outbreak management plans

Data, code and results for the study: An individual-based spatial epidemiological model for the spread of plant diseases.
"Performance of outbreak management plans for emerging plant diseases: the case of almond leaf scorch caused by *Xylella fastidiosa* in mainland Spain"
by Cendoya et al., 2023.

### Dataset

**coords.csv** contains the coordinates (x, y) of the almond trees in the study area. They have been generated with a row spacing of 7x7 m in the identified plots with this crop type. Also includes the 1ha grid cell in which they are located and the initial state (t0) of the disease. 

### Simulations of outbreak management plans

**model_functions.py** contains the functions used for the simulation of disease spread and the implementation of management plans including surveillance and control.

**run_simulations.py** runs the functions to generate the simulations for the baseline scenario without interventions and with each combination of monitoring and control.

### Summary results

**summary_results_scripts** folder contains the scripts used to process the results and obtain the interpretable databases.

### Data and Code citation

[![DOI](https://zenodo.org/badge/726494586.svg)](https://zenodo.org/doi/10.5281/zenodo.10251506)


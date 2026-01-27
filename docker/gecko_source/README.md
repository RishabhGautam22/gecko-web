### GECKO-A wiki
A detailed wiki is avalaible on GECKO-a gitlab at the following URL : https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/gecko-a/-/wikis/home

### Structure of the repository

##### DATA

Contains all the input data such as experimental rate constants, inorganic chemical schemes, forced reactions, etc.

##### INPUT

Contains input files that can be changed: (1) cheminput.dat is the file where you write the precursors for which you want to generate the gas phase oxidation mechanism and (2) gecko.nml is the namelist file where you can change some parameters such as the maximum number of generation, the saturation vapor pressure threshold, etc.

##### LIB

Contains all the GECKO-A code in fortran 90 files.

##### OBJ

This is the compilation directory. Type *make clean* to clean the folder and then *make* to compile the code.

##### RUN

This the directory where you run the model using gecko.sh script. Results are in OUT folder.


### Compile the code
In the OBJ folder, clean the previous compilation if needed with *make clean*, and compile the code with *make*.

### Run the code
To configure the mechanism, edit *gecko.nml* and *cheminput.dat* in INPUT. Then go in the RUN directory and run *gecko.sh*



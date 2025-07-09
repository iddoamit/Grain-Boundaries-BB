**************************************************
Introduction
**************************************************

This package contains a script that simulates a structure of three grains and two grain boundaries, with grain size variation

**************************************************
Pre-requisites
**************************************************

For the simulation to run, you will need to have:
1. A working installation of Python3, with the following packages:
	a. NumPy - install using python3 -m pip install numpy
	b. Pandas - install using python3 -m pip install pandas
	c. Subprocess - Native to Python3
	d. ltspice - install using python3 -m pip install ltspice
	e. SciPy - install using python3 -m pip install scipy
	f. os - Native to Python3
2. An installation of LTSpice - You will need to note down to yourself the location of the executable (scad3.exe) file
3. Microsoft Office - specifically Excel
4. Download the two simulation files "GBLIB.py" and "params.xlsx" from this repository 

**************************************************
Simulation parameters
**************************************************

There are two categories for the simulation parameters. 

--------------------------------------------------
Materials parameters - locked
--------------------------------------------------
These parameters are in the "params.xlsx" file, within the locked "Materials" tab. The parameters are:
1. Nominal doping (N) - 3e16 cm^-3
2. Relative permittivity (er) - 11.7
3. Surface density of trapped charges (NT) - 3.5e11 cm^-2
4. Temperature (T) - 300 K
5. Cross sectional area (A) - 4e-10 cm^2
6. Nominal grain length - 5e-5 cm

--------------------------------------------------
Simulation parameters - user controlled (unlocked)
--------------------------------------------------
These parameters are in the "params.xlsx" file, within the locked "Simulation" tab. The parameters are:
1. Grain length variation - given as a fraction of the nominal value (e.g., a value of 0.1 means 5e-6 cm). Grain sizes will be sampled from a normal distribution with the stated standard deviation
2. Number of simulation iterations - recommended 1000 or more
3. File name - the simulation netlist and output files will be saved using this name ("filename.net" and "filename.raw"). The output file will use this name too ("filename.xlsx") so it is important to change it between iterations.
4. LTSpice path - The location of the ltspice executable on your machine. It is normally set to the Mac default, but needs to be changed for other computers.

**************************************************
Running a simulation
**************************************************

To run a simulation: 
1. Copy the two files ("GBLIB.py" and "params.xlsx") to a folder of your choice
2. Edit the "Simulation" tab in "params.xlsx":
	a. Enter the path to the LTSpice executable ("scad3.exe") - This needs to be done only once
	b. Enter the grain length variation
	c. Enter the number of iterations
	d. Enter a file name for this simulation
3. Save the "params.xlsx" file. Do not change its name.
4. Open a shell / console and run the Python file (e.g., ">> python3 GBLIB.py")
5. When the simulation has ended, a message will appear on the screen with the location of the output file.

**************************************************
Simulation output
**************************************************

The simulation will output a Microsoft Excel file with two tabs:
1. Current tab - the first column is the applied voltage. The other columns are the calculated currents (in A) for each iteration
2. Resistance tab - the first column is the applied voltage. The other columns are the calculated resistance (in Ohm) for each iteration

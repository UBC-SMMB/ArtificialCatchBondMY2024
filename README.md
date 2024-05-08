# System Requirements
OS: Windows 10+

Software: MATLAB 2023b+

# Installation
1. Download the contents of this zip file.
2. Add the folder + subfolders to your working directory in MATLAB.

## Generating the figures used in the paper
1. Clear your workspace.
2. Import the contents of FigureData.mat.
3. Run the ArtificialCatchBondFigures.m script. The subsections of this script should be runnable in any order.

## Running the Analytical Simulation
1. Clear your workspace.
2. Open runCatchBondSim.m. Set parameters in the first subsection:

- number of sequences to generate and average (nSeq)
- sequence (sLength array containing the DNA lengths you would like to probe; minimum length is 7 bp)
- salt (struct containing sodium, potassium, and magnesium concentrations in M)
  
4. Run subsection 1. Depending on how many sequences/parameters you are testing, this could take several minutes.
5. Run subsection 2. This will create an intersections array containing all sequence pairs which have the jaw and hook unzipping in the correct order.

## Running the Monte Carlo Simulation
1. If you have just run the analytical simulation, skip to step 5.
2. Clear your workspace.
3. Import the contents of IntersectionsData.mat.
4. Open runCatchBondSim.m.
5. Scroll down to the heading "3. Monte Carlo Simulation to get catch bond rupture times"

There are 3 modes: 
- non-linear force ramp (this is what was used in the paper). It simulates our specific optical tweezers trap stiffness and our specific dsDNA handle length, which leads to non-linear force ramps.
- linear force ramp.
- force clamp.

Set the initial parameters (commented) and run the section which corresponds to the mode you wish to use.


   

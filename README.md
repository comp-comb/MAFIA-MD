# MAFIA-MD: Molecular Arrangements and Fringe Identification and Analysis from Molecular Dynamics
To Find out the number of cyclic rings from molecular dynamics simulation. 

# Directory structure
  a) parent directory

      - mafiamd.py: the main python code
      - Makefile: Makefile for Linux-like systems
      - requirements.yml: requirement file to build the conda environment
      - LICENSE: terms and condition of the license for this code
      - README.md: this "readme" file
      - callgraph.png: the callgraph of the code

  b) external_tool : this directory contains the external tool xyz2mol from https://github.com/jensengroup/xyz2mol

      -xyz2mol.py: xyz2mol source code
      -LICENSE: terms and condition of the license for xyz2mol

  c) input : this directory contains the example input files
  
      -set1_validation_demo : example input trajectory files for ring detection and chemical analysis segment of MAFIA-MD

          -fabricated.xyz 
          -fabricated2.xyz
          -real_MD.xyz

      -set2_fringe_analysis_demo: example input trajectory file for fringe analysis segment of MAFIA-MD

          -1769atoms.xyz
  d) output : a blank directory for saving the output files from MAFIA-MD

  e) ancillary_script : contains a bash script for splitting a large trajectory file containing multiple timesteps into a number of trajectory files containg one timestep each. 

        - splitting.sh : the bash script
        - splitting_demo.xyz: a large trajectory file containg multiple timesteps. 

# Requirements
rdkit is required for the chemical representations. Conda is required to install rdkit. 

# Installation guide
1. install conda : (large installation)

        https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

    or Miniconda: (minimal installation)

        https://docs.conda.io/en/latest/miniconda.html#

2. create a virtual environment using the Makefile:

    a) Linux or MacOS systems: 

        make
    b) for Windows systems, use the Anaconda prompt installed during anaconda or miniconda installation and execute the following in the desired directory:
    
        conda env create -f requirements.yml -p MAFIAMD

3. activate the virtual environment:

        conda activate ./MAFIAMD

4. run the code: 

        python mafiamd.py
  
5. to deactivate the virtual environment: 
        
        conda deactivate

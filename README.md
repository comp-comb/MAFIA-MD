# MAFIA-MD: Molecules, Aromatics and Fringe Identification and Analysis from Molecular Dynamics
To Find out the number of aromatic rings from molecular dynamics simulation. 

# Requirements
rdkit is required for the chemical representations. Conda is required to install rdkit. 

# WorkFlow
1. install conda : (large installation)

  * https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

  or Miniconda: (minimal installation)
  https://docs.conda.io/en/latest/miniconda.html#

2. create a virtual environment using the Makefile:

  * make

3. activate the virtual environment:

  * conda activate ./MAFIAMD

4. run the code: 

  * python mafiamd.py
  
5. to deactivate the virtual environment: 
  * conda deactivate

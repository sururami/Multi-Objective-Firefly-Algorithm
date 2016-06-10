This is the Linux free software distribution of the I-PAES algorithm proposed in the following papers (please cite both articles):


V. Cutello, G. Narzisi, G. Nicosia,A Multi-Objective Evolutionary Approach to the Protein Structure Prediction Problem ,Journal of the Royal Society Interface, Royal Society Publications London, 3(6):139-151, 2006.

V. Cutello, G. Narzisi, G. Nicosia,A Class of Pareto Archived Evolution Strategy Algorithms using Immune inspired Operators for Ab-Initio Protein Structure Prediction,Third European Workshop on Evolutionary Computation and Bioinformatics, EvoWorkshops 2005, EvoBio 2005,30 March - April 2005, Lausanne, Switzerland.Springer, LNCS 3449:54-63, 2005.



The distribution contains three sub-directories:
 
- src/        #source code directory
- bin/        #TINKER executables directory
- instances/  #protein intances directory


1) Compilation and installation:

To compile the software type the command "make" inside the directory "src".

I-PAES code uses some external routines fron the TINKER Molecular Modeling Package: 

- analyze
- protein
- xyzpdb

and the force field parameter set of CHARMM (version 27) energy
function:

-charmm27.prm

You can download these files directly from the TINKER web-site 
(http://dasher.wustl.edu/tinker/).

Routine analyze, protein and xyzpdb must be saved inside the directory "bin".
Force field parameter file (charmm27.prm) must be saved inside directory "bin/params".

NOTE: to compile the software you need installed on your machine the gcc compiler for C.

2) Ouput data files:

The program generates output files of the results inside 
a new directory with name given in input
The directory cointains:
1.parameters.dat
*         This file prints the values of all the parameters given 
*         in input to the algorithm

2.archive_statistics.dat
*        This file maintains statistical information  about the archive a 
*        for all the iterations (archive length, mean objective values, 
*        minimal energy, ecc).

3.current_solution_statistics.dat
*        This file maintains information about the total energy and the 
*        values of the two objectives for the current solution.

4.final archive directory        
*        The algorithm creates a new directory with the informations about the 
*        final archive found. For each solution a Protein Data Bank (PDB) 
*        file is also created.

3) Input data files: 

Only one input file is needed by the algorithm.
The file must contain the amonacid sequence of the 
protein (in standard format: PHE, ASN, MET, ecc) and the predicted 
secondary or supersecondary structure for each residue.
One line for residue is needed, an example follows:

PHE U
ASN U
MET H
GLN H
CYS H
GLN E
ARG E
ARG H
PHE H

4) Test example:

The software is distributed with a test input files 
for protein 1ZDD (inside instances directory).

To run the test simulation on 1ZDD use the bash-script file 
"simulation.sh" from the root directory of the software.

Complite description about the input parameters of the algorithm can be found in the "i-paes.c"
source file in the "src" directory.



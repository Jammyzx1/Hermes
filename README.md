# Hermes
Structural analysis tool for normal modes distortion and molecular dynamics trajectories. Solvent shell cutting from MD trajectory. Feature reordering for solvated system feature reordering for FERBUS.

README for Hermes Version 1.2

Documentation last updated 22/09/16

Hermes is designed in a modular format for easy upkeep and modification. Hermes was principally an analysis tool when it was first written, allowing statistical data based on the sampling and kriging models to be easily calculated. Additional statistics
can easily be added if they are thought to be useful. Hermes is now also capable of:

- Statistical analysis of molecular geometries over DL POLY trajectory or TYCHE distortion.
- Cutting out the first solvation for a solute or atom of a solute from DL POLY trajectory 
- Reordering and reducing of solvated systems features for Kriging with FEREBUS  

Hermes requires an input file for which examples and explanations are given below. Each task has additional 
required input file(s), these are listed in the task descriptions in this file. Hermes has an additional 
directory of pre and post processing scripts written in python. These can prepare the other input files where 
necessary. There is also an R sub-directory, which contains an automatic plotting script for R for use after Task 1.

We provide a list of and description of each *.f90 routine below, this should be updated between versions. 
Main files, which all begin main*, are the control files these call the appropriate routines to carry out the tasks.
Module files, which all begin module*, contain subroutines, functions and constants used in the program.

The Makefile assumes the use of ifort please change this to the compiler you have (gfortran, f95, pg95 .....)

_________________________________________________________________________________________________

MAINS
_____

main-Hermes.f90
---------------
This file is main control file it requests the input file from the user (usually named input.hermes) and reads the
options. Having read the options the correct task is selected and appropriate arguments read before calling the
main subroutine related to that task


main-structural-analysis.f90
----------------------------
This file is the main control file for the structural analysis over a DL POLY trajectory or a distortion set from
Tyche. It requires an additional file explained below under Task 1 which is generated using a python script provided (file labelled.xyz).
This subroutine analyses the changes in bond length(BL), bond angle(BA) and dihedral angles(DA) over an MD trajectory or a Tyche
distortion set. The full set of possible bonds, angles and diheahrals is found by a call to the subroutine bonding 
(main-structural-analysis.f90; line 71). This list is stored in a file structure_input.txt and can be edited by the user if only a subset of interest 
are required. To do this choose to exit not carry on at the question Do you wish to continue or do you want to look at the input file the that has been automatically created?

The routines to calculated BL, BA and DA are in module-subroutines.f90 and module-functions.f90

Six output files are generated two for each BL, BA and DA:
BONDS.csv: a full output of the BL for each bond specified in structure_input.txt followed by summary statistics on each bond 
(note statistics like RMSD are calculated relative to the first example)
BONDS_summary.csv: Contains only the statistics for each BL in the file structure_input.txt
ANGLES.csv: a full output of the BA for each bond specified in structure_input.txt followed by summary statistics on each bond 
(note statistics like RMSD are calculated relative to the first example)
ANGLES_summary.csv: Contains only the statistics for each BA in the file structure_input.txt
DIHEADRAL.csv: a full output of the DA for each bond specified in structure_input.txt followed by summary statistics on each bond 
(note statistics like RMSD are calculated relative to the first example)
DIHEDRAL_summary.csv: Contains only the statistics for each DA in the file structure_input.txt

Additionally, three directories are created one for BL (BONDS), BA(ANGLES) and DA(DIHEDRALS), these contain the information needed to use the
automatic plotting scripts from the R directory. Place the R script in the parent directory of these directories and run the R script as explained below.

Please note this set of routines WILL OVERWRITE EXISTING FILES


main-rdf-f90
-------
This file is the control file for cutting out the first solvation shell of a whole solute molecule. The first solvation shell is defined as any water solvent oxygen within 
the maximum first solvation shell of the constituent atoms. i.e. we create and radial distribution function (rdf) for each atom in the molecule and take the maximum distance from the 
rdf origin to the first solvation shell peak over all of the constituent atoms to be the cutoff for the first solvation shell. 

The code contains a limited library of atom types. These are defined from line 69 to 102 in version 1.1. If your molecule contains new atom types not in this list YOU NEED TO ADD THEM. 
If you add any please do the following 
Line 69; increase this to the total number of atom type definitions there will be once you have added you new definitions
Line 103; begin adding your new definitions to the array atom_types the second index should be 1 for the atom type label ('CH3')
Line 104; Linked to line 103 add your new definitions length to the array atom_types the second index should be 2 for the atom type label length enter as a character string ('2')

The main output file is called HISTORYCUT.xyz this is a stacked xyz file of the solute and the first solvation shell water solvent. Additional outputs are given
rdf; a file of the rdf data for plotting
history.xyz; an xyz file of the original histroy file


main-microsolvation.f90
-----------------------
This file is the control file for cutting out the first solvation shell of one atom of a solute. The first solvation shell is defined using an rdf, which has it origin on an atom in the solute. 
The first solvation shell is then distance between the origin and the first peak in the rdf. This is a specific use for investigating hydration effects local to one atom in a solute. This code
has seen the least development since its original version as it is very specific. 

On line 289 is an if statement all atom types in the solute must be listed here, if your solute atom type is not listed here please add it. This is an inefficient and ugly loop which will be the focus
of development if work in this direction continues.


main-multipole-prediction-stats.f90
-----------------------------------
This file is the only file required for to calculated statistics based on a kriging model for predictions using AIMAll multipoles. 
The file gives an output of of prediction_evaluation_metrics.txt


main-iqa-prediction-stats.f90
-----------------------------------
This file is the only file required for to calculated statistics based on a kriging model for predictions using AIMAll IQA terms. 
The file gives an output of of prediction_evaluation_metrics.txt

main-densnodes_auto.f90
-----------------------------------
This file caluclates the solvent density in a particular region as a weighted sum accounting for denisty in neighbouring areas. The resolution parameter determines the size of the volume elements and the overlap determine the amount of influence neightbouring areas denisty has. Nodes are placed at the regions of highest denisty.

main-waterorderauto.f90
-----------------------------------
Here the code orders the water features to represent them in relation to the node positions. This reduces the conformational space of the water as for example water 1 will always the closest water molecule to node 1 dispite the fact that the particular moelcule labelled as water one originally may now be the cloest water molecule to node 7. 

main-anova_auto.f90
-----------------------------------
Code performs an ANOVA to define which features are most important to describing that atom

main-noderank_auto.f90
-----------------------------------
Now the features are reordered and reduced based on importance.

main-solorder.f90
-----------------------------------
Reorders the solvent molecules.

MODULES
_______

module-bonding.f90
------------------
This code is quick but not fool proof, the code determines the molecule structure i.e. what is bonded together and hence what are valid angles and dihedrals. The criteria are 
for heavy atoms (not H atoms) are considered bonded if the distance is > 0.4A but < 2.0A. If it is a  H then it is bonded if the distance is > 0.4A but < 1.2A (A = Angstroms).
The code is quite messy but fairly efficient. The process is done step wise and built up from bonds to determine valid angles and dihedrals. Initially the input and output files are opened.
Then the code determines how many lines to read of a stacked xyz file to cover one geometry. The code then calculates the separating distances between a given atom and all other atoms in the solute.
This is done for all atoms one at a time in the solute. An array of this information is stored and printed in a log file. The array is sorted by the number of bonds and to insist that H is placed at the top
as this will only ever have one bond in reality this is done via a bubble sort algorithm. The code then runs through the ordered list and determines whether valid bonds are made between pairs of atoms.
valid bonds are then printed to the output file used in main-structural-analysis.f90 for task 1. The same process but over three and four atoms respectively is then carried out to find valid angles and dihedrals.

The out file is structure_input.txt and has the following format

For a bond starts with a capital B then lists the two atom names in the bond
B H1 C1

For an angle starts with capital A then lists the three atom names in the angle
A H1 C1 C2

For a dihedral angle starts with capital D then lists the four atom names in the dihedral
D H1 C1 C2 C3


module-subroutines.f90
----------------------
This module contains all subroutines used generically in various main codes


module-functions.f90
--------------------
This module contains all functions used generically in the various main codes and subroutines


module-constants.f90
--------------------
This module contains various constant parameters used by some functions and subroutines

__________________________________________________________________________________________________________________


The Python scripts are :
-------------------------------------------------
(see directory Python-processing-scripts)

Prepare_HISTORY_for_Hermes.py    - This should be used before task 1 is run if you wish to analyses a DL POLY run. 
                                   It should be run in same directory as the DL POLY HISTORY file and will produce
                                   a labelled.xyz file which is the additional input for task 1 in Hermes.

Prepare_Tyche_gjfs_for_Hermes.py - This should be used before task 1 is run if you wish to analyses a set of Tyche  
                                   distortions. It should be run in same directory as the GJF s from Tyche and 
                                   will produce a labelled.xyz file which is the additional input for task 1 in
                                   Hermes.

STACKED-XYZ-TO-GJF.py            - Will convert a stacked xyz file from task 2 or task 3 in to a set of gjfs

STACKED-XYZ-TO-GJF-SORT.py       - Will convert a stacked xyz file from task 2 or task 3 in to a set of gjfs and
                                   sort them in to directories dependent on the number of water solvent molecules given

The R scripts are :
-------------------------------------------------
(see directory R_scripts_plots)

MDplot.R  - This script will automatically plot the data from task 1 in R. The commands to run this
            script in R are ;

1. set path to the directory where you have run Hermes <br />
`path <- 'path/to/Hermes'` <br />
`setwd(path)`
        
2. load the script <br />
`source("MDplot.R")`
	
3. run the script <br /> 
`MDplot.R(path)`
 

-------------------------------------------------

## FORMAT
###### task is always the first parameter word and defines the calculation type
###### tasks can be defined by a word or the task number
######Task 1 word structure number 1
######Task 2 word microsolv number 2
###### Task 3 word solvshell number 3
###### Task 4 word statskrig number 4
###### The file always has a last line which states `END` for the last line of the input file

------------------------------------------------

## FOR TASK 1 
###### Options are :
###### filename - list of bonds angles and dihedrals to calculate or auto to automatically generate this input
###### atom - an atom to centre the coordinates on for clearer viewing in jmol

REQUIRES A labelled.xyz FILE GENERATED FROM A DL POLY HISTORY FILE USING Prepare_HISTORY_for_Hermes.py FROM Python-processing-scripts SUB-DIRECTORY OR
 A labelled.xyz FILE GENERATED FROM TYCHE GJF FILES USING Prepare_Tyche_gjfs_for_Hermes.py FROM Python-processing-scripts SUB-DIRECTORY.
If creating a labelled.xyz file from Tyche's gjfs run the python script in the same directory as the gjfs then copy the labelled.xyz file to where you are running Hermes.
An R script is provided which can produce plots of  all data on bonds, bond angles and dihedrals deviation automatically. Just run the R script in the directory where you have
 run Hermes and provide the path to that directory as the argument to the function. See readme in the R script directory for R specific commands.

------------------------------------------------

task structure <br />
filename auto <br />
atom C1 <br />
END



task 1 <br />
filename auto <br />
atom C1 <br />
END



task 1 <br />
filename inputbonds.txt <br />
atom N1 <br />
END

------------------------------------------------
## FOR TASK 2
###### Options are :
###### Otype - The oxygen label of the solvent water molecule
###### solute - The atom label of the solute to calculate the rdf from
###### solvent charge sites - The number of solvent charge sites eg TIP3P = 3 TIP5P = 5
###### solute number - The number of atoms in the solute

 REQUIRES HISTORY FILE FROM DL POLY IN THE SAME DIRECTORY

------------------------------------------------

task microsolv <br />
Otype OTP <br />
solute name OC <br />
solvent charge sites 3 <br />
solute number 15 <br />
END



task 2 <br />
Otype OTP <br />
solute name OC <br />
solvent charge sites 3 <br />
solute number 15 <br />
END



task microsolv <br />
Otype OP3 <br />
solute name CH <br />
solvent charge sites 5 <br />
solute number 15 <br />
END

------------------------------------------------

## FOR TASK 3
###### Options are :
###### Otype - The oxygen label of the solvent water molecule
###### solute filename - The filename of a file which contains a list of solute atom types to calculate the rdf from 
###### ( it will calculate an rdf from each and take the maximum value as the first solvation shell cut off for the whole molecule)
######solvent charge sites - The number of solvent charge sites eg TIP3P = 3 TIP5P = 5
######solute number - The number of atoms in the solute

 REQUIRES HISTORY FILE FROM DL POLY IN THE SAME DIRECTORY

------------------------------------------------

task solvshell <br />
Otype OTP <br />
solute filename solute <br />
solvent charge sites 3 <br />
solute number 15<br />
END


task 3 <br />
Otype OTP <br />
solute filename solute_atoms.inp <br />
solvent charge sites 3 <br />
solute number 15 <br />
END


task 3 <br />
Otype OP3 <br />
solute filename sol.in <br />
solvent charge sites 5 <br />
solute number 15 <br />
END

------------------------------------------------

## FOR TASK 4
######Options are :
###### multipole - To analysing data after kriging multipoles
###### IQA - To analysing data after kriging IQA self and interaction terms

 REQUIRES FINPUT FILE AND PREDICTION IN THE SAME DIRECTORY

------------------------------------------------

task statskrig <br />
multipole  <br />
END <br />

task statskrig <br />
IQA <br />
END <br />

task 4 <br />
IQA <br />
END

-------------------------------------------------

## FOR TASK 5
######Options are :
###### Training_set_name - Name of the training set file
###### No_Solute_atoms - The number of solute atoms, it is assumed the solute atoms are always the first defined
###### No_waters - Number of water molecule
###### Resolution - Used to define the volume the density is calculate for. Suggested default 0.05
###### Overlap - Number of neighbouring denisty regions to consider in weighted sum for denisty definition
###### Feat_var_thresh - Treshold to consider a feature important. Suggested default 4
###### No_bins - Number of bins. Suggested default 6

 REQUIRES FINPUT FILE AND TRAINING SET FILE IN THE SAME DIRECTORY EXPECTS TO BE RUN IN A FEREBUS ATOM FOLDER

------------------------------------------------

task featorder <br />
Training_set_name       C11_TRAINING_SET.txt <br />
No_Solute_atoms         13 <br />
No_waters               26 <br />
Resolution              0.05 <br />
Overlap                 4 <br />
Feat_var_thresh         0.1 <br />
No_bins                 6 <br />
END <br />

task 5 <br />
Training_set_name       C11_TRAINING_SET.txt <br />
No_Solute_atoms         13 <br />
No_waters               26 <br />
Resolution              0.05 <br />
Overlap                 4 <br />
Feat_var_thresh         0.1 <br />
No_bins                 6 <br />
END <br />

task featorder <br />
Training_set_name       O7_TRAINING_SET.txt <br />
No_Solute_atoms         13 <br />
No_waters               26 <br />
Resolution              0.1 <br />
Overlap                 2 <br />
Feat_var_thresh         0.1 <br />
No_bins                 6 <br />
END <br />

-------------------------------------------------

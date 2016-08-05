#!/usr/bin/python

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group      ! 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import os 
from re import *
import glob 

def labeller(element,cO,cH,cC,cN) :
        if element == 'C' :
                cC = cC + 1
                element = element + str(cC)
        elif element == 'H' :
                cH = cH + 1
                element = element + str(cH)
        elif element == 'O' :
                cO = cO + 1
                element = element + str(cO)
        elif element == 'N' :
                cN = cN + 1
                element = element + str(cN)
        else :
                print('unknown element passed to labeller')
        return(element,cO,cH,cC,cN)

# User input and variables 
atom_number = input('Please enter the number of atoms ')
startline = input('Please enter the line the number the molecular coordinates begin on in the gjf files from tyche ') # The line number to start reading from gjfs 
end_line = startline + atom_number                                   # The last line in the gjf to read
#file_list = []
geom = 0
Tyche_message = 'No cell size structure from Tyche'
product = open('labelled.xyz','w')

# Write the opening lines of the labelled.xyz
product.write('gjf from Tyche\n')
product.write('Atom Count = ' + str(atom_number) + '\n')
product.write(Tyche_message + '\n')

# Write the distorted molecules postions 
file_list = next(os.walk(os.getcwd()))[2]

for file in sorted(glob.glob('*.gjf')) :
	print file
	line_number = 0
	cO = 0
	cH = 0
	cC = 0
	cN = 0
	atcount = 0
	geom = geom + 1
	with open(file,'r') as tychegjf :
		for line in tychegjf :
			line_number = line_number + 1
			if line_number == startline - 1 :
				product.write('\n----------- ' + str(geom) + ' ----------\n')
			elif line_number >= startline and line_number < end_line :
				words = line.split()
				element = words[0]
				lab_element = labeller(element,cO,cH,cC,cN)
				element = lab_element[0]
				cO = lab_element[1]
				cH = lab_element[2]
				cC = lab_element[3]
				cN = lab_element[4]
				atcount = atcount + 1
				if atcount > atom_number :
					print('WARNING the atom count of structure ', geom, ' is great than the atom count input ', atom_number)
				line = element + '\t' + words[1] + '\t' + words[2] + '\t' + words[3] + '\n'
				product.write(line)
			elif line_number == end_line + 1 :
				product.write(Tyche_message)
	tychegjf.close()
product.write('\nEND\n')
product.close()


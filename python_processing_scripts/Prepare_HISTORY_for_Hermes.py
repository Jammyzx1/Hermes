#!/usr/bin/python

# To calculate bond length, angle and dihedral angle distortion from a set of stacked xyz files.

from re import *
import math, os, subprocess

# Functions
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
		print('unknown element passed to labeller', element)
	return(element,cO,cH,cC,cN)

# Step 1 make a labelled xyz file

#------- Get the user input -------

boxsize = str(raw_input('Please, enter the length of one of the box edges this code assumes the box is cubic - '))
watO = str(raw_input('Please enter the water oxygen atom type (eg OTP) if any or NNAA if not- '))
watH = str(raw_input('Please enter the water hydrogen atom type (eg H3P) if any or NNAA if not- '))
name = 'labelled.xyz'

#------- Set the regex -------

pattern='timestep' 
#pattern='        0.0000000000        0.0000000000       ' + boxsize
#pattern='        0.0000000000        0.0000000000       50.0000000000'
print('regular expression pattern will be - ',pattern)
#watO = 'OT'
#watH = 'H3P'

#------- Open the files -------

filein = open('HISTORY','r')
fileout = open( name,'w')
temporyfile = open('tempory.out','w')
print('Opened files, ', filein.name, fileout.name)

i = 0
j = 0	
store = 0
countatom = 0
atcount = 0
atomcount = 0
cO = 0
cH = 0
cC = 0
cN = 0
for line in filein :
	i = i + 1
	if  i == 1 :
		title = line
		fileout.write(title)

	startline = match(pattern,line)

	if startline :
		j = j + 1
		fileout.write('Cubic cell of size ' + boxsize)
		fileout.write('\n----------- ' + str(j) + ' ----------\n')
		# element counts initalisation below
		cO = 0
		cH = 0
		cC = 0
		cN = 0
		atomcount = atcount
		atcount = 0

	valid = match('[CHNOTW1]',line) # Add to this list any other valid values in the atom names you want included in the analysis
	waterO = match(watO,line)
	waterH = match(watH,line)

	if waterO or waterH :
		countatom = countatom + 1
	else :
		if valid :
			element = line[0]
			store = i + 1
			lab_element = labeller(element,cO,cH,cC,cN)
			element = lab_element[0]
			cO = lab_element[1]
			cH = lab_element[2]
			cC = lab_element[3]
			cN = lab_element[4]
			atcount = atcount + 1
		else :
			if i == store :
				scientific_notation = line.find('E')
				if (scientific_notation != -1) :
					coor = line.split()	
					Efoundx = coor[0].find('E')
					Efoundy = coor[1].find('E')
					Efoundz = coor[2].find('E')
					if (Efoundx != -1) :
						x = ('% 12.10f' % float(coor[0]))
					else :
						x = coor[0]
					if (Efoundy != -1) :
						y = ('% 12.10f' % float(coor[1]))
					else :
						y = coor[1]
					if (Efoundz != -1) :
						z = ('% 12.10f' % float(coor[2]))
					else :
						z = coor[2]
					#coords = '    '.join(map(str,(
					fileout.write(element + '    ' + x + '       ' + y + '       ' + z + '\n') 
					if j == 1 :
						temporyfile.write(element + '    ' + x + '       ' + y + '       ' + z + '\n')
				else :
					coord = line
					fileout.write(element + coord)
					if j == 1 :
						temporyfile.write(element + coord)
fileout.write('END')
filein.close()
fileout.close()
temporyfile.close()
stratomcount = str(atomcount)
command = "sed -i '2i Atom Count : %s' labelled.xyz" % (stratomcount) # Mac's might have a problem with this line
subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

#------- display tempory file -------
temporyfile = open('tempory.out','r')
print('A new file has been created called ',fileout.name)
print('Below is a copy of the unique atom identifiers that have been assigned.\n')
print('You need to provide an input file using these atom assignments\n\n')
print('For bond lengths the format is : B atom identifier 1 atom identifer 2\n\n')
print('For angles the format is : A atom idenifier 1 atom identifier 2 atom identifer 3\n\n')
print('For dihedreal angles the format is : D atom identifer 1 atom identifer 2 atom idenifier 3 atom identifier 4\n\n')
print('The last line of the file must always end with the word END')
for line in temporyfile :
	print(line)

temporyfile.close()
os.remove('tempory.out')
# END OF PREPROCESSING LOOP

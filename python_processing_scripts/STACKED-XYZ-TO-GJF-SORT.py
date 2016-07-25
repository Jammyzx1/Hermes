#!/usr/bin/python

from re import *
import string
import array
import os
import shutil
import glob

pattern = 'Timestep'
line_counter = 0
system_counter = 0
level_of_theory = '# SP B3LYP/6-311+G** NOSYMM OUTPUT=wfn 6d 10f PUNCH=DERIVATIVE'
atom_count = 0
atomsno = 0
solventno = 0
startline = 0
endline = 0
write_to_file = 0
list_size = -1 
solvent_list = []

filein = str(raw_input('Please enter the file name containing the stacked XYZ file from MD-to-XYZ - '))
molecule = str(raw_input('Please enter the molecules name - '))
atomsno = int(raw_input('Please enter the number of solute atoms - '))
solventno = int(raw_input('Please enter the number of solvent atoms per solvent molecule - '))

input=open(filein, 'r')
for line in input :
	line_counter += 1
#	print(line_counter) # debug

	if write_to_file == 1 :
		if line_counter <= endline :
			outline = ' ' + line
			output.write(outline)
		else :
			output.write('\n\n')
			output.write(fileout + '.wfn\n\n')
			write_to_file = 0
			output.close()

	if line_counter == endline + 1 :
		atom_count = int(string.strip(line))
                solvent_water_count = (atom_count - atomsno)/solventno
		unique = solvent_list.count(solvent_water_count)
		if unique == 0 :
			solvent_list.append(solvent_water_count)
			list_size += 1

#               print(solvent_water_count) # debug

	
	startline = match(pattern,line)
	if startline :
		startline = line_counter
		endline = startline + atom_count
#		print('in startline sc = 0', startline, endline) # debug
		system_counter += 1
		fileout = molecule + '-' + str(system_counter) + '-' + str(solvent_water_count) + '.gjf'
		output=open(fileout,'w')
#		output.write('%nproc = 12\n')
		output.write(level_of_theory)
		output.write(' \n\n')
		output.write(fileout)
		output.write(' \n\n')
		output.write(' 0 1\n')
		write_to_file = 1
				
input.close()

current = os.getcwd()
print('The number of solvents molecules are ; ', solvent_list)

for item in solvent_list :
	num_solvent = item
#	print('item ; ',item, ' num_solvent ; ', num_solvent) # debug
#	print(solvent_list.index(2)) # debug
	dirname = str(num_solvent) + '_solvent_molecules'
	try:
		os.stat(dirname)
	except:
		os.mkdir(dirname)

	extname = '*-' + str(num_solvent) + '.gjf'
	files = glob.iglob(os.path.join(current, extname))
	for file in files:
		if os.path.isfile(file):
			shutil.copy2(file, dirname)

allfiles = glob.iglob(os.path.join(current, '*.gjf'))
try:
	os.stat('All_gjfs')
except:
	os.mkdir('All_gjfs')

for file in allfiles :
	shutil.move(file, 'All_gjfs')

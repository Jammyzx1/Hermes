#!/usr/bin/python

from re import *
import string

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

filein = str(raw_input('Please enter the file name containing the stacked XYZ file from MD-to-XYZ - '))
molecule = str(raw_input('Plase enter the molecules name - '))
atomsno = int(raw_input('Please enter the number of solute atoms - '))
solventno = int(raw_input('Please enter the number of solvent atoms per solvent molecule - '))

input=open(filein, 'r')
for line in input :
	line_counter += 1
#	print(line_counter) # debug

	if write_to_file == 1 :
		if line_counter <= endline :
			output.write(line)
		else :
			output.write('\n\n')
			output.write(fileout + '.wfn\n\n')
			write_to_file = 0
			output.close()

	if line_counter == endline + 1 :
		atom_count = int(string.strip(line))
                solvent_water_count = (atom_count - atomsno)/solventno
#               print(solvent_water_count) # debug

	
	startline = match(pattern,line)
	if startline :
		startline = line_counter
		endline = startline + atom_count
#		print('in startline sc = 0', startline, endline) # debug
		system_counter += 1
		fileout = molecule + '-' + str(system_counter) + '-' + str(solvent_water_count) + '.gjf'
		output=open(fileout,'w')
		output.write('%nproc = 12\n')
		output.write(level_of_theory)
		output.write(' \n\n')
		output.write(fileout)
		output.write(' \n\n')
		output.write('0 1\n')
		write_to_file = 1
				
input.close()

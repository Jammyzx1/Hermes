import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot(title, labx, laby, xdata, ydata, fname):
	"""
	A plotting function for the bond, angles and dihedral data. 
	"""

	# delete the column title from the lists
	del xdata[0]
	del ydata[0] 
	
	# map data from strings to float
	xdata = map(float, xdata)
	ydata = map(float, ydata)
	
	# plot scatter
	fnam = fname + ".png"
	fig = plt.figure(figsize=(20, 10))
	plt.subplot(1, 2, 1)
	plt.scatter(xdata, ydata, s=80, facecolors='none', edgecolors='k')
	
	# histogram
	plt.subplot(1, 2, 2)
	plt.hist(ydata, bins='auto', histtype="bar", color="k", orientation="horizontal")

	plt.savefig(fnam)
	plt.close()

def bonds(lfin):
	"""
	Get the bond data from HERMES
	"""

	# loop over all files in the csv files in the current directory 
	for ent in lfin:
		print "Current file being analysed {} .....".format(ent)
		title = ent.split(".")[0].strip() + "Bond Length"
		labx = "Time Steps (1000's)"
		laby = "Bond Length (A)"
		xdata = [elt.split(",")[3].strip() for elt in open(ent, "r")] 
		ydata = [elt.split(",")[-1].strip() for elt in open(ent, "r")] 
		plot(title, labx, laby, xdata, ydata, ent.split(".")[0].strip())


def angles(lfin):
	"""
	Get the angle data from HERMES
	"""

	# loop over all files in the csv files in the current directory 
	for ent in lfin:
		print "Current file being analysed {} .....".format(ent)
		title = ent.split(".")[0].strip() + "Bond Angle"
		labx = "Time Steps (1000's)"
		laby = "Bond Angle (Degrees)"
		xdata = [elt.split(",")[4].strip() for elt in open(ent, "r")] 
		ydata = [elt.split(",")[-1].strip() for elt in open(ent, "r")] 
		plot(title, labx, laby, xdata, ydata, ent.split(".")[0].strip())

def dihedrals(lfin):
	"""
	Get the diheadral data from HERMES
	"""

	# loop over all files in the csv files in the current directory 
	for ent in lfin:
		print "Current file being analysed {} .....".format(ent)
		title = ent.split(".")[0].strip() + "Dihedral Angle"
		labx = "Time Steps (1000's)"
		laby = "Dihedral Angle (degrees)"
		xdata = [elt.split(",")[5].strip() for elt in open(ent, "r")] 
		ydata = [elt.split(",")[-1].strip() for elt in open(ent, "r")] 
		plot(title, labx, laby, xdata, ydata, ent.split(".")[0].strip())

def main():
	"""
	Function to run bond, angle and dihedral angle analyses.
	"""

	count = 1
	
	# Get user input the path
	path = raw_input(str("Please enter the path to running directory of HERMES "))

	# loop over bonds, angles and dihedrals
	for index in range(1, 4, 1):
		if index == 1:
			pathbonds = os.path.join(path, "BONDS_plot")
			print "\nPath to bonds analysis; {}".format(pathbonds)
			os.chdir(pathbonds)
			files = [f for f in os.listdir('.') if os.path.isfile(f) and os.path.splitext(f)[1] == ".csv"]
			bonds(files)
		elif index == 2:
			pathangles = os.path.join(path, "ANGLES_plot")
			print "Path to angles analysis; {}".format(pathangles)
			os.chdir(pathangles)
			files = [f for f in os.listdir('.') if os.path.isfile(f) and os.path.splitext(f)[1] == ".csv"]
			angles(files)
		elif index == 3:
			pathd = os.path.join(path, "DIHEDRAL_plot")
			print "Path to the dihedrals analysis; {}".format(pathd)
			os.chdir(pathd)
			files = [f for f in os.listdir('.') if os.path.isfile(f) and os.path.splitext(f)[1] == ".csv"]
			dihedrals(files)

	print "FINISHED"

if __name__ == "__main__":
	main()	

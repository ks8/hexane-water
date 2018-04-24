"Script to perform post-production density profile calculation, by Kirk Swanson"
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys, argparse
from math import *
import numpy as np
import time
import copy, os
from ast import literal_eval as make_tuple
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import math
import csv
import pandas as pd 
import scipy
from scipy.optimize import curve_fit


def indices_containing_substring(the_list, substring):
    indices = []
    for i, s in enumerate(the_list):
        if substring in s:
            indices.append(i)
    return indices

# Upload the density profile data
trajectory = open('interface-dens.xyz', 'r')
contents_traj = trajectory.readlines()
trajectory.close()




# Empty files of orientations
output = open('orientations_hexane.txt', 'w')
output.close()

# Loop through the production configurations that we have (200 of them)
for k in range(200):

	print(k)


	# First timestep 
	index_start = indices_containing_substring(contents_traj, "24600")[300 + k]

	bounds_xlo = contents_traj[index_start + 1].split()[2]
	bounds_xlo = float(bounds_xlo[1:(len(bounds_xlo) - 1)])

	bounds_ylo = contents_traj[index_start + 1].split()[3]
	bounds_ylo = float(bounds_ylo[0:(len(bounds_ylo) - 1)])

	bounds_zlo = contents_traj[index_start + 1].split()[4]
	bounds_zlo = float(bounds_zlo[0:(len(bounds_zlo) - 1)])

	bounds_xhi = contents_traj[index_start + 1].split()[6]
	bounds_xhi = float(bounds_xhi[1:(len(bounds_xhi) - 1)])

	bounds_yhi = contents_traj[index_start + 1].split()[7]
	bounds_yhi = float(bounds_yhi[0:(len(bounds_yhi) - 1)])

	bounds_zhi = contents_traj[index_start + 1].split()[8]
	bounds_zhi = float(bounds_zhi[0:(len(bounds_zhi) - 1)])


	z_len = bounds_zhi - bounds_zlo
	y_len = bounds_yhi - bounds_ylo
	x_len = bounds_xhi - bounds_xlo

	volume = y_len*x_len*(z_len/200.0)

	orientation_segments_hexane = []

	# Loop through the histogram bins
	for i in range(200):

		print(i)


		total_Szz = 0.0
		vector_num = 0.0

		for j in range(500):

			# Compute Szz four times
			r1 = np.array([float(contents_traj[index_start + 2 + 14600 + 20*j].split()[1]), float(contents_traj[index_start + 2 + 14600 + 20*j].split()[2]), float(contents_traj[index_start + 2 + 14600 + 20*j].split()[3])])
			r2 = np.array([float(contents_traj[index_start + 2 + 14600 + 20*j + 1].split()[1]), float(contents_traj[index_start + 2 + 14600 + 20*j + 1].split()[2]), float(contents_traj[index_start + 2 + 14600 + 20*j + 1].split()[3])])
			r3 = np.array([float(contents_traj[index_start + 2 + 14600 + 20*j + 7].split()[1]), float(contents_traj[index_start + 2 + 14600 + 20*j + 7].split()[2]), float(contents_traj[index_start + 2 + 14600 + 20*j + 7].split()[3])])
			r4 = np.array([float(contents_traj[index_start + 2 + 14600 + 20*j + 9].split()[1]), float(contents_traj[index_start + 2 + 14600 + 20*j + 9].split()[2]), float(contents_traj[index_start + 2 + 14600 + 20*j + 9].split()[3])])
			r5 = np.array([float(contents_traj[index_start + 2 + 14600 + 20*j + 13].split()[1]), float(contents_traj[index_start + 2 + 14600 + 20*j + 13].split()[2]), float(contents_traj[index_start + 2 + 14600 + 20*j + 13].split()[3])])
			r6 = np.array([float(contents_traj[index_start + 2 + 14600 + 20*j + 15].split()[1]), float(contents_traj[index_start + 2 + 14600 + 20*j + 15].split()[2]), float(contents_traj[index_start + 2 + 14600 + 20*j + 15].split()[3])])

			if r1[2] < (bounds_zlo + (i + 1)*z_len/200.0) and r1[2] > (bounds_zlo + i*z_len/200.0):

				vec = np.subtract(r1, r3)
				theta = math.acos((np.dot(vec, np.array([0, 0, 1])))/(np.linalg.norm(vec)))
				Szz = (1.0/2.0)*(3*pow(math.cos(theta), 2) - 1)
				total_Szz += Szz
				vector_num += 1.0

			if r2[2] < (bounds_zlo + (i + 1)*z_len/200.0) and r2[2] > (bounds_zlo + i*z_len/200.0):

				vec = np.subtract(r2, r4)
				theta = math.acos((np.dot(vec, np.array([0, 0, 1])))/(np.linalg.norm(vec)))
				Szz = (1.0/2.0)*(3*pow(math.cos(theta), 2) - 1)
				total_Szz += Szz
				vector_num += 1.0

			if r3[2] < (bounds_zlo + (i + 1)*z_len/200.0) and r3[2] > (bounds_zlo + i*z_len/200.0):

				vec = np.subtract(r3, r5)
				theta = math.acos((np.dot(vec, np.array([0, 0, 1])))/(np.linalg.norm(vec)))
				Szz = (1.0/2.0)*(3*pow(math.cos(theta), 2) - 1)
				total_Szz += Szz
				vector_num += 1.0

			if r4[2] < (bounds_zlo + (i + 1)*z_len/200.0) and r4[2] > (bounds_zlo + i*z_len/200.0):

				vec = np.subtract(r4, r6)
				theta = math.acos((np.dot(vec, np.array([0, 0, 1])))/(np.linalg.norm(vec)))
				Szz = (1.0/2.0)*(3*pow(math.cos(theta), 2) - 1)
				total_Szz += Szz
				vector_num += 1.0


		if vector_num > 0.0:
			orientation_segments_hexane.append(total_Szz/vector_num)

		else:
			orientation_segments_hexane.append(-1000.0)



	output = open('orientations_hexane.txt', 'a')
	for i in range(len(orientation_segments_hexane)):
		if i == (len(orientation_segments_hexane) - 1):
			output.write(str(orientation_segments_hexane[i]))
		else:
			output.write(str(orientation_segments_hexane[i]) + '\t')
	output.write('\n')
	output.close()





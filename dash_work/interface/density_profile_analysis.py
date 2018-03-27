""" File for plotting average density profiles of interfacial systems, by Kirk Swanson"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys, argparse
from math import *
import numpy as np
import time
import copy, os
import json
import random
import csv
import pandas as pd 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Upload the density profile data
density_dat = pd.read_csv('densities.txt', sep='\t', error_bad_lines=False, header=None)
density_dat = np.asarray(density_dat)

# Compute the average for each bin
average_profile = np.mean(density_dat, axis=0)

# Plot the data
plt.plot(average_profile)
plt.savefig('test.png')



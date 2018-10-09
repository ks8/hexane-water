"""Script for analysis of Standard Error and block size of q-TIP4P/F simulation output data, by Kirk Swanson"""
# Import modules
import sys
import argparse
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
from scipy import stats 

# Functions to parse arguments
def create_parser():

	# Create parser and add arguments
	parser = argparse.ArgumentParser(description='Read data files')
	parser.add_argument('-DASH_file', dest='DASH_file', default=None, help='CSV file containing DASH data')
	parser.add_argument('-OpenMM_file', dest='OpenMM_file', default=None, help='CSV file containing OpenMM data')
	parser.add_argument('-SE_analysis', dest='SE_analysis', default=True, help='Boolean for performing SE analysis')
	parser.add_argument('-index_start', dest='index_start', default=0, help='Starting index for data analysis')
	parser.add_argument('-index_end', dest='index_end', default=1000, help='Starting index for data analysis')
	parser.add_argument('-analysis_col', dest='analysis_col', default=2, help='Data column for analysis')
	parser.add_argument('-NUM_blocks', dest='NUM_blocks', default=5, help='Nnumber of blocks for SE calculation')
	parser.add_argument('-filename', dest='filename', default='test', help='Name for output files')

	return parser

# Function to convert arguments
def convert_args(args):

	# Files dictionary
	files={}
	files['DASH_file'] = args.DASH_file
	files['OpenMM_file'] = args.OpenMM_file

	# Options dictionary
	options = {}
	options['SE_analysis'] = args.SE_analysis
	options['index_start'] = args.index_start
	options['index_end'] = args.index_end
	options['analysis_col'] = args.analysis_col
	options['NUM_blocks'] = args.NUM_blocks
	options['filename'] = args.filename

	return files, options

# Function to convert a string to a boolean for the argparse options
def str2bool(string):
    if string.lower() == 'true':
        return True
    elif string.lower() == 'false':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Function to process the data file
def process_datafile(files, options):

	# Load files and options
	DASH_file = str(files['DASH_file'])
	OpenMM_file = str(files['OpenMM_file'])
	SE_analysis = str2bool(options['SE_analysis'])
	index_start = int(options['index_start'])
	index_end = int(options['index_end'])
	analysis_col = int(options['analysis_col'])
	NUM_blocks = int(options['NUM_blocks'])
	filename = str(options['filename'])

	# Load data file
	DASH_data = np.loadtxt(DASH_file, delimiter=',')
	OpenMM_data = np.loadtxt(OpenMM_file, delimiter=',')

	# Leave out the first 100 ps as equilibration
	DASH_data = DASH_data[index_start:index_end, :]
	OpenMM_data = OpenMM_data[index_start:index_end, :]

	# Compute total number of samples
	num_samples = len(DASH_data[:, 0])

	if SE_analysis:
		# Compute SE values for different block sizes and plot
		SE_values = []
		for b in range(num_samples/2):

			blocksize = b + 1
			print(blocksize)

			num_SE_blocks = int(float(num_samples)/float(blocksize))
			DASH_averages = list(np.zeros(num_SE_blocks))
			OpenMM_averages = list(np.zeros(num_SE_blocks))

			for i in range(num_SE_blocks):

				DASH_averages[i] = np.mean(DASH_data[i*blocksize:(i + 1)*blocksize, analysis_col])
				OpenMM_averages[i] = np.mean(OpenMM_data[i*blocksize:(i + 1)*blocksize, analysis_col])

			SE_values.append([b, stats.sem(DASH_averages), stats.sem(OpenMM_averages)])

		SE_values = np.asarray(SE_values)
		plt.subplot(2, 1, 1)
		plt.plot(SE_values[:, 0], SE_values[:, 1], 'bo', label='Block Averaging', markersize=2)
		plt.ylabel('DASH Density SE')
		plt.title('Block Averaging Analysis')
		plt.subplot(2, 1, 2)
		plt.plot(SE_values[:, 0], SE_values[:, 2], 'bo', label='Block Averaging', markersize=2)
		plt.ylabel('OpenMM Density SE')
		plt.subplots_adjust(left=0.2)
		plt.xlabel('Block size')
		plt.savefig('block_averaging_'+filename+'.png')

	# Analyze data with a specified block size
	num_SE_blocks = NUM_blocks 
	DASH_averages = list(np.zeros(num_SE_blocks))
	OpenMM_averages = list(np.zeros(num_SE_blocks))

	for i in range(num_SE_blocks):

		DASH_averages[i] = np.mean(DASH_data[i*int(num_samples/num_SE_blocks):(i+1)*int(num_samples/num_SE_blocks), analysis_col])
		OpenMM_averages[i] = np.mean(OpenMM_data[i*int(num_samples/num_SE_blocks):(i+1)*int(num_samples/num_SE_blocks), analysis_col])

	print('DASH_mean: '+str(np.mean(DASH_averages)))
	print('DASH_SE (2*): '+str(2.0*stats.sem(DASH_averages)))
	print('OpenMM_mean: '+str(np.mean(OpenMM_averages)))
	print('OpenMM_SE (2*): '+str(2.0*stats.sem(OpenMM_averages)))


# Main function
def main():

	# Process arguments
	parser = create_parser()
	args = parser.parse_args()
	files, options = convert_args(args)

	# Analyze data
	process_datafile(files, options)

# Run the program
main()

	# # Discover the right block size
	# SE_values = []

	# for b in range(2000):

	# 	print(b+1)

	# 	num_SE_blocks = float(num_samples)/float(b+1)
	# 	num_SE_blocks = int(num_SE_blocks)

	# 	# Split the data into the scpecified number of blocks and compute average profiles
	# 	average_profiles_water = list(np.zeros(num_SE_blocks))
	# 	average_profiles_hexane = list(np.zeros(num_SE_blocks))

	# 	for i in range(num_SE_blocks):

	# 		average_profiles_water[i] = np.mean(density_dat_water[i*int(num_samples/num_SE_blocks):int((i+1)*num_samples/num_SE_blocks),:], axis=0)
	# 		average_profiles_hexane[i] = np.mean(density_dat_hexane[i*int(num_samples/num_SE_blocks):int((i+1)*num_samples/num_SE_blocks),:], axis=0)

	# 	# Compute bulk densities and width values for each data block
	# 	bulk_water_densities = list(np.zeros(num_SE_blocks))
	# 	bulk_hexane_densities = list(np.zeros(num_SE_blocks))
	# 	intrinsic_widths = list(np.zeros(num_SE_blocks))
	# 	thermal_widths = list(np.zeros(num_SE_blocks))
	# 	h_water = list(np.zeros(num_SE_blocks))
	# 	h_hexane = list(np.zeros(num_SE_blocks))

	# 	for i in range(num_SE_blocks):

	# 		par_init = np.array([popt_water[0], popt_hexane[0], popt_water[1]])
	# 		fit_test, cov_test = leastsq(my_residuals, par_init, args=(xdata[bulk_water_zhi_bin:bulk_hexane_zlo_bin], average_profiles_water[i][bulk_water_zhi_bin:bulk_hexane_zlo_bin], average_profiles_hexane[i][bulk_water_zhi_bin:bulk_hexane_zlo_bin], np.mean(average_profiles_water[i][bulk_water_zlo_bin:bulk_water_zhi_bin]), np.mean(average_profiles_hexane[i][bulk_hexane_zlo_bin:bulk_hexane_zhi_bin])))

	# 		bulk_water_densities[i] = np.mean(average_profiles_water[i][bulk_water_zlo_bin:bulk_water_zhi_bin])
	# 		bulk_hexane_densities[i] = np.mean(average_profiles_hexane[i][bulk_hexane_zlo_bin:bulk_hexane_zhi_bin])
	# 		intrinsic_widths[i] = box_len_av*(fit_test[1] - fit_test[0])
	# 		thermal_widths[i] = box_len_av*fit_test[2]
	# 		h_water[i] = fit_test[0]
	# 		h_hexane[i] = fit_test[1]

	# 	SE_values.append(list([b, stats.sem(bulk_water_densities), stats.sem(bulk_hexane_densities), stats.sem(intrinsic_widths), stats.sem(thermal_widths)]))

	# SE_values = np.asarray(SE_values)
	# plt.subplot(2, 1, 1)
	# plt.plot(SE_values[:, 0], SE_values[:, 1], 'bo', label='Block Averaging', markersize=2)
	# plt.ylabel('Water Density SE (g/${cm}^3$)')
	# plt.title('Block Averaging Analysis')
	# plt.subplot(2, 1, 2)
	# plt.plot(SE_values[:, 0], SE_values[:, 2], 'bo', label='Block Averaging', markersize=2)
	# plt.ylabel('Hexane Density SE (g/${cm}^3$)')
	# plt.subplots_adjust(left=0.2)
	# plt.xlabel('Block size')
	# plt.savefig('density_profile_block_averaging_1_'+filename+'.png')

	# plt.clf()
	# plt.cla()

	# plt.subplot(2, 1, 1)
	# plt.plot(SE_values[:, 0], SE_values[:, 3], 'bo', label='Block Averaging', markersize=2)
	# plt.ylabel('Intrinsic Width SE ($\AA$)')
	# plt.title('Block Averaging Analysis')
	# plt.subplot(2, 1, 2)
	# plt.plot(SE_values[:, 0], SE_values[:, 4], 'bo', label='Block Averaging', markersize=2)
	# plt.ylabel('Thermal Width SE ($\AA$)')
	# plt.xlabel('Block size')
	# plt.savefig('density_profile_block_averaging_2_'+filename+'.png')
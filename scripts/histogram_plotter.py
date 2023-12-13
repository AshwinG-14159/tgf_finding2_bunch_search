import numpy as np
import matplotlib.pyplot as plt


binning = 0.005

my_type = 1
ind_quad_thres = 300
num_quads = 4

path_results_csv = 'results_csv'
list_of_tgfs = [107.2690991,75.16806856,97.28008224,53.85038408,157.2401033,160.0727875,106.9133546,101.7770566,113.2917237,86.54987548,89.15999654,90.15747339,62.72428582,62.05877933,117.1723738,17.41574538,91.13095645,70.70809884,135.2766278,54.57677413,92.0578866,56.33161923,57.95203473,131.5604911,119.4041283,42.97344527,78.10797909,118.788337,120.1986017,41.49320457]

file = f'{path_results_csv}/angle_output_bunch_{my_type}_{ind_quad_thres}_{num_quads}_{binning}.csv'
data = np.genfromtxt(file, delimiter=',')
plt.hist(data[:,0], color = 'r', bins=50, label = f'bunch cross {ind_quad_thres}')
plt.hist(list_of_tgfs, color = 'b', bins = 50, label = 'eyeballed list')
plt.title('Histogram of tgf detections')
plt.xlabel('Earth Theta')
plt.ylabel('number of detections')
plt.legend()

plt.savefig(f'plots/histogram/histogram_{my_type}_{ind_quad_thres}_{num_quads}_{binning}.png')




from run_experiment import createFolder, run_experiment
import subprocess as sp
import pandas as pd
import numpy as np
import pickle
from matplotlib import pyplot as plt
import os

path = '/Users/bjornjuel/projects/Renzo_AA/mabe_2sensor/'
os.chdir(path)

experiment_path = 'Experiments/kyoko/'
cmd = [experiment_path + 'experiment_kyoko.cfg',
       experiment_path + 'experiment_record_activity.cfg',
       experiment_path + 'experiment_record_jory_genome.cfg']
datafile = path + 'Experiments/kyoko/version1_LOD_data.csv'
genomefile = path + 'Experiments/kyoko/version1_LOD_organisms.csv'
activityfile = path + 'Experiments/kyoko/markov_IO_map.csv'
TPMjoryfile = path + 'Experiments/kyoko/new_snapshot_data_0.csv'
data = []
genome = []
activity = []
TPM_jory = []

runs = 5

for r in list(range(0,runs)):
	print(['run number ' + str(r)])
	# Running experiment
	run_experiment(cmd[0])
	# Getting data from file
	data.append(pd.read_csv(datafile))
	genome.append(pd.read_csv(genomefile))
	# Running recorders
	run_experiment(cmd[1])
	run_experiment(cmd[2])
	print(['success!!!'])
	# Getting recordings from file
	activity.append(pd.read_csv(activityfile))
	# TPM_jory.append(pd.read_csv(TPMjoryfile))

	with open('Experiments/kyoko/version1_LOD_data.pkl', 'wb') as f:
	    pickle.dump(data, f)

	with open('Experiments/kyoko/version1_genome.pkl', 'wb') as f:
	    pickle.dump(genome, f)

	with open('Experiments/kyoko/version1_activity.pkl', 'wb') as f:
	    pickle.dump(activity, f)

	# with open('Experiments/deterministic/TPM_jory.pkl', 'wb') as f:
	#     pickle.dump(TPM_jory, f)

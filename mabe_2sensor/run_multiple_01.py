from run_experiment import createFolder, run_experiment
import subprocess as sp
import pandas as pd
import numpy as np
import pickle
from matplotlib import pyplot as plt
import os

os.chdir('/Users/Renzo/Documents/mabe_tpm')

cmd = ['Experiments/deterministic/experiment_deterministic.cfg',
'Experiments/deterministic/experiment_record_genome.cfg',
'Experiments/deterministic/experiment_record_activity.cfg']
datafile = '/Users/Renzo/Documents/mabe_tpm/Experiments/deterministic/deterministic_LOD_data.csv'
genomefile = '/Users/Renzo/Documents/mabe_tpm/Experiments/deterministic/deterministic_LOD_organisms.csv'
activityfile = '/Users/Renzo/Documents/mabe_tpm/Experiments/deterministic/markov_IO_map.csv'
TPMjoryfile = '/Users/Renzo/Documents/mabe_tpm/Experiments/deterministic/new_snapshot_data_0.csv'
data = []
genome = []
activity = []
TPM_jory = []

runs = 30

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
	# Getting recordings from file
	activity.append(pd.read_csv(activityfile))
	# TPM_jory.append(pd.read_csv(TPMjoryfile))

	with open('Experiments/deterministic/LOD_data.pkl', 'wb') as f:
	    pickle.dump(data, f)

	with open('Experiments/deterministic/genome.pkl', 'wb') as f:
	    pickle.dump(genome, f)

	with open('Experiments/deterministic/activity.pkl', 'wb') as f:
	    pickle.dump(activity, f)

	# with open('Experiments/deterministic/TPM_jory.pkl', 'wb') as f:
	#     pickle.dump(TPM_jory, f)

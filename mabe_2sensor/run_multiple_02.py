from run_experiment import createFolder, run_experiment
import subprocess as sp
import pandas as pd
import numpy as np
import pickle
from matplotlib import pyplot as plt
import os

os.chdir('/Users/bjornjuel/projects/Renzo_AA/mabe_2sensor/')

cmd = ['Experiments/deterministic_2sensor_task3/experiment_2sensor_task3.cfg',
'Experiments/deterministic_2sensor_task3/experiment_record_activity.cfg',
'Experiments/deterministic_2sensor_task3/experiment_record_jory_genome.cfg']
datafile = '/Users/bjornjuel/projects/Renzo_AA/mabe_2sensor/Experiments/deterministic_2sensor_task3/LOD_data.csv'
genomefile = '/Users/bjornjuel/projects/Renzo_AA/mabe_2sensor/Experiments/deterministic_2sensor_task3/LOD_organisms.csv'
activityfile = '/Users/bjornjuel/projects/Renzo_AA/mabe_2sensor/Experiments/deterministic_2sensor_task3/markov_IO_map.csv'
TPMjoryfile = '/Users/bjornjuel/projects/Renzo_AA/mabe_2sensor/Experiments/deterministic_2sensor_task3/jory_snapshot_data_0.csv'
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
	TPM_jory.append(pd.read_csv(TPMjoryfile, error_bad_lines=False))

	with open('Experiments/deterministic_2sensor_task3/LOD_data.pkl', 'wb') as f:
	    pickle.dump(data, f)

	with open('Experiments/deterministic_2sensor_task3/genome.pkl', 'wb') as f:
	    pickle.dump(genome, f)

	with open('Experiments/deterministic_2sensor_task3/activity.pkl', 'wb') as f:
	    pickle.dump(activity, f)

	with open('Experiments/deterministic_2sensor_task3/TPM_jory.pkl', 'wb') as f:
	    pickle.dump(TPM_jory, f)

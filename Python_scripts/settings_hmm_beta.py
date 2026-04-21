
"""
===========
Settings file
===========

Configuration parameters for beta event characterization.
"""
# Task
task = "restEO" #leftCKCevoked or restEO

# Session
sessions = ["01", "02"]
session = '01'

# Model parameters:
lfreq = 13
hfreq = 30
sfreq = 200

n_labels = 450 # left: 226, right: 224, both: 450
lag = 8
group_id = "group_model_1"
ch_type = 'source-ave'
side = "both"

# State mapping
state_mapping = {'low': 'state3',
                 'high': 'state1',
                 'gam': 'state2',
                 'bl': 'state4'}

# Task related parameters
task_parameters = {'leftCKCevoked':{
                    'n_of_tp': 701, # Number of timepoints in the epoch
                    'n_of_bl_tp': 100, # Baseline length in timepoints
                    'epo_tmin': -0.5, # Epoch start
                    'epo_tmax':3,
                    'stim_channel':['STI001']}} # Epoch end

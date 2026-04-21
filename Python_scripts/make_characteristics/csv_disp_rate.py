"""
Make CSV file of the dispersion and rate.
CSV has (subjects as rows and parcels values in the columns)
"""
import numpy as np
import scipy
import pandas as pd
import os
import scipy
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from functions_make_characteristics import cut_ones_into_list
from config import fname
from settings_hmm_beta import (state_mapping, n_labels, group_id, sessions, sfreq, task)

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Sorted mapping
sorted_map = sorted(state_mapping.items(), key=lambda x: x[1])

# Make columns for the csv
columns = ["subject","session","state_id"]
for sens in range(0,n_labels):
    columns.append('sens' + str(sens+1))

# make empty dataframe
save_dict_rate = pd.DataFrame(columns=columns)
save_dict_dispersion = pd.DataFrame(columns=columns)

# Loop over sessions
for j, session in enumerate(sessions):
    # Loop over the subjects
    for i, row in df_subjects.iterrows():
        subject = row["subject"]
        print(subject)

        hmm_dual_path = fname.hmm_dual(
            subject = subject,
            ses = session,
            job_id = group_id,
            task = task)

        # If the dual path exists, go through the sensors
        if os.path.exists(hmm_dual_path):
            # Load the dual model
            hmm_dual_dict = scipy.io.loadmat(hmm_dual_path)
            hmm_subject = hmm_dual_dict["output"]

            # Loop then over the states and parcels
            for s_key_i, s_key in enumerate(sorted_map):
                dict_append_rate = {'subject': subject,'session':session,'state_id':s_key[0]}
                dict_append_dispersion = {'subject': subject,'session':session,'state_id':s_key[0]}

                for sens_i in range(0,450):
                    sens_id = "sens" + str(sens_i + 1)
                    vmap_sens = hmm_dual_dict["output"][sens_id][0][0]["viterbi"][0][0]

                    
                    vmap_multiply = (vmap_sens == (s_key_i+1)).astype(int)
                    vmap_multiply = np.squeeze(vmap_multiply)[:]

                    # separate zeros and ones from the vmap_multiply
                    ones, zeros = cut_ones_into_list(vmap_multiply)

                    # Calculate the waiting times
                    wt = np.array([len(zeros[i]) for i in range(len(zeros))])

                    # Calculate dispersion
                    dispersion = np.mean(wt) / np.std(wt)
                    dict_append_dispersion[sens_id] = dispersion

                    # Calculate rate
                    rate = len(ones) / (len(vmap_multiply) / sfreq)
                    dict_append_rate[sens_id] = rate

                save_dict_rate.loc[len(save_dict_rate)] = dict_append_rate
                save_dict_dispersion.loc[len(save_dict_dispersion)] = dict_append_dispersion

save_dict_rate.to_csv(fname.feature_csv(feature = 'rate',job_id=group_id, task=task))
save_dict_dispersion.to_csv(fname.feature_csv( feature = 'dispersion',job_id=group_id, task=task))


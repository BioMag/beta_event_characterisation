"""
Make CSV file of the fractional occupancies.
CSV has (subjects as rows and parcels values in the columns)
"""
import pandas as pd
import scipy
import os
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from config import fname
from settings_hmm_beta import (state_mapping, n_labels, task, group_id, sessions)

# Read the subjects
df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Sorted mapping
sorted_map = sorted(state_mapping.items(), key=lambda x: x[1])

# Make columns for the csv
columns = ["subject","session","state_id"]
for sens in range(0,n_labels):
    columns.append('sens' + str(sens+1))

# make empty dataframe
save_dict = pd.DataFrame(columns=columns)

for j, session in enumerate(sessions):
    for i, row in df_subjects.iterrows():
        subject = row["subject"]
        print(subject)

        # Path for the fractional occupancy
        FO_path = fname.FO_dual(
            subject = subject,
            ses = session,
            job_id = group_id,
            task = task)
        
        if os.path.exists(FO_path):
            # Load the FO map
            dual_FO = scipy.io.loadmat( FO_path )

            # Loop then over the states and parcels
            for s_key_i, s_key in enumerate(sorted_map):
                dict_append = {'subject': subject,'session':session,'state_id':s_key[0]}
                for sensor in range(0, n_labels):
                    sens_id = "sens" + str(sensor + 1)
                    FOs_sensor_val = dual_FO["FOs_subject"][sensor][0][0][s_key_i]
                    dict_append[sens_id] = FOs_sensor_val
                save_dict.loc[len(save_dict)] = dict_append

save_dict.to_csv(fname.feature_csv(feature = 'FO',job_id=group_id, task=task))

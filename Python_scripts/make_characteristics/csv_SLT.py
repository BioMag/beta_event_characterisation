"""
Make CSV file of the state life times.
CSV has (subjects as rows and parcels values in the columns)
"""
import pandas as pd
import scipy
import os
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.append(parent_dir)

from config import fname
from settings_hmm_beta import (state_mapping, n_labels, group_id, task, sessions)

df_subjects = pd.read_csv("subject_text_files/list_of_subjects.txt", names=["subject"])

# Define basics
measures = ["mean", "median", "std", "max", "robust_max"]

# Sorted mapping
sorted_map = sorted(state_mapping.items(), key=lambda x: x[1])

# Make columns for the csv
columns = ["subject","session","state_id","measure"]
for sens in range(0,n_labels):
    columns.append('sens' + str(sens+1))

# make empty dataframe
save_dict = pd.DataFrame(columns=columns)

# Loop over sessions
for j, session in enumerate(sessions):
    # Loop over the subjects
    for i, row in df_subjects.iterrows():
        subject = row["subject"]
        print(subject)

        SLT_path = fname.SLT_dual(
            subject = subject,
            ses = session,
            job_id = group_id,
            task = task)
        
        if os.path.exists(SLT_path): 
            
            # Load the state lifetimes
            dual_SLT = scipy.io.loadmat( SLT_path )
            SLTs_subject = dual_SLT["output_SLT"]

            for s_key_i, s_key in enumerate(sorted_map):
                for m_i, measure in enumerate(measures):
                    dict_append = {'subject': subject,'session':session,'state_id':s_key[0],'measure':measure}

                    for sens_i in range(0,n_labels):
                        sens_id = "sens" + str(sens_i + 1)

                        # Take the measure (calculated with Matlab)
                        SLT_sensor_val = SLTs_subject[sens_id][0][0][s_key[1]][0][0][measure][0][0]
                        if SLT_sensor_val.shape == (1,1):
                            dict_append[sens_id] = SLT_sensor_val[0][0]
                        else:
                            dict_append[sens_id] = 0
                    save_dict.loc[len(save_dict)] = dict_append

save_dict.to_csv(fname.feature_csv(feature = 'SLT',job_id=group_id, task=task))



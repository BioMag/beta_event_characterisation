# Beta event dynamics

Git repository for HMM based beta event characterization. Organization and content of the codes are very similar with repository: https://github.com/BioMag/beta_event_dynamics

Structure of this repository is as follows:
```
beta_event_characterization
│
└─── Python_scripts
│  └─── A_manuscript_results -> This consists all plots and statistical testings
│    └─── intra_class_correlation -> Codes to estimate and visualize the ICC for Event amplitude, FO, rate and SLT-
│      └─── ICC_AE.py
│      └─── ICC_FO.py
│      └─── ICC_rate.py
│      └─── ICC_SLT.py
│    └─── multisite_stability -> Testing the differences between the two cohorts (from different sites). Contains also the statistcs and visualization.
│      └─── ttest_map_AE.py
│      └─── ttest_map_FO.py
│      └─── ttest_map_rate.py
│      └─── ttest_map_SLT.py
│    └─── plot_spatial_characteristics -> Plot how the states FO and amplitude distribute across the brain
│      └─── plot_event_amplitude_bands.py
│      └─── plot_event_amplitude_states.py
│      └─── plot_FO_states.py
│    └─── spatial_sup_reb_difference_task.py -> Visualize how the states and freuency bands modulate during task
│      └─── plot_evoked_gamma_spatial.py
│      └─── plot_induces_response_spatial.py
│
└─── make_characteristics  -> Make csv of characteristics during rest, required for spatial mapping
│  └─── csv_event_amplitudes.py
│  └─── csv_FO.py
│  └─── csv_disp_rate.py
│  └─── csv_SLT.py
│  └─── csv_commands.sh
│
└─── subject_text_files
│  └─── list_of_subjects.txt
│  └─── list_of_subjects_site2.txt
│
└─── settings_hmm_beta.py
└─── process_raw_data -> Preprocessing of the data and source transformation
│  └─── ica_annotations
│    └─── make_annotations.py
│    └─── make_ica.py
│  └─── 1_covariance.py
│  └─── 02_forward.py
│  └─── 03_inverse_matfile_rest.py
│  └─── 03_inverse_matfile_task.py
└─── functions_preprocessing.py
└─── functions_make_characteristics.py
└─── fnames.py
└─── config.py

└─── Matlab_scripts
│  └───dual_model -> Individual subject and parcel specific model training
│    └───get_slt.m
│    └─── get_data_paths.m
│    └─── make_dual_model.m
│    └─── get_options_dual.m
│    └─── bash_dual_run.sh
│    └─── list_of_subjects.txt
│  └─── group_model -> group model training
│    └─── func_make_group_model.m
│    └─── list_of_subjects.txt
│    └─── make_group_model.m
│    └─── get_data_paths.m
│    └─── get_options_group.m
│    └─── bash_group_run.sh

└─── group_model_1.mat -> This is the pretrained group model presented in the manuscript
```
Python scripts contains basically the MEG data preprocessing steps and data transformation from the sensor level to the source level. In addition the Python scipts contains the plotting routines and implementations for the statistical testing.

Matlab codes includes the scripts to make HMM group model and the indivudual trained models.

Data is assumed to be organized as follows:

```

└─── MRI -> MRI images of the subject
|  └─── subject-01
└─── MEG_derivatives -> Processed files of the subjects
|  └─── subject-01
|    └─── ses-01
|      └─── noise_cov
|      └─── forward
|      └─── HMM_output
|      └─── stcs
|      └─── annot
|      └─── inverse
|      └─── ica
|      └─── data_mat
|      └─── MEG .fif files
|  └─── characteristics_csvs -> CSVs of the characteristics
|    └─── EA_csv
|    └─── SLT_csv
|    └─── dispersion_csv
|    └─── rate_csv
|    └─── FO_csv
└─── HMM_derivatives -> Output from the HMM analysis
|  └─── subject-01
|    └─── ... SLT, hmm, FO etc. output files
|____group_model_1.mat -> group model

```

The codes were developed using Python version 3.10.14 and mne-python version 1.7.0.

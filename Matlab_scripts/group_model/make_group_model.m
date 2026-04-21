%Add edited path
disp('Starting the process')

% Set identificator for the group model and load the group model
group_id = "group_model_1"
k = 4 %Number of states

% Get paths defined in get_data_paths()
[HMM_path, group_model_path, data_path, matlab_scripts_path]  = get_data_paths();

% Get the HMM functions
addpath(genpath(strcat(HMM_path, 'HMM-MAR'))) % Biomag-path

% Load the subjects
subjects_path =  strcat(matlab_scripts_path, '/list_of_subjects.txt')
subjects = importdata(subjects_path)

% Set sensortype and downsampling frequency
sessions = {'01'}
N_subjects = length(subjects)
N_sessions = length(sessions)
fs = 200
lfreq = 13
hfreq = 30
task='restEO'

%Make the data array
X_array = {}
T_array = {}

%Load the data for each subject and for the session 01
for sub_index = 1:N_subjects
    for ses_index = 1:N_sessions

        processed_path = strcat(data_path, subjects{sub_index}, '/ses-', sessions{ses_index},'/');
        X_name = strcat(processed_path,'/data_mat/', subjects{sub_index}, '_lfreq-',string(lfreq), '_hfreq-', string(hfreq),'_task-',task,'_mat.mat');

        X = load(X_name);
        X = X.x;

        sz = size(X);
        T = sz(1);

        X = reshape(X,[],1);
        X_array{end+1,1} = X;
        T_array{end+1,1} = zeros(1,sz(2))+T;

    end
end
disp('Data loaded')

% Run the HMM group-level model
[hmm_tde,gamma_tde,xi_tde,vi_tde,spectra_tde,options] = func_group_concat_sens(X_array,T_array,fs,k,lfreq,hfreq);

% Set the outputs to be saved
tde_hmm_r = hmm_tde;
gamma_r = gamma_tde;
Vipath_r = vi_tde;
freq_r = spectra_tde.state(1).f(:);
spect_tde_r = spectra_tde;

%Save the group model
save_path = strcat(group_model_path, '/',group_id ,'_73.mat');
save_path2 = strcat(group_model_path, '/',group_id ,'.mat');

save(save_path,"tde_hmm_r","gamma_r","Vipath_r","freq_r","spect_tde_r","T_array", "options","-v7.3")
save(save_path2,"tde_hmm_r","gamma_r","Vipath_r","freq_r","spect_tde_r","T_array", "options")


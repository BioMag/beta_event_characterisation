%Add edited path
clear all
disp('Starting the process')

% Get paths defined in get_data_paths()
[HMM_path, group_model_path, data_path, matlab_scripts_path]  = get_data_paths();

% Get the HMM functions
addpath(genpath(strcat(HMM_path, 'HMM-MAR'))) % Biomag-path

% Set identificator for the group model and load the group model
group_id = "job_source_2"
group_hmm = load(strcat(group_model_path,'sensors_concat_group_',group_id,'.mat')).tde_hmm_r;

%save_hmm = '/projects/HMM-beta/HMM_beta_sara/processed/group/sensors_concat/group_hmm.mat'
%save(save_hmm,"group_hmm")

% Set important parameters
sessions = {'01'}
N_sessions = length(sessions)

task = 'restEO'
fs = 200
lfreq = 13
hfreq = 30
k = 4 %Number of states
dual_id = 'group_model_1'

for task = ['restEO', 'leftCKCevoked']
    for ses_index=1:N_sessions
        % Read the subjects   
        subjects_path = strcat(matlab_scripts_path, 'list_of_subjects.txt');
        subjects = importdata(subjects_path);
        N_subjects = length(subjects);
        
        %Make empty data arrays
        X_array = {}
        T_array = {}
        subjects_incl = {}
     
        % Load the data for each subject
        for sub_index = 1:N_subjects
            processed_path = strcat(data_path, subjects{sub_index}, '/ses-', sessions{ses_index},'/');
            X_name = strcat(processed_path,'/data_mat/', subjects{sub_index}, '_lfreq-',string(lfreq), '_hfreq-', string(hfreq),'_task-',task ,'_mat.mat');
            
            if isequal(exist(X_name),2)
                X = load(X_name);
                X = X.x;
                sz = size(X);
                T = sz(1);
        
                X = num2cell(X,1);
                X_array{end+1,1} = X;
                T_array{end+1,1} = num2cell(zeros(1,sz(2))+T ,1);
                subjects_incl{end+1,1} = subjects{sub_index};
            end
        end
    
        % Correct the subjects and N_subjects if data was not found from
        % everyone
        subjects = subjects_incl
        N_subjects = length(subjects)
        
        % Load the options
        [options] = get_options_dual(fs,k,lfreq,hfreq);
        
        % Then iterate over subjects and parcels separately
        for sub_index = 1:N_subjects
            disp(sub_index)
            X_subject = X_array{sub_index};
            T_subject = T_array{sub_index};
    
            output = struct();
            output_SLT = struct();
            FOs_subject = {};
            output_TPM = struct();
    
            n_sens_labels = size(X_subject);
    
            for sensor = 1:n_sens_labels(2)
                sens_id = append('sens', int2str(sensor));
                X_sensor = X_subject{sensor};
                T_sensor = T_subject{sensor};
                
                [hmm_elem,Gamma_elem,vpath_elem, Xi_elem, LL_elem, Xt_elem] = hmmdual(X_sensor,T_sensor,group_hmm);
                spectra_dual = hmmspectramt(X_sensor,T_sensor,Gamma_elem,options);
                spectra_original = hmmspectramt(X_sensor,T_sensor,Gamma_elem,options);
        
                output.(sens_id).spectra = spectra_dual;
                output.(sens_id).spectra_orig = spectra_original;
                output.(sens_id).viterbi = vpath_elem;
                output.(sens_id).hmm = hmm_elem;
                output.(sens_id).gamma = Gamma_elem;
    
                % Get the SLT, FO and TPM to save some time (loading dual
                % takes time)
                output_SLT = get_slt(vpath_elem,T_sensor,sens_id, options, output_SLT);
                frac_oc = getFractionalOccupancy(Gamma_elem,T_sensor,[],2);
                FOs_subject{end+1,1} = frac_oc;    
                output_TPM.(sens_id).tpm = hmm_elem.P ;
    
            end
            
            save_path = strcat(data_path, subjects{sub_index}, '/ses-', sessions{ses_index},'/HMM_output/');
            
            % Define the names for the saved files
            save_file = strcat(save_path, '/dual_hmm_',dual_id,'_task_',task,'.mat' );
            save_path_TPM = strcat(save_path, '/TPM_parameters_',dual_id,'_task_',task,'.mat');
            save_path_SLT = strcat(save_path , '/SLT_parameters_',dual_id,'_task_',task,'.mat');
            save_path_FO = strcat(save_path, '/FO_',dual_id,'_task_',task,'.mat');
    
            if isequal(exist(save_path, 'dir'),7)
                save(save_file,'output')
                save(save_path_TPM,'output_TPM')
                save(save_path_SLT,'output_SLT');
                save(save_path_FO,'FOs_subject');
            else
                mkdir(save_path)
                save(save_file,'output')
                save(save_path_TPM,'output_TPM')
                save(save_path_SLT,'output_SLT');
                save(save_path_FO,'FOs_subject');
            end
        end
    end
end


function [op] = get_options_dual(fs,k,lfreq,hfreq)
% Get options for individual level HMM models

    options = struct();
    options.K = k;

    %For TDE we need to specify
    options.order=0;
    options.covtype='full';
    options.zeromean=1;
    options.dropstates = 0;
    options.Fs = fs;
    options.standardise = 1;
    options.standardise_pc = options.standardise;
    
    %Use the stochastic inference scheme
    % In the dual we don't want to use the stochastic scheme

    %Specifically for TDE
    options.embeddedlags=-8:8;
    options.useParallel = 0;
    options.pca = 6;

    %Non-parametric spectral information
    %We set the options for the spectral estimation:
    options.fpass = [lfreq hfreq]; % frequency range we want to look at
    options.tapers = [4 7]; % internal multitaper parameter
    options.win = 5 * options.Fs; % window length, related to the level of detail of the estimation;
    
    op = options;
    
end


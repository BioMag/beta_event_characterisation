function [op] = get_options_group(fs,k,lfreq,hfreq)

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
    options.BIGNinitbatch = 5;
    options.BIGNbatch = 5;
    options.BIGtol = 1e-7;
    options.BIGcyc = 500;
    options.BIGundertol_tostop = 5;
    options.BIGdelay = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;

    %Specifically for TDE
    options.embeddedlags=-8:8;
    options.useParallel = 0;
    options.pca = 6;

    %Non-parametric spectral information
    %We set the options for the spectral estimation:
    options.fpass = [lfreq hfreq]; % frequency range we want to look at, in this case between 1 and 40 Hertzs.
    options.tapers = [4 7]; % internal multitaper parameter
    options.win = 5 * options.Fs; % window length, related to the level of detail of the estimation;
    
    op = options;
    
end


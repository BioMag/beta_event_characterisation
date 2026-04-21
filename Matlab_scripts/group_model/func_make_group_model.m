function [hmm,gamma,xi,vi,spectra,op] = func_make_group_model(X_array,T_array,fs,k,lfreq,hfreq)

    [options] = get_options_group(fs,k,lfreq,hfreq);
    
    [hmm_tde,Gamma_tde,Xi_tde,Vipath_tde] = hmmmar(X_array,T_array,options);

    spectra_tde = hmmspectramt(X_array,T_array,Gamma_tde,options);

    hmm = hmm_tde;
    gamma = Gamma_tde;
    xi = Xi_tde;
    vi = Vipath_tde;
    spectra = spectra_tde;
    op = options;
    
end


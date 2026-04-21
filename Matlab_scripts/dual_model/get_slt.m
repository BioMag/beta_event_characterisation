function [output_SLT] = get_slt(dual_vi,T_array,sens_id, options, output_SLT_in)
%GET_SLT Calculate the state life time characteristics

    output_SLT = output_SLT_in;
    
    slt = getStateLifeTimes(dual_vi,T_array,options,10); %Threshold is 10 datapoints which means 50ms
    
    output_SLT.(sens_id).state1.mean = mean(slt{1});
    output_SLT.(sens_id).state2.mean = mean(slt{2});
    output_SLT.(sens_id).state3.mean = mean(slt{3});
    output_SLT.(sens_id).state4.mean = mean(slt{4});
    
    output_SLT.(sens_id).state1.median = median(slt{1});
    output_SLT.(sens_id).state2.median = median(slt{2});
    output_SLT.(sens_id).state3.median = median(slt{3});
    output_SLT.(sens_id).state4.median = median(slt{4});
    
    output_SLT.(sens_id).state1.std = std(slt{1});
    output_SLT.(sens_id).state2.std = std(slt{2});
    output_SLT.(sens_id).state3.std = std(slt{3});
    output_SLT.(sens_id).state4.std = std(slt{4});
    
    output_SLT.(sens_id).state1.max = max(slt{1});
    output_SLT.(sens_id).state2.max = max(slt{2});
    output_SLT.(sens_id).state3.max = max(slt{3});
    output_SLT.(sens_id).state4.max = max(slt{4});
    
    n_robust_1 = round(0.05*length( slt{1} ));
    n_robust_2 = round(0.05*length( slt{2} ));
    n_robust_3 = round(0.05*length( slt{3} ));
    n_robust_4 = round(0.05*length( slt{4} ));
    
    output_SLT.(sens_id).state1.robust_max = mean( maxk(slt{1},n_robust_1) );
    output_SLT.(sens_id).state2.robust_max = mean( maxk(slt{2},n_robust_2) );
    output_SLT.(sens_id).state3.robust_max = mean( maxk(slt{3},n_robust_3) );
    output_SLT.(sens_id).state4.robust_max = mean( maxk(slt{4},n_robust_4) );

end


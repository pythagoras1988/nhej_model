function time = compute_next_time_pre_ligation(stateObj,rate_constant_pre_ligation)
    %rate_constant_pre_ligation = 1./[2100,180000];
    r = rand(); 
    pre_ligation = stateObj.stateArr; 
    ind = find(pre_ligation==1);
    switch ind
        case 1 % Initial synapse complex 
            time = 1/(rate_constant_pre_ligation(1) + 2*rate_constant_pre_ligation(2)) * log(1/r);
        otherwise
            time = inf; 
    end
    % Output time will be in seconds
end
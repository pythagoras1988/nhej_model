function time = compute_next_time_z2(stateObj,rate_constant_z2)
    %rate_constant_z2 = 1./[120,10,2100,180000];
    r = rand(); 
    z2 = stateObj.stateArr; 
    ind = find(z2==1);
    switch ind
        case 1 % Initial synapse complex 
            time = 1/rate_constant_z2(1) * log(1/r);
        case 2 % Repairing dirty end in synapse
            time = 1/(2*rate_constant_z2(2) + rate_constant_z2(4)) * log(1/r);
        otherwise
            time = inf;
    end
    % Output time will be in seconds
end
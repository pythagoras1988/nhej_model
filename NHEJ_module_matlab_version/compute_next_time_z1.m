function time = compute_next_time_z1(stateObj,rate_constant_z1)
    %rate_constant_z1 = 1./[2100,180000];
    r = rand(); 
    z1 = stateObj.stateArr; 
    ind = find(z1==1);
    switch ind
        case 1 % Initial synapse complex 
            time = inf;
        otherwise
            time= inf; 
    end
    % Output time will be in seconds
end
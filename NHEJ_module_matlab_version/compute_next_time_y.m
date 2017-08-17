function time = compute_next_time_y(stateObj,rate_constant_y)
    %rate_constant_y  = 1./[0.1,0.2,3,3,500,1,60];
    y = stateObj.stateArr;
    r = rand(); 
    ind = find(y==1);
    switch ind
        case 1 % Raw end
            time = 1/rate_constant_y(1) * log(1/r); 
        case 2 % Opened end
            time = 1/sum(rate_constant_y([2,3])) * log(1/r); 
        case 3 % repairing opened end 
            time = 1/rate_constant_y(4) * log(1/r);  
        case 4 % ku70/80 attached
            time = 1/sum(rate_constant_y([5,6])) * log(1/r);  
        case 5 % DNA-PKcs attached 
            time = 1/rate_constant_y(7) * log(1/r);  
        otherwise
            time = inf;
    end
    % Output time will be in seconds
end
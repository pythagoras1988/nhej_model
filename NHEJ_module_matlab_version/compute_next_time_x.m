function time = compute_next_time_x(stateObj,rate_constant_x)
    r = rand(); 
    x = stateObj.stateArr; 
    ind = find(x==1);
    switch ind
        case 1 % Raw end
            time = 1/rate_constant_x(1) * log(1/r); 
        case 2 % Opened end
            time = 1/rate_constant_x(2) * log(1/r); 
        case 3 % ku 70/80 attached
            time = 1/sum(rate_constant_x([3,4])) * log(1/r);  
        otherwise
            time = inf; 
    end
    % Output time will be in seconds
end
function out_state = compute_next_state_x(stateObj,rate_constant_x,dt) 
    r1 = rand(); 
    x = stateObj.stateArr; 
    ind = find(x==1);
    switch ind
        case 1 % Raw end
            if r1<rate2prob(rate_constant_x(1),dt)
                stateObj.stateArr = [0,1,0,0];
                out_state = stateObj;
                return;
            end
        case 2 % Opened end
            if r1 < rate2prob(rate_constant_x(2),dt)
                stateObj.stateArr = [0,0,1,0];
                out_state = stateObj;
                return;
            end
        case 3 % ku 70/80 attached
            if r1 < rate2prob(rate_constant_x(3),dt)
                stateObj.stateArr = [0,1,0,0]; 
                out_state = stateObj;
                return; 
            elseif r1 >= rate2prob(rate_constant_x(3),dt) && ...
                r1 <= (rate2prob(rate_constant_x(3),dt) + rate2prob(rate_constant_x(4),dt))
                stateObj.stateArr = [0,0,0,1];
                out_state = stateObj;
                return;
            end
    end
    out_state = stateObj; 
end
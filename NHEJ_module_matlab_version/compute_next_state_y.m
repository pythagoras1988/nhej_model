function out_state = compute_next_state_y(stateObj,rate_constant_y,dt) 
    r1 = rand(); 
    y = stateObj.stateArr; 
    ind = find(y==1);
    switch ind
        case 1 % Raw end
            if r1<rate2prob(rate_constant_y(1),dt)
                stateObj.stateArr = [0,1,0,0,0];
                out_state = stateObj;
                return;
            end
        case 2 % Opened end
            if r1 < rate2prob(rate_constant_y(2),dt)
                stateObj.stateArr = [0,0,1,0,0];
                out_state = stateObj;
                return;
            elseif r1 >= rate2prob(rate_constant_y(2),dt) && ... 
                    r1 < rate2prob(rate_constant_y(2),dt)+ rate2prob(rate_constant_y(3),dt)
                stateObj.stateArr = [0,0,0,1,0];
                out_state = stateObj;
                return;                
            end
        case 3 % repairing opened end
            if r1 < rate2prob(rate_constant_y(4),dt)
                stateObj.stateArr = [0,1,0,0,0]; 
                out_state = stateObj;
                return;
            end
        case 4 % ku70/80 attached
            if r1 < rate2prob(rate_constant_y(5),dt)
                stateObj.stateArr = [0,1,0,0,0]; 
                out_state = stateObj;
                return;
            elseif r1 >= rate2prob(rate_constant_y(5),dt) && ...
                r1 < rate2prob(rate_constant_y(5),dt) + rate2prob(rate_constant_y(6),dt)
                stateObj.stateArr = [0,0,0,0,1];
                out_state = stateObj;
                return;
            end
        case 5 % DNA-PKcs attached
            if r1 < rate2prob(rate_constant_y(7),dt)
                stateObj.stateArr = [0,1,0,0,0]; 
                out_state = stateObj;
                return;
            end
    end
    out_state = stateObj; 
end
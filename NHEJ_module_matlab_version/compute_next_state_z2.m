function out_state = compute_next_state_z2(stateObj,rate_constant_z2,dt) 
    r1 = rand(); 
    z2 = stateObj.stateArr; 
    ind = find(z2==1);
    switch ind
        case 1 % Initial synapse complex
            if r1<rate2prob(rate_constant_z2(1),dt)
                stateObj.stateArr = [0,1,0];
                out_state = stateObj;
                return;
            end
        case 2 % Repairing dirty end in synapse
            if r1 < rate2prob(rate_constant_z2(2),dt)
                stateObj.stateArr = [0,1,0]; 
                out_state = stateObj;
                return;
            elseif r1 >= rate2prob(rate_constant_z2(2),dt) && ... 
                    r1 < 2*rate2prob(rate_constant_z2(2),dt)
                stateObj.stateArr = [0,0,1]; 
                out_state = stateObj;
                return;
            elseif r1 >= 2*rate2prob(rate_constant_z2(2),dt) && ... 
                    r1 < 2*rate2prob(rate_constant_z2(2),dt)+rate2prob(rate_constant_z2(4),dt)
                stateObj.stateArr = [1]; % Go to REsidual
                stateObj.stateIndex = 5;
                out_state = stateObj;
                return;
            end      
    end
    out_state = stateObj; 
end
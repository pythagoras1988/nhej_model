function out_state = compute_next_state_pre_ligation(stateObj,rate_constant_pre_ligation,dt) 
    r1 = rand(); 
    pre_ligation = stateObj.stateArr; 
    ind = find(pre_ligation==1);
    switch ind
        case 1 % Initial synapse complex
            if r1<rate2prob(rate_constant_pre_ligation(1),dt)
                stateObj.stateArr = [0,1]; %complete repair
                out_state = stateObj;
                return;
            elseif r1>=rate2prob(rate_constant_pre_ligation(1),dt) && ...
                    r1<rate2prob(rate_constant_pre_ligation(1),dt) + rate2prob(rate_constant_pre_ligation(2),dt)
                stateObj.stateArr = [0,1,0,0,0];
                stateObj.stateIndex = 2;
                out_state = stateObj;
                return;
            elseif r1>=rate2prob(rate_constant_pre_ligation(1),dt) + rate2prob(rate_constant_pre_ligation(2),dt) && ...
                r1<rate2prob(rate_constant_pre_ligation(1),dt) + 2*rate2prob(rate_constant_pre_ligation(2),dt)
                stateObj.stateArr = [1];
                stateObj.stateIndex = 5;
                out_state = stateObj;
                return;
            end
    end
    out_state = stateObj; 
end
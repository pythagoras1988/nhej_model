function out_state = compute_next_state_z1(stateObj,rate_constant_z1,dt) 
    r1 = rand(); 
    z1 = stateObj.stateArr; 
    ind = find(z1==1);
    switch ind
        case 1 % Initial synapse complex 
            stateObj.stateArr = [1];
            out_state = stateObj;
            return;
    end
    out_state = stateObj; 
end
function error_message(state)
    % Check for no more than 2 equal synapse ID
    synapse_ID_dat = zeros(1,ceil(length(state)/2));
    
    for k = 1:length(state)     
        if state(k).synapse_ID ~= -1
            synapse_ID_dat(state(k).synapse_ID) = synapse_ID_dat(state(k).synapse_ID) + 1;
        end
    end 
    if any(synapse_ID_dat>2) 
        error('More than 2 similar synapse ID.... '); 
    end
    
    % Check for correct state and state Index 
    error_flag = 0; 
    for  k = 1:length(state) 
        switch state(k).state.stateIndex 
            case 1
                if length(state(k).state.stateArr)~= 4 
                    error_flag = 1; 
                end
            case 2
                if length(state(k).state.stateArr)~= 5                   
                    error_flag = 1; 
                end
            case 3 
                if length(state(k).state.stateArr)~= 1 
                    error_flag = 1; 
                end
            case 4
                if length(state(k).state.stateArr)~= 3 
                    error_flag = 1; 
                end
            case 5
                if length(state(k).state.stateArr)~= 1 
                    error_flag = 1; 
                end
            case 6
                if length(state(k).state.stateArr)~= 2
                    error_flag = 1; 
                end
        end
    end
    if error_flag 
        error('State and state Index mismatch, please rectify!')
    end
end

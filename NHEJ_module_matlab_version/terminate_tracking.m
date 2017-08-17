function output = terminate_tracking(state) 
    switch state.stateIndex
        case 5
            if isequal([1],state.stateArr)
                output = 1; 
                disp('Terminate Tracking for state res...'); 
                return;
            end
        case 6
            if isequal([0,1],state.stateArr)
                output = 1; 
                disp('Terminate Tracking for state z2...'); 
                return;
            end
    end
    output = 0; % If none of the conditions are met
end
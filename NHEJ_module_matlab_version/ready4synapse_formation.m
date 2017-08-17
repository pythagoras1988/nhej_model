function output = ready4synapse_formation(state) 
    switch state.stateIndex
        case 1
            if isequal([0,0,0,1],state.stateArr)
                output = 1; 
                return; 
            end
        case 2
            if isequal([0,0,0,0,1],state.stateArr)
                output = 1; 
                return; 
            end
    end
    output = 0; 

end
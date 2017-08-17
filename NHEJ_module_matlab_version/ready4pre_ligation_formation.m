function output = ready4pre_ligation_formation(state) 
    switch state.stateIndex
        case 3
            if isequal([1],state.stateArr)
                output = 1; 
                return; 
            end
        case 4
            if isequal([0,0,1],state.stateArr)
                output = 1; 
                return; 
            end
    end
    output = 0; 
end
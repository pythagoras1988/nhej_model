function state = mat2state(mat) 
    len = length(mat); 
    switch len 
        case 4 
            state = 'x'; 
        case 5 
            state = 'y'; 
        case 2 
            state = 'z1'; 
        case 
end
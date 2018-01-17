function output_mat = remove_repeat_damage(input_mat) 
    % input_mat = N x 9 matrix 
    % This function always delete the second and subsequent repeating
    % damage basepairs
    if isempty(input_mat)
        output_mat =  input_mat;
    end    
     
    len = length(input_mat(:,1)) ; 
    tmpMat = zeros(len,1);
    for k = 2:len 
        if isequal(input_mat(k,1:3),input_mat(k-1,1:3)) 
            tmpMat(k) = 1;
            disp('-------------------------------------------------------'); 
            disp('One repeating data deleted.....'); 
            disp('-------------------------------------------------------');
        end
    end
    input_mat(logical(tmpMat),:) = []; 
    output_mat = input_mat;
end
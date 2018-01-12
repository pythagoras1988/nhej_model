function out_mat = consolidate_direct_dmg(in_mat) 
    % Use this code only when full edep file is used for energy deposition.
    % It adds up all the contribution of direct damage to similar basepair 
    if isempty(in_mat) 
        out_mat = in_mat; 
        return;
    end
    
    temp = in_mat(1,:); 
    tmpMat = zeros(length(in_mat(:,1)),1);
    temp_ind = 1; 
    
    for k = 2:length(in_mat(:,1))
        if isequal(in_mat(k,1:3),temp(1,1:3)) && (in_mat(k,7)==0) 
            in_mat(temp_ind,9) = in_mat(temp_ind,9) + in_mat(k,9);
            tmpMat(k) = 1; 
            disp('-------------------------------------------------------'); 
            disp('MERGED.....'); 
            disp('-------------------------------------------------------');
        else 
            temp = in_mat(k,:); 
            temp_ind = k; 
        end
    end
    in_mat(logical(tmpMat),:) = []; 
    out_mat = in_mat; 
end
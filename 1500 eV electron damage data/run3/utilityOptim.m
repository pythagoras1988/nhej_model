function val = utilityOptim(mat) 
    [lenX lenY] = size(mat); 
    reqMember   = min(lenX,lenY); 
    totMember   = sum(sum(mat==1));
    [ind1 ind2] = find(mat==1); 
    
    sizeInd = size(ind1); 
    
    if sizeInd(2)>sizeInd(1) 
        indMat = [ind1',ind2'];
    else
        indMat = [ind1,ind2];
    end
    maxCheck = 0; 
    for k = 1:2^totMember-1
        bitRep = de2bi(k);
        tmp = indMat;
        if sum(bitRep)<=reqMember 
            tmp(~bitRep',:) = [] ;
            check = sameGroup(tmp); 
            if check && length(tmp(:,1))>maxCheck 
                val = tmp; 
                maxCheck = length(tmp(:,1)); 
            end
        end
    end     
end

function val = sameGroup(mat) 
    [len1,len2] = size(mat);
    if length(unique(mat(:,1)))<len1 || length(unique(mat(:,2)))<len1
        val =  0; 
    else 
        val =  1; 
    end
end

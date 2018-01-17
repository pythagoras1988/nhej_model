clear all; close all; clc; 

files = dir('*.txt'); 
count = 1; 
for k = 1:length(files) 
    if (1)
        if files(k).name(1) == 'd' 
            dat = load(files(k).name); 
            tmp = dat(:,7)==1; 
            dat(~tmp,:) = []; 
            if ~isempty(dat) 
                numChem(count) = dat(1,9);
                numChemDmg(count) = length(dat(:,1));
                count = count + 1; 
            end
        end
    else
         if files(k).name(1) == 'c' && files(k).name(end-4)=='d'
            dat = load(files(k).name); 
            tmp = unique(dat(:,2));  
            numChem(count) = length(tmp);
            count = count + 1; 
         end
    end
    if rem(k,20) == 0
        k
    end
end

mean(numChem) 
std(numChem) 
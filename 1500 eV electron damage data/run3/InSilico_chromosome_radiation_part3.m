clear all; close all; clc; 

%% ***********************************************************************
% This will determine the number of direct vs indirect strands breaks as
% well as DSB and SSB 
%
% Date: 18/09/12016
%%%***********************************************************************

%% -----------------------------------------------------------------------
% File I/O and matrix conditioning
%%%-----------------------------------------------------------------------

txtFile = 1; 
newDataType = 1; 
simpleMethod = 0; 
file = dir('*.txt'); 
len = length(file)-1; 
nbParticle = 0; 
eDepTot = 0 ; 
count2   = 1; 
count1  = 0 ;
mass = 1.0 * (1.5^3) * 10^-12 * 10^-3; 
chromFactor = 1; 
edepFile = importdata('edepMaster.txt'); 
m = 1;

while(m<len)

    if ~txtFile 
        load('500_eV_electron.mat'); 
        dat = fullDmgMat;
    else
        dat = load(file(m).name);
        eDep(m) = edepFile(1+str2num(file(m).name(10:end-4)));
        eDep(m) = eDep(m) / (chromFactor * mass) * 1.6 *10^-19; 
        eDepTot = eDepTot + eDep(m);
    end 

    %Condition the data
    if newDataType 
        tmp = []; 
        Ethresh = 16.5; 
        for k = 1:length(dat(:,1)) 
            if dat(k,7)==0 && dat(k,9)<Ethresh
                tmp(k) = 1;
            end 
        end
        if ~isempty(tmp)
            tmp = logical(tmp); 
            dat(tmp,:) = [];
        end
    end


    totalBreak      = length(dat(:,1)); 
    directBaseDmg   = sum(dat(:,7)==1 & dat(:,2)==0);
    indirectBaseDmg = sum(dat(:,7)==1 & dat(:,2)==0);
   

    %% -----------------------------------------------------------------------
    % Process SSB Vs DSB break 
    %%%-----------------------------------------------------------------------

    % Insertion sort 
    for k = 2:totalBreak
        tmp   = dat(k,:); 
        count = k; 
        while count > 1 && tmp(1)<dat(count-1,1) 
            dat(count,:) = dat(count-1,:); 
            count = count - 1; 
        end
        dat(count,:) = tmp; 
    end
    
    dat(dat(:,2)==0,:) = []; 
  dat(dat(:,8)>3.2 & dat(:,7)==0,:) = []; 
  %dat(dat(:,8)>2.8 & dat(:,7)==1,:) = []; 
   % dat(dat(:,8)>3.2,:) = []; 
    %dat(dat(:,7)==1,:) = []; 
    
    directBreak     = sum(dat(:,7)==0 & dat(:,2)==1);
    indirectBreak   = sum(dat(:,7)==1 & dat(:,2)==1);

    % DBSCAN algorithm 
    epsilon = 10; 
    MinPts  = 2; 
    X       = dat(:,1); 
    [IDX, isnoise]=DBSCAN(X,epsilon,MinPts);
    %plot(dat(:,1),IDX,'ro');
    clustered_lesion = max(IDX); 

    % Cluster Analysis for SSB and DSB 
    SSB_isolated = sum(IDX==0); 
    tmp = dat(IDX==0,:);
    SSB_direct   = sum(tmp(:,7)==0 & tmp(:,2)==1);
    SSB_indirect = sum(tmp(:,7)==1 & tmp(:,2)==1);
    DSB_direct   = 0; 
    DSB_indirect = 0; 
    
    if (max(IDX) ~= 0 & ~isempty(IDX))
        for k = 1:max(IDX) 
            tmp     = dat(IDX==k,:);            
            break_A = sum(tmp(:,3) == 65); 
            break_B = sum(tmp(:,3) == 66);
            DSB(k)  = min(break_A,break_B); 
            A_mat   = tmp(tmp(:,3) == 65,:); 
            B_mat   = tmp(tmp(:,3) == 66,:); 
            
            
            if DSB(k) >= 1
                for l = 1:break_A
                    for ll = 1:break_B 
                        clusterMat(l,ll) = abs(A_mat(l,1) - B_mat(ll,1));
                    end
                end
                clusterMat = clusterMat < epsilon; 
                if sum(clusterMat(:)) > 0 
                    if simpleMethod
                        DSB(k) = 1; 
                        [indX indY] = find(clusterMat == 1);
                        DSB_direct   = DSB_direct + int32(A_mat(indX(1),7)==0) + int32(B_mat(indY(1),7)==0);
                        DSB_indirect = DSB_indirect + int32(A_mat(indX(1),7)==1) + int32(B_mat(indY(1),7)==1);
                    else
                        outMat = utilityOptim(clusterMat); 
                        DSB(k) = length(outMat(:,1)); 
                        DSB_direct   = DSB_direct + sum(int32(A_mat(outMat(:,1),7)==0) + int32(B_mat(outMat(:,2),7)==0));
                        DSB_indirect = DSB_indirect + sum(int32(A_mat(outMat(:,1),7)==1) + int32(B_mat(outMat(:,2),7)==1));
                    end
                else
                    DSB(k) = 0; 
                end
            end
            
            % clusterDmgType: direct x indirect 
            clusterDmgType(k,1) = sum(tmp(:,4)==0);
            clusterDmgType(k,2) = sum(tmp(:,4)==1);
            SSB(k)   = break_A + break_B - 2*DSB(k); 
            clear clusterMat;
        end
    else 
        SSB = 0;
        DSB = 0;
        clustered_lesion = 0;
    end

    %% -----------------------------------------------------------------------
    % Summary printout 
    %%%-----------------------------------------------------------------------

    SSBtot = SSB_isolated;% + sum(SSB); 
    DSBtot = sum(DSB); 
%     fid = fopen('Summary_results.txt','w'); 
%     fprintf(fid,'Number of Clustered Lesion = %d \n', clustered_lesion); 
%     fprintf(fid,'Direct SSB break = %d \n', directBreak); 
%     fprintf(fid,'Indirect SSB break = %d \n', indirectBreak); 
%     fprintf(fid,'Total number of SSB = %d \n',SSBtot); 
%     fprintf(fid,'Total number of DSB = %d \n',DSBtot);
% 
%     fclose(fid); 
        totalData(1,m) = clustered_lesion; 
        totalData(2,m) = SSB_direct; %directBreak - DSB_direct; 
        totalData(3,m) = SSB_indirect; %indirectBreak - DSB_indirect; 
        totalData(4,m) = SSBtot; 
        totalData(5,m) = DSBtot; 
        totalData(6,m) = eDep(m); 
    
    if count1 == 0 
        masterData(1,count2) = 0;
        masterData(2,count2) = 0;
        masterData(3,count2) = 0;
        masterData(4,count2) = 0;
        masterData(5,count2) = 0;
        masterData(6,count2) = 0;
    end
        
    if eDepTot <1.05 
        masterData(1,count2) = masterData(1,count2) + clustered_lesion; 
        masterData(2,count2) = masterData(2,count2) + SSB_direct; %directBreak - DSB_direct; 
        masterData(3,count2) = masterData(3,count2) + SSB_indirect; %indirectBreak - DSB_indirect; 
        masterData(4,count2) = masterData(4,count2) + SSBtot; 
        masterData(5,count2) = masterData(5,count2) + DSBtot; 
        masterData(6,count2) = masterData(6,count2) + eDep(m); 
        count1 = count1 + 1; 
        m = m + 1; 
    else 
        masterData(1:5,count2) = masterData(1:5,count2) / chromFactor; 
        eDepTot = 0; 
        count2 = count2 + 1; 
        m = m - floor(count1/2); 
        count1 = 0;
    end
    clear tmp DSB SSB; 
    disp(num2str(m)); 

end

daltonPerBP = 650; 
numBP = 1575*80808;
if (0)
    unitFactor = 1575*80808/10^9; % in per Gray per Gbp
else
    unitFactor = (daltonPerBP * numBP) * 10^-11; % in per gray per dalton x 10^-11
end

disp(['Number of direct SSB = ' num2str(mean(masterData(2,1:end-1))/unitFactor) ', std = ' num2str(std(masterData(2,1:end-1))/unitFactor/sqrt(count2-1))]); 
disp(['Number of indirect SSB = ' num2str(mean(masterData(3,1:end-1))/unitFactor) ', std = ' num2str(std(masterData(3,1:end-1))/unitFactor/sqrt(count2-1))]);
disp(['Number of SSB = ' num2str(mean(masterData(4,1:end-1))/unitFactor) ', std = ' num2str(std(masterData(4,1:end-1))/unitFactor/sqrt(count2-1))]);
disp(['Number of DSB = ' num2str(mean(masterData(5,1:end-1))/unitFactor) ', std = ' num2str(std(masterData(5,1:end-1))/unitFactor/sqrt(count2-1))]);


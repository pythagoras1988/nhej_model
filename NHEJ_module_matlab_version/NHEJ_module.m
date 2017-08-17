clear all; close all; clc; 

%%%***********************************************************************
% This code aims to simulate and model the NHEJ repair process in cell
% according to the network of biochemical interactions. A stochastic
% approach will be used instead of an ODE approach
%
% Author : Higgsino 
% Date   : 3/1/2017 
%%%***********************************************************************

%% -----------------------------------------------------------------------
% Load data , I/O 
%%%-----------------------------------------------------------------------d

% Incorporate 23 chromosomes into NHEJ model for a complete nucleus model
chromosome_positions = importdata('chromosome_positions_final.txt');
gene_complex_threshold = 10; % Number of base pair, within which the damage is regarded as complex
nucl_radius = 40000; % in term of angstrom
count = 1;

for p = 1:length(chromosome_positions)
    dat    = importdata(strcat('damageMat',num2str(p),'.txt')); 
    offset = chromosome_positions(p,:) - 7500*[1 1 1];

    %% ------------------------------------------------------------------------
    % Preprocess the data to determine location and complexity of double strand
    % break 
    %%%------------------------------------------------------------------------

    % Insertion sort 
    for k = 2:length(dat(:,1))
        tmp   = dat(k,:); 
        count1 = k; 
        while count1 > 1 && tmp(1)<dat(count1-1,1) 
            dat(count1,:) = dat(count1-1,:); 
            count1 = count1 - 1; 
        end
        dat(count1,:) = tmp; 
    end

    dat(dat(:,8)>3.2 & dat(:,7)==0,:) = []; % Limit energy deposition to within van der waal radius 

    tmp = []; 
    Ethresh = 16.5; 
    for k = 1:length(dat(:,1)) 
        if dat(k,7)==0 && dat(k,9)<Ethresh
            tmp(k) = 1;
        end 
    end
    if ~isempty(tmp)
        tmp = logical(tmp); 
        dat(tmp,:) = []; % Energy > Ethresh to qualify as breaks
    end 

    % Create an array which contain strand break damage only! 
    dat_strand = dat;
    dat_base   = dat; 
    dat_strand(dat_strand(:,2)==0,:) = []; 
    dat_base(dat_base(:,2)==1,:)     = []; 

    % Classify DSB and SSB with DBSCAN algo
    epsilon = 10; 
    MinPts  = 2; 
    X       = dat_strand(:,1); 
    [IDX, isnoise]=DBSCAN(X,epsilon,MinPts);
    clustered_lesion = max(IDX); 

    SSB_isolated = sum(IDX==0); 
    DSB = []; % initialize
 
    for k = 1:max(IDX) 
        tmp     = dat_strand(IDX==k,:);            
        break_A = sum(tmp(:,3) == 65); 
        break_B = sum(tmp(:,3) == 66);
        DSB(k)  = min(break_A,break_B); 

        % Create an structure that contains all the information on the double strand
        % breaks : DSB_mat
        assert(DSB(k)<=1); 

        if DSB(k)==1 
            DSB_mat(count).info_A = tmp(tmp(:,3)==65,:); 
            DSB_mat(count).info_B = tmp(tmp(:,3)==66,:); 
            DSB_mat(count).posA = DSB_mat(count).info_A(1,4:6) + offset; 
            DSB_mat(count).posB = DSB_mat(count).info_B(1,4:6) + offset; 
            gene_ind_A = DSB_mat(count).info_A(1);
            gene_ind_B = DSB_mat(count).info_B(1);       
            gene_ind_min = min(gene_ind_A, gene_ind_B); 
            gene_ind_max = max(gene_ind_A, gene_ind_B); 
            DSB_mat(count).chrom_ID = p; 

            tmpMat_min = sum((dat(:,1)<gene_ind_min) & (dat(:,1)> gene_ind_min-gene_complex_threshold));
            tmpMat_max = sum((dat(:,1)>gene_ind_max) & (dat(:,1)< gene_ind_max+gene_complex_threshold));

            if gene_ind_A <gene_ind_B 
                DSB_mat(count).complex_A = (tmpMat_min>0) ;
                DSB_mat(count).complex_B = (tmpMat_max>0) ;
            else 
                DSB_mat(count).complex_A = (tmpMat_max>0) ;
                DSB_mat(count).complex_B = (tmpMat_min>0) ;
            end               
            count = count + 1;
        end
    end

    disp(['Number of chromosome is ', num2str(p), '...']);
    disp(['Number of Double strand breaks is ', num2str(sum(DSB)), '...']);
end
 
%% -----------------------------------------------------------------------
% NHEJ module Constructor
%%%-----------------------------------------------------------------------

% Define parameters for the stochastic simulation 
dt        = 0.5; % in seconds 
totalTime = 4; %in hours 
totalTime = totalTime * 3600; %Convert to seconds
numDSB    = sum(DSB); 
dump_count = 1; 
D1 = 170*100; % Diffusion constant in angstrom^2/s
D2 = D1 / 10; 

% Define the chemical species 
% x: The state of the clean end 
x = []; 
x(1) = 0; % Raw broken end 
x(2) = 0; % Opened end
x(3) = 0; % ku70/80 attached
x(4) = 0; % XL attached

% y: The state of the dirty end 
y = []; 
y(1) = 0; % Raw Broken end
y(2) = 0; % opened end 
y(3) = 0; % repairing opened end 
y(4) = 0; % ku70/80 attached
y(5) = 0; % DNA-PKcs attached 

% z1: post synapse clean end 
z1 = []; 
z1(1) = 0; % Initial synapse complex 

% z2: post synapse dirty end 
z2 = []; 
z2(1) = 0; % Initial synapse complex
z2(2) = 0; % Repairing dirty end in synapse
z2(3) = 0; % before ligation and complete repair of misc damage

% res: Residual breaks, irreparable
res = []; 
res(1) = 0; 

% pre_ligation: State just right before ligation and only formed when both
% ends have been thoroughly repaired
pre_ligation = []; 
pre_ligation(1) = 0; % Pre ligation state
pre_ligation(2) = 0; % Fully repaired

% Define all the rate constants
rate_constant_x  = 1./[0.1,3,500,120]; 
rate_constant_y  = 1./[0.1,0.2,3,3,500,1,60];
rate_constant_z1 = []; % Currently no rate constant for z1,
rate_constant_z2 = 1./[120,10,2100,180000]; % last one to residual break or post-synaptic failure scenario
rate_constant_pre_ligation = 1./[2100,180000];

% Define all arrays and parameters relating to synapse
synapse_prob = 0.95; 
synapse_dist = 250; % in terms of angstrom
synapse_ID_all = [0]; % Contain all synapse ID

%% Create container maps for all state 
% keySet   = {'x', 'y', 'z1', 'z2', 'res','pre_ligation'}; 
% valueSet = {x,y,z1,z2,res,pre_ligation}; 
% stateObj = containers.Map(keySet,valueSet); 
% stateObj('stateIndex') = 0; % 1=x,2=y,3=z1,4=z2,5=res,6=pre_ligation

%% Create structure for all state 
stateObj.stateIndex = 0; % 1=x,2=y,3=z1,4=z2,5=res,6=pre_ligation
stateObj.stateArr = []; 

%% -----------------------------------------------------------------------
% NHEJ module Stochastic program
%%%-----------------------------------------------------------------------     
 
% Initialize all DSB
for k = 1:length(DSB_mat)
    if DSB_mat(k).complex_A 
        stateObj.stateArr = [1,0,0,0,0];
        stateObj.stateIndex = 2; 
        DSB_mat(k).state_A = stateObj; % initialize to y state 
    else
        stateObj.stateArr = [1,0,0,0];
        stateObj.stateIndex = 1; 
        DSB_mat(k).state_A = stateObj;  % initalize to x state
    end
    
    if DSB_mat(k).complex_B 
        stateObj.stateArr = [1,0,0,0,0];
        stateObj.stateIndex = 2; 
        DSB_mat(k).state_B = stateObj; % initialize to y state 
    else
        stateObj.stateArr = [1,0,0,0];
        stateObj.stateIndex = 1; 
        DSB_mat(k).state_B = stateObj;  % initalize to x state
    end
end

% redefine the dsb matrix
count = 1; 
for k = 1:length(DSB_mat) 
    dsb_masterdata(count).pos = DSB_mat(k).posA;
    dsb_masterdata(count).complex = DSB_mat(k).complex_A;
    dsb_masterdata(count).state = DSB_mat(k).state_A;
    dsb_masterdata(count).genomic_index = (count+1)/2; % Store ID of original breaks 
    dsb_masterdata(count).synapse_ID = -1;
    dsb_masterdata(count).chrom_ID   = DSB_mat(k).chrom_ID ;
    count = count + 1; 
    dsb_masterdata(count).pos = DSB_mat(k).posB;
    dsb_masterdata(count).complex = DSB_mat(k).complex_B;
    dsb_masterdata(count).state = DSB_mat(k).state_B;
    dsb_masterdata(count).genomic_index = count/2; % Store ID of original breaks 
    dsb_masterdata(count).synapse_ID = -1;
    dsb_masterdata(count).chrom_ID   = DSB_mat(k).chrom_ID ;
    count = count + 1; 
end
clear DSB_mat; 

% Gillespie algorithm 
fid = fopen('nhej_log.txt','w');
currTime = 0; 
saveCounter = 0; 
while currTime < totalTime 
    % Debugging and error checking phase
    error_message(dsb_masterdata); 
    
    % Termination condition
    % Description: DSB that has reached state 5 or 6 and finished
    % terminating will be stored in dump_masterdata
    tmp = zeros(1,length(dsb_masterdata)); 
    for k = 1:length(dsb_masterdata) 
        if terminate_tracking(dsb_masterdata(k).state)
            tmp(k) = 1; 
            dump_masterdata(dump_count) = dsb_masterdata(k); 
            dump_count = dump_count + 1; 
        end
    end
    dsb_masterdata(logical(tmp)) = []; 
    clear tmp;
   
    % compute minimum time from all DSB to determine next time step 
    next_time_step = inf; 
    for k = 1:length(dsb_masterdata);                
        switch dsb_masterdata(k).state.stateIndex 
            case 1 
                dt = compute_next_time_x(dsb_masterdata(k).state,rate_constant_x); 
            case 2
                dt = compute_next_time_y(dsb_masterdata(k).state,rate_constant_y); 
            case 3 
                dt = compute_next_time_z1(dsb_masterdata(k).state,rate_constant_z1); 
            case 4
                dt = compute_next_time_z2(dsb_masterdata(k).state,rate_constant_z2); 
            case 6
                dt = compute_next_time_pre_ligation(dsb_masterdata(k).state,rate_constant_pre_ligation);
        end             
        if dt<next_time_step 
            next_time_step = dt; 
        end
    end
    
    % The choice of time here depends on the time scale of the diffusion;
    % To be changed when necessary ( To be updated ....)
    if (next_time_step == inf)
        next_time_step = 0.2; % in seconds
    end
    
    disp(['Next Time Step is ' num2str(next_time_step) ' seconds...']);
    
    for k = 1:length(dsb_masterdata); 
        % Compute next state for all DSB
        switch dsb_masterdata(k).state.stateIndex 
            case 1 
                dsb_masterdata(k).state = compute_next_state_x(dsb_masterdata(k).state,rate_constant_x,next_time_step); 
            case 2
                dsb_masterdata(k).state = compute_next_state_y(dsb_masterdata(k).state,rate_constant_y,next_time_step); 
            case 3 
                dsb_masterdata(k).state = compute_next_state_z1(dsb_masterdata(k).state,rate_constant_z1,next_time_step);
            case 4
                dsb_masterdata(k).state = compute_next_state_z2(dsb_masterdata(k).state,rate_constant_z2,next_time_step);
            case 6
                dsb_masterdata(k).state = compute_next_state_pre_ligation(dsb_masterdata(k).state,rate_constant_pre_ligation,next_time_step);
        end   
    end
    
    synapse_ID_tmp = [0]; 
    synapse_ID_ind = [0]; 
    for  k = 1:length(dsb_masterdata)     
        % Compute next position based on diffusion; Requires more editing
        if dsb_masterdata(k).synapse_ID == -1
            % ooOoo... This is for dsb state that has not formed a synapse
            % at all (isolated) ...ooOoo
            dsb_masterdata(k).pos = diffusion_model(dsb_masterdata(k).pos, D1, next_time_step, nucl_radius); 
            
        elseif ~any(synapse_ID_tmp == dsb_masterdata(k).synapse_ID)
            % ooOoo... This is for dsb state with synapse but is first
            % encountered in this for loop ...ooOoo
            dsb_masterdata(k).pos = diffusion_model(dsb_masterdata(k).pos, D2, next_time_step, nucl_radius);
            synapse_ID_tmp = [synapse_ID_tmp dsb_masterdata(k).synapse_ID]; % Store synapse_ID that has occured
            synapse_ID_ind = [synapse_ID_ind k];% Store which DSB shares the same synapse
            
        else
            %ooOoo... This is for dsb state with synapse but has the same synapse ID as a previous dsb state 
            % identified in the loop. The position of this state will be
            % equal to that of the same synapse ID ...ooOoo
            tmpInd = find(synapse_ID_tmp == dsb_masterdata(k).synapse_ID);
            dsb_masterdata(k).pos = dsb_masterdata(synapse_ID_ind(tmpInd)).pos;  
            %disp(num2str(dsb_masterdata(k).pos)); % for debugging use
            %disp(num2str(dsb_masterdata(synapse_ID_ind(tmpInd)).pos)); % for debugging use
        end
    end
    
    %-------------------------------------------------------------------
    % Synapse formation condition    
    %-------------------------------------------------------------------
    
    % ooOoo................ooOoo....................ooOoo................ 
    % Create array containing all positions of DSB ends that is ready to
    % form the synapse
    count = 1; 
    dsb_synapse_end_pos =zeros(1,4);
    for m = 1:length(dsb_masterdata) 
        if ready4synapse_formation(dsb_masterdata(m).state)
            dsb_synapse_end_pos(count,1:3) = dsb_masterdata(m).pos; % record position
            dsb_synapse_end_pos(count,4)   = m; % record index 
            count = count + 1; 
        end
    end
    % ooOoo................ooOoo....................ooOoo................
    
    % Check if any synapse will be formed and modify the state accordingly  
    m = 1;
    while(m<=length(dsb_synapse_end_pos(:,1))-1)
       for mm = m+1:length(dsb_synapse_end_pos(:,1))
           dist = norm(dsb_synapse_end_pos(m,1:3)-dsb_synapse_end_pos(mm,1:3)); 
           if dist < synapse_dist && rand() <0.95 
               %---------------Synapse formation!--------------------------
               % Update synapse ID for master data
               n  = dsb_synapse_end_pos(m,4); 
               nn = dsb_synapse_end_pos(mm,4);  
               dsb_masterdata(nn).synapse_ID = max(synapse_ID_all) + 1; % Create new synapse ID for newly formed synapse
               dsb_masterdata(n).synapse_ID  = max(synapse_ID_all) + 1; 
               synapse_ID_all = [synapse_ID_all max(synapse_ID_all)+1];  % Update master synapse ID 
               % Update state for master data
               if dsb_masterdata(n).state.stateIndex == 1
                   dsb_masterdata(n).state.stateIndex = 3; 
                   dsb_masterdata(n).state.stateArr = [1];
               else 
                   dsb_masterdata(n).state.stateIndex = 4; 
                   dsb_masterdata(n).state.stateArr = [1,0,0];
               end
               
               if dsb_masterdata(nn).state.stateIndex == 1
                   dsb_masterdata(nn).state.stateIndex = 3; 
                   dsb_masterdata(nn).state.stateArr = [1];
               else 
                   dsb_masterdata(nn).state.stateIndex = 4; 
                   dsb_masterdata(nn).state.stateArr = [1,0,0];
               end
               % Identify correct or incorrect synapse formation 
               if dsb_masterdata(n).genomic_index == dsb_masterdata(nn).genomic_index 
                   disp('------------------------------------------------------------'); 
                   disp(['| Synapse formed with index ' num2str(dsb_masterdata(n).synapse_ID)]);
                   disp(['| Correct rejoining for genomic index ' num2str(dsb_masterdata(nn).genomic_index)]); 
                   disp('------------------------------------------------------------'); 
               else
                   disp('------------------------------------------------------------'); 
                   disp(['| Synapse formed with index ' num2str(dsb_masterdata(n).synapse_ID)]);
                   disp('| Incorrect rejoining.....'); 
                   disp('------------------------------------------------------------'); 
               end
               
               % Replace the position of the synapse by the midpoint positions of both ends
               tmp_midpoint = (dsb_masterdata(nn).pos + dsb_masterdata(n).pos)/2;
               dsb_masterdata(nn).pos = tmp_midpoint; 
               dsb_masterdata(n).pos  = tmp_midpoint; 
               
               % Remove data after synapse formation to prevent repeated sampling    
               dsb_synapse_end_pos(mm,:) = []; 
               break;
               %-----------------------------------------------------------
           end
       end
       m = m + 1; 
    end
    % ooOoo................ooOoo....................ooOoo................
    
    %% --------------------------------------------------------------------
    % Formation of pre-ligation complex
    %%%--------------------------------------------------------------------
    for m = 1:length(dsb_masterdata) 
       if (ready4pre_ligation_formation(dsb_masterdata(m).state)) && (dsb_masterdata(m).synapse_ID ~= -1)
           for mm = m+1:length(dsb_masterdata) 
               if (dsb_masterdata(mm).synapse_ID == dsb_masterdata(m).synapse_ID) && ...
                       (ready4pre_ligation_formation(dsb_masterdata(mm).state)) 
                   % Converting both ends to "completed repair" state 
                   dsb_masterdata(m).state.stateArr = [1,0]; % First end
                   dsb_masterdata(m).state.stateIndex = 6; 
                   
                   dsb_masterdata(mm).state.stateArr= [1,0]; % Second complementary end
                   dsb_masterdata(mm).state.stateIndex = 6; 
                   break;
               end
           end
       end
    end
    
    
    %% --------------------------------------------------------------------
    % NHEJ Destructor
    %%%--------------------------------------------------------------------
   
    currTime = currTime + next_time_step; 
    disp(['*********** Current Time is ' num2str(currTime) ' seconds..... *********']);
    pause(0.05);
    
    if ~isempty(dsb_masterdata) && rem(saveCounter,8) == 0
        fprintf(fid,'%f \t',currTime);
        for k = 1:length(dsb_masterdata)
            fprintf(fid,'%d \t',dsb_masterdata(k).state.stateIndex);
        end
        fprintf(fid,'\n');
    end
    
    % If no more dsb ends to step forward in time, terminate the simulation
    if isempty(dsb_masterdata); 
        disp('Finished simulation! All ends are rejoined'); 
        save('remaining_dsb_masterdata.mat','dsb_masterdata');
        save('Finished_dsb_masterdata.mat','dump_masterdata');
        break;        
    end
    saveCounter = saveCounter  + 1;
end
fclose(fid)
save('remaining_dsb_masterdata.mat','dsb_masterdata');
save('Finished_dsb_masterdata.mat','dump_masterdata');






 
 
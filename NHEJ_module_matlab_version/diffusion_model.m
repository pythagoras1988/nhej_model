%% This is one of the main module of the IRSGS program. It determines
% the diffusion of the dsb ends under different conditions. The
% different condition can be activated using a condition_flag. The
% conditions are as follows 
%
% condition_flag=1: Simple implemention with reflective BC 
% condition_flag=2: with Subdiffusion character
% condition_flag=3: with constrained spatial diffusion due to
%                   attachment to centromere
%
% Return: new calculated positions 
%
% Date: 24/1/2017 
%%%----------------------------------------------------------------------

function newPosition = diffusion_model(currPosition,diff_coeff,dt,nucl_radius) 
    % Generate new step size based on smoluchowski equation 
    stepSize = sqrt(2*diff_coeff*dt)*randn(1,3); 
    tmpPosition = currPosition + stepSize; 
    
    % Implementation of reflective boundary condition
    %release_dist = 50; % in angstrom 
    if norm(tmpPosition-nucl_radius*[1,1,1])<nucl_radius 
        newPosition = tmpPosition; 
        return; 
    else  
        step = linspace(0.01,1,20);
        step = fliplr(step); 
        for k = 1:length(step) 
            if norm(currPosition + step(k)*stepSize-nucl_radius*[1,1,1])<nucl_radius
                newPosition  = currPosition + step(k)*stepSize; 
                release_dist = nucl_radius - norm(currPosition + step(k)*stepSize-nucl_radius*[1,1,1]);
                disp(['Distance from boundary = ' num2str(release_dist) ' Angstrom']); 
                return; 
            end
        end 
        disp('Reflecting boundary implemented!');  
    end
    newPosition = currPosition;
end

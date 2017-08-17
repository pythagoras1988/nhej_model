clear all; close all; clc; 

dat = importdata('nhej_log.txt'); 
time = dat(:,1); 

for k = 1:length(dat(:,1))
    num_simple_dsb(k) = sum(dat(k,:) == 1); 
    num_complex_dsb(k) = sum(dat(k,:) == 2); 
    num_completed(k) = sum(dat(k,:) == 6) + sum(isnan(dat(k,:)));
end

time = time'/60;% in minutes
plot(time,num_simple_dsb,'b-'); 
hold on;
plot(time,num_complex_dsb,'r-'); 
plot(time,num_completed,'k-'); 
xlabel('Time / mins'); 
ylabel('Quantity'); 
legend('Simple breaks','Complex breaks', 'completed'); 
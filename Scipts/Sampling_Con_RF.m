
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%   This script generates a sample of extreme events of NTR peak and RF peak combinations where RF is
%   over the selecetd threshold. 
%
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear 
clc

load('***Path***\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat');
load('***Path***\POT_Events_for_each_RF_Accumulation_time.mat');


j=18; % 18 hr accumulation is considered
%% Sampling for the TC events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

POT_RF = events(j).POT;

%%
Maximum_NTR_events = struct;
Time_NTR= POT_NTR_and_NTR_timeseries.Time;
NTR = POT_NTR_and_NTR_timeseries.NTR;

%%
for i = 1:length(POT_RF)
    Time_RF = POT_RF(i,1);
    B=datevec(Time_RF); % Need to proceed to the datevec and again date num cz some issue
    Time_RF=datenum(B);
    [index_1] = find(Time_NTR(:,1)==Time_RF);
    [Max_NTR, index] = max(NTR(index_1-36:index_1+36)); % taking the 3 dday window
    index = index+index_1-35; % need to correct the index since 'index' was from the selected vector

          
    Maximum_NTR_events(j).Accumulation(i).Time_RF = Time_RF;
    Maximum_NTR_events(j).Accumulation(i).POT_RF = POT_RF(i,2);


    Maximum_NTR_events(j).Accumulation(i).Time_NTR = Time_NTR(index,1);
    Maximum_NTR_events(j).Accumulation(i).Max_NTR = Max_NTR;
end

%%
Maximum_NTR_events_of_POT_RF = Maximum_NTR_events(j).Accumulation;
save('Maximum_NTR_events_for_each_POT_RF_for_18hr_RF_acc',"Maximum_NTR_events_of_POT_RF")

%%


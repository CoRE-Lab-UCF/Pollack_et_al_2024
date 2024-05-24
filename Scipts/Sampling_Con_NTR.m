%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%   This script generates a sample of extreme events of NTR peak and RF peak combinations where NTR is
%   over the selecetd threshold. 
%
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

clear 
clc


load('***Path***\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat');
load('***Path***\Hourly_accumulation_Bias_Corrected_RF_data_Philli_airport.mat');

POT_NTR = POT_NTR_and_NTR_timeseries.POT;

%%
Maximum_RF_events = struct;
for j=1:48
    for i = 1:length(POT_NTR)
        Time_NTR = POT_NTR(i,1);
        [index_1] = find (Data(:,1)== Time_NTR);
        [Max_RF, index_2] = max(Data(index_1-36:index_1+36,j+1)); % taking the 3 dday window
        index_2 = index_2+index_1-35; % correcting the index to to find the row in 'Data'
    

        Maximum_RF_events(j).Accumulation(i).Time_NTR = Time_NTR;
        Maximum_RF_events(j).Accumulation(i).POT_NTR = POT_NTR(i,2);
              
        Maximum_RF_events(j).Accumulation(i).Time_RF = Data(index_2,1);
        Maximum_RF_events(j).Accumulation(i).Max_Rainfall = Max_RF;

    end
end
%%
save('Maximum_RF_events_for_each_POT_NTR',"Maximum_RF_events")

%%


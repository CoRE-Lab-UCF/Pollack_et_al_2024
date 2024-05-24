%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script creates an array of hourly accumulated RF for 1-48 hours
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc
close all

Pr = load("Bias Corrected MS at Philli airport_With_AORC_of_GC.mat");
Pr=Pr.Bias_corrected_MS(:,2);


%% Creating the date vector
t1 = datetime(1901,01,01,01,00,00,'Format','yyyy MM dd mm');
t2 = datetime(2021,12,31,23,00,00,'Format','yyyy MM dd mm');
t=t1:hours(1):t2;
t=t';
time=datenum(t);

%% Creating accumulation data for 1 to 48 hours

Data = NaN(length(Pr),49);% crating a dummy vetor
for acc = 1:48
    nos = length(Pr)-acc+1;
    
    for i = 1 : nos
        Data(i+acc-1,acc+1)= sum(Pr(i:i+acc-1));
    end

end
%% Saving the data
Data(:,1)= time;
save('Hourly_accumulation_Bias_Corrected_RF_data_Philli_airport','Data')






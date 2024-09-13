%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script does the bias-correction for the Measured data set based on basin-averaged AORC data
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% When calculating the threshold for the modeled data set matching the number of Dry hours with the measured data set, the hours with trace precipitation may needed to be accounted. 

clear
close all
clc

%Loading Dta
Pr_Measured = load("**Path**\Phlli_airport_RF_data_1901_to_2021.mat"); % Loading the measured data file [time, hourly RF]
Pr_Measured=Pr_Measured.Philli_airport_RF_data_1901_to_2021;
AORC=load("**Path**\AORC_9_grid_loc_Gloucester.mat");% Loading the measured data file [time, hourly RF(grid1), hourly RF(grid2),...hourly RF(grid9)]

%Taking the basin average rainfall
Pr_AORC =mean(AORC.prcpdata(:,2:end),2);
Time_AORC=AORC.prcpdata(:,1);



%% Creating the date vector
t1 = datetime(1979,02,01,00,00,00,'Format','yyyy MM dd HH mm'); %% insert the starting date manually / NOTE THAT AORC STATRTS WITH 02ND OF JAN
t2 = datetime(2021,12,31,23,00,00,'Format','yyyy MM dd HH mm'); %% Enter the end date manually
 
t=t1:hours(1):t2; %create the time series
t=t';


%% Finding the threshold
%finding the index of begining of AORC data
index=find(Pr_Measured(:,1)==datenum(t1));

MS_Sorted = sort(Pr_Measured(index:end,2),'ascend');
AORC_Sorted = sort(Pr_AORC,'ascend');

index2 = find(MS_Sorted > 0, 1); % find the location of the first element larger than 0
Threshold = AORC_Sorted (index2);

%% Creating very Small rainfall in to dry hours 

MOD_PR_AORC = Pr_AORC;
MOD_PR_AORC(MOD_PR_AORC < Threshold)= 0; % Make all the values below the threshold, zero before fitting the distribution

%% fitting to a gamma Distribution
C=[];
Pr_Measured_1=Pr_Measured(:,2); 

gamma_MS=fitdist(Pr_Measured_1(index:end),'gamma');% We take only the overlapping period to fit the distribution
shp_par_MS = gamma_MS.a;
scl_par_MS = gamma_MS.b;

gamma_AORC=fitdist(MOD_PR_AORC,'gamma');
shp_par_AORC = gamma_AORC.a;
scl_par_AORC = gamma_AORC.b;

y = cdf('gamma',Pr_Measured_1,shp_par_MS,scl_par_MS); % Find the CDF of the measured data
Bias_corrected_MS(:,2) = gaminv(y,shp_par_AORC,scl_par_AORC); % Find the inverse values of the measured data

Bias_corrected_MS(:,1)=Pr_Measured(:,1);

% Saving the data file
save("Bias Corrected MS at Philli airport_With_AORC_of_GC.mat","Bias_corrected_MS");
% for St. Petersburg: save("Bias Corrected MS at St_Pete_With_AORC.mat","Bias_corrected_MS");


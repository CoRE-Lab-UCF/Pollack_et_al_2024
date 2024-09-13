%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script is the main program for calculating return periods and
%   scaling historical events to given return periods
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
clc


% ET_NTR = vector of POT Non-TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
% TC_NTR = vector of POT TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
% ET_RF = vector of POT Non-TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
% TC_RF = vector of POT TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
% Thres_NTR = single numeric vector of NTR threshold
% Thres_RF = single numeric vector of RF threshold
% n_years = Total number of years
% Q_RP = Vector of Return periods of interest
% l_b_NTR = Lover bound of desctritized NTR space for combining two populations;
% U_b_NTR = Upper bound of desctritized NTR space for combining two populations;
% l_b_RF = Lover bound of desctritized RF space for combining two populations;
% U_b_RF = Upper bound of desctritized RF space for combining two populations;
% RL_NTR, RL_RF = calculated combinde return level vectors for given return
% perios in the vector Q_RP


load("***Path***\ETC_events_conditioning_POT_NTR.mat");
ET_NTR=table2array(struct2table(ET_Cyclones(18).Event));
load("***Path***\TC_events_conditioning_POT_NTR.mat");
TC_NTR=table2array(struct2table(Tropical_Cyclones(18).Event));

clearvars Tropical_Cyclones ET_Cyclones

load("***Path***\ETC_events_conditioning_POT_RF.mat");
ET_RF=table2array(struct2table(ET_Cyclones(18).Event));
load("***Path***\TC_events_conditioning_POT_RF.mat");
TC_RF=table2array(struct2table(Tropical_Cyclones(18).Event));

Thres_NTR = 0.634; 
Thres_RF = 35.648; 
n_years = 114; 
Q_RP=[1,2,5,10,15,20,25,50,75,100,200,500];
l_b_NTR = Thres_NTR; 
U_b_NTR = 3.7; 
l_b_RF = Thres_RF; 
U_b_RF = 175;

% Running the following function to estimate the return periods
[RL_NTR, RL_RF]=Uni_Return_level_calc(ET_NTR,TC_NTR,ET_RF,TC_RF, Thres_NTR, Thres_RF,n_years,Q_RP,l_b_NTR,U_b_NTR,l_b_RF,U_b_RF);


clearvars -except RL_RF RL_NTR figure


%% Generating Scaled Rainfall Fields
Hour_no=28; % for the purpose of plotting
Acc=18; % The selected rainfall accumulation
Cycl_name ={'Sandy','Irene','Floyd'}; % Three TC events are considered
% Selected event = Irene (rank 2)
% For sandy; The NTR rank is 2 ( The second largest event based on peak NTR)

Rank=2;
Ex_events_TC_row = load("***Path***\TC_events_con_NTR_from_1991.mat");
Ex_events_TC= Ex_events_TC_row.Tropical_Cyclones(Acc).Event;
Ex_events_TC=table2array(struct2table(Ex_events_TC));
Ex_events_TC=sortrows(Ex_events_TC,2,"descend");
Max_NTR = Ex_events_TC(Rank,2);
Peak_hour =Ex_events_TC(Rank,1);
Max_RF = Ex_events_TC(Rank,4);
%Max_RF_hour = Ex_events_TC(Rank,3);

% loading NTR data
NTR=load("***Path***\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat");
% Loading the RF fields time series over the Gloucester City
RF_Field = load("***Path***\AORC_Gloucester_city_cliped_data.mat");
% importing basin average AORC Data 
AORC = load("***Path***\Hourly_accumulation_RF_data_Gloucester_City_AORC.mat");
AORC=AORC.Data;



% Finding the hourly NTR and rainfall field of 3 days of events
Time=NTR.POT_NTR_and_NTR_timeseries.Time; % Time of the peak NTR
ind_peak_ntr = find(Time==Peak_hour);
Day_3_ntr = NTR.POT_NTR_and_NTR_timeseries.NTR(ind_peak_ntr-36:ind_peak_ntr+35); % Time series of NTR for 3days around the peak
Day_3_ntr=[Time(ind_peak_ntr-36:ind_peak_ntr+35) Day_3_ntr];

D3_tide=NTR.POT_NTR_and_NTR_timeseries.Tide(ind_peak_ntr-36:ind_peak_ntr+35); % Time series of Tide for 3days around the peak
D3_WL=NTR.POT_NTR_and_NTR_timeseries.WL_raw(ind_peak_ntr-36:ind_peak_ntr+35); % Time series of WL for 3days around the peak

TTT=RF_Field.AORC_GC_rectang_clip.Time;
latt=RF_Field.AORC_GC_rectang_clip.lattitude;
lon=RF_Field.AORC_GC_rectang_clip.Longitude;
ind_peak_RF = find(TTT==Peak_hour);
RF_field_event = RF_Field.AORC_GC_rectang_clip.PRCP(:,:,ind_peak_RF-36:ind_peak_RF+35);

BA_AORC_hourly = [AORC(:,1) AORC(:,2)];
ind_peak_RF = find(BA_AORC_hourly(:,1)==Peak_hour);
BA_AORC_event_hourly = BA_AORC_hourly(ind_peak_RF-36:ind_peak_RF+35,1:2);


%% Scaling Rainfall fields and NTR time series of the event


for i=1:length(RL_NTR)    
    
    Scale_RF(i,1) = RL_RF(2,i)/Max_RF;
    Design_RF_field_scaled(i).event = RF_field_event(:,:,:).*Scale_RF(i,1);
    BA_RF_scaled(i).event = [BA_AORC_event_hourly(:,1), BA_AORC_event_hourly(:,2)*Scale_RF(i,1)];
    
    % Adjusting the NTR datum
    % find the lowest value
    min_NTR = min(Day_3_ntr(:,2));
    Day_3_ntr_datum_adj = Day_3_ntr(:,2)+abs(min_NTR);
    
    % calculating the scaling factor
    Scale_NTR(i,1) = (RL_NTR(2,i)+abs(min_NTR))/(Max_NTR+abs(min_NTR));
    Design_NTR_scaled(i).event = [Day_3_ntr(:,1) (Day_3_ntr_datum_adj.*Scale_NTR(i,1)-abs(min_NTR))];
  
end


Data_raw = [Day_3_ntr(:,1) Day_3_ntr(:,2)];
T=datetime(datevec(Data_raw(:,1)));


%% Fitting distribution to the WL moving average for the last 18.6 years 
st=18.6*365*24;
Mov_avg = NTR.POT_NTR_and_NTR_timeseries.MovA_WL(end-st:end);
Mov_avg(isnan(Mov_avg)==1, :) = [];
pd=fitdist(Mov_avg,'Normal');
x_values = -0.4:0.05:1;
y = pdf(pd,x_values);
figure
plot(x_values,y)
histfit(Mov_avg,40)
% 
WL_u=pd.mu+pd.sigma;
WL_l=pd.mu-pd.sigma;
WL_mean=pd.mu;


%% Creating Water levels combining Tide, MSL level and Scaled NTR

Design_WL=[];
Design_MSL=WL_mean; % use the 
TT=NTR.POT_NTR_and_NTR_timeseries.Time;
ind_peak = find(TT==Peak_hour);
Day_3_Tide = NTR.POT_NTR_and_NTR_timeseries.Tide(ind_peak-36:ind_peak+35);
Day_3_WL_detrend = NTR.POT_NTR_and_NTR_timeseries.WL_Det(ind_peak-15*24:ind_peak+15*24,2);
Design_event=struct;

for i=1:length(RL_NTR)
    wl=WL_mean; %mean(Day_3_WL_detrend); % or can use the WL mean (the mean of last 18 years)
    for j=1:length(RL_NTR)
        ntr= Design_NTR_scaled(j).event;
        Design_WL(:,1) = Day_3_Tide+wl+ntr(:,2);

        Design_rf= Design_RF_field_scaled(j).event;

        Design_event(j).event.Water_level = Design_WL;
        Design_event(j).event.Rainfall = Design_rf;
        Design_event(j).event.RP = RL_NTR(1,j);
        Design_event(j).event.NTR = ntr(:,2);
        Design_event(j).event.Tide = Day_3_Tide;
    end
    save((['Design_events_',Cycl_name{Rank}]),'Design_event');
end


    
%%  Defining the starting and ending time manually 
% These hours were checked and selected manually to optimize the
% computational time of the numerical modeling

st = 14; % Starting hour
endd = 50 ; % ending hour 


for i=1:length(RL_NTR)
    Design_events(i).WL = Data.Design_event(i).event.Water_level(st:endd);
    Design_events(i).Rainfall = Data.Design_event(i).event.Rainfall(:,:,st:endd);
    Design_events(i).RP = Data.Design_event(i).event.RP;
    Design_events(i).Tide = Data.Design_event(i).event.Tide(st:endd);
    Design_events(i).NTR = Data.Design_event(i).event.NTR(st:endd);
end


save((['Design_events_Cropped_',Cycl_name{Rank}]),'Design_events');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script is the main program for calculating return periods with bootstrapping and
%   scaling historical events to given return levels
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc
load("***Path***\ETC_events_conditioning_POT_NTR.mat");
ET_NTR=table2array(struct2table(ET_Cyclones(18).Event));
load("***Path***\TC_events_conditioning_POT_NTR.mat");
TC_NTR=table2array(struct2table(Tropical_Cyclones(18).Event));

clearvars Tropical_Cyclones ET_Cyclones

load("***Path***\ETC_events_conditioning_POT_RF.mat");
ET_RF=table2array(struct2table(ET_Cyclones(18).Event));
load("***Path***\TC_events_conditioning_POT_RF.mat");
TC_RF=table2array(struct2table(Tropical_Cyclones(18).Event));


% ET_NTR = vector of POT Non-TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
% TC_NTR = vector of POT TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
% ET_RF = vector of POT Non-TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
% TC_RF = vector of POT TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
% Thres_NTR = single numeric vector of NTR threshold
% Thres_RF = single numeric vector of RF threshold
% n_years = Total number of years
% Q_RP = Vector of Return periods of interest
% l_b_NTR = Lover bound of discretized NTR space for combining two populations;
% U_b_NTR = Upper bound of discretized NTR space for combining two populations;
% l_b_RF = Lover bound of discretized RF space for combining two populations;
% U_b_RF = Upper bound of discretized RF space for combining two populations;
% RL_NTR, RL_RF = calculated combined return level vectors for a given return
% perios in the vector Q_RP
% n_sim = number of bootstrapped samples needed

[RL_NTR_com,RL_RF_com, RL_TC_NTR, RL_ETC_NTR, RL_TC_RF, RL_ETC_RF,RP_NTR,RP_RF]=Uni_Return_level_calc_with_boostrap_new(n_sim,ET_NTR,TC_NTR,ET_RF,TC_RF, Thres_NTR, Thres_RF,n_years,Q_RP,l_b_NTR,U_b_NTR,l_b_RF,U_b_RF);



% making the infinite  numers finite numbers

RL_NTR_com(RL_NTR_com==inf)=10^30;
RL_RF_com(RL_RF_com==inf)=10^30;

% RL_NTR_com_I=RL_NTR_com;
% RL_RF_com_I=RL_RF_com;

for m=1:n_sim
    [C, i]=unique(RL_NTR_com(m+1,:),'last');
    RL_NTR_com_I(m,:)=interp1(C,RL_NTR_com(1,i),RP_NTR);

    [D, j]=unique(RL_RF_com(m+1,:),'last');
    RL_RF_com_I(m,:)=interp1(D,RL_RF_com(1,j),RP_RF);
end

%%


NTR_Com_50=prctile(RL_NTR_com_I(2:end,:),50);
NTR_TC_50=prctile(RL_TC_NTR,50);
NTR_ETC_50=prctile(RL_ETC_NTR(2:end,:),50);

RF_Com_50=prctile(RL_RF_com_I(2:end,:),50);
RF_TC_50=prctile(RL_TC_RF,50);
RF_ETC_50=prctile(RL_ETC_RF,50);

NTR_Com_95=prctile(RL_NTR_com_I(2:end,:),95);
NTR_TC_95=prctile(RL_TC_NTR,95);
NTR_ETC_95=prctile(RL_ETC_NTR,95);

RF_Com_95=prctile(RL_RF_com_I(2:end,:),95);
RF_TC_95=prctile(RL_TC_RF,95);
RF_ETC_95=prctile(RL_ETC_RF,95);

NTR_Com_05=prctile(RL_NTR_com_I(2:end,:),5);
NTR_TC_05=prctile(RL_TC_NTR,5);
NTR_ETC_05=prctile(RL_ETC_NTR,5);

RF_Com_05=prctile(RL_RF_com_I(2:end,:),5);
RF_TC_05=prctile(RL_TC_RF,5);
RF_ETC_05=prctile(RL_ETC_RF,5);


% Calculating the RLs of NTR for confidence intervals
RL_NTR_Q_95 = interp1(RP_NTR1,NTR_Com_95,Q_RP);
RL_NTR_Q_05 = interp1(RP_NTR2,NTR_Com_05,Q_RP);

% Calculating the RLs of RF for confidence intervals
RL_RF_Q_95 = interp1(RP_RF1,RF_Com_95,Q_RP);
RL_RF_Q_05 = interp1(RP_RF2,RF_Com_05,Q_RP);




%% Creating Scaled Rainfall Field
Hour_no=28;
Acc=18; % The selected rainfall accumulation
Cycl_name ={'Sandy','Irene','Floyd'};
% Selected event = Irene
% For sandy; The NTR rank is 2 ( The largest events based on peak NTR)

Rank=2;
Ex_events_TC_row = load("***Path***\TC_events_con_NTR_from_1991.mat");
Ex_events_TC= Ex_events_TC_row.Tropical_Cyclones(Acc).Event;
Ex_events_TC=table2array(struct2table(Ex_events_TC));
Ex_events_TC=sortrows(Ex_events_TC,2,"descend");
Max_NTR = Ex_events_TC(Rank,2);
Peak_hour =Ex_events_TC(Rank,1);
Max_RF = Ex_events_TC(Rank,4);
%Max_RF_hour = Ex_events_TC(Rank,3);

% loading NTR and RF data files
NTR=load("***Path***\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat");
RF_Field = load("***Path***\AORC_Gloucester_city_cliped_data.mat");
% importing basin average AORC Data 
AORC = load("***Path***\Hourly_accumulation_RF_data_Gloucester_City_AORC.mat");
AORC=AORC.Data;



% Find the hourly NTR and rainfall field of 3 days of events

Time=NTR.POT_NTR_and_NTR_timeseries.Time;
ind_peak_ntr = find(Time==Peak_hour);
Day_3_ntr = NTR.POT_NTR_and_NTR_timeseries.NTR(ind_peak_ntr-36:ind_peak_ntr+35);
Day_3_ntr=[Time(ind_peak_ntr-36:ind_peak_ntr+35) Day_3_ntr ];

D3_tide=NTR.POT_NTR_and_NTR_timeseries.Tide(ind_peak_ntr-36:ind_peak_ntr+35);
D3_WL=NTR.POT_NTR_and_NTR_timeseries.WL_raw(ind_peak_ntr-36:ind_peak_ntr+35);

TTT=RF_Field.AORC_GC_rectang_clip.Time;
latt=RF_Field.AORC_GC_rectang_clip.lattitude;
lon=RF_Field.AORC_GC_rectang_clip.Longitude;
ind_peak_RF = find(TTT==Peak_hour);
RF_field_event = RF_Field.AORC_GC_rectang_clip.PRCP(:,:,ind_peak_RF-36:ind_peak_RF+35);

BA_AORC_hourly = [AORC(:,1) AORC(:,2)];
ind_peak_RF = find(BA_AORC_hourly(:,1)==Peak_hour);
BA_AORC_event_hourly = BA_AORC_hourly(ind_peak_RF-36:ind_peak_RF+35,1:2);



%% Scaling RF and NTR of the event (95% upper margin)


% for 95th Percentile
for i=1:length(RL_RF_Q_95)    
    
    Scale_RF(i,1) = RL_RF_Q_95(1,i)/Max_RF;
    Design_RF_field_scaled(i).event = RF_field_event(:,:,:).*Scale_RF(i,1);
    BA_RF_scaled(i).event = [BA_AORC_event_hourly(:,1), BA_AORC_event_hourly(:,2)*Scale_RF(i,1)];
    
    % Adjusting the NTR datum
    % find the lowest value
    min_NTR = min(Day_3_ntr(:,2));
    Day_3_ntr_datum_adj = Day_3_ntr(:,2)+abs(min_NTR);
    
    % calculating the scaling factor
    Scale_NTR(i,1) = (RL_NTR_Q_95(1,i)+abs(min_NTR))/(Max_NTR+abs(min_NTR));
    Design_NTR_scaled(i).event = [Day_3_ntr(:,1) (Day_3_ntr_datum_adj.*Scale_NTR(i,1)-abs(min_NTR))];
  
end


Data_raw = [Day_3_ntr(:,1) Day_3_ntr(:,2)];
T=datetime(datevec(Data_raw(:,1)));


%% Fitting distribution to the WL moving average for the last 18 years 
st=18.61*365*24;
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


%% Creating Water levels combining Tide, MA_water level and Scaled NTR
i=1;
Design_WL=[];
Design_MSL=WL_mean; % use the 
TT=NTR.POT_NTR_and_NTR_timeseries.Time;
ind_peak = find(TT==Peak_hour);
Day_3_Tide = NTR.POT_NTR_and_NTR_timeseries.Tide(ind_peak-36:ind_peak+35);
Design_event=struct;

for i=1:length(RL_RF_Q_95)
    wl=Design_MSL; 
    for j=1:length(RL_RF_Q_95)
        ntr= Design_NTR_scaled(j).event;
        Design_WL(:,1) = Day_3_Tide+wl+ntr(:,2);

        Design_rf= Design_RF_field_scaled(j).event;

        Design_event(j).event.Water_level = Design_WL;
        Design_event(j).event.Rainfall = Design_rf;
        Design_event(j).event.RP = Q_RP(1,j);
        Design_event(j).event.NTR = ntr(:,2);
        Design_event(j).event.Tide = Day_3_Tide;
    end
    save((['Design_events_95_CI_U_',Cycl_name{Rank}]),'Design_event');
end



    %% Select the starting and ending time manually
    st = 14;
    endd = 50 ; %these hours were checked and selected manually
for i=1:length(RL_NTR_Q_95)

    Design_events(i).WL = Data.Design_event(i).event.Water_level(st:endd);
    Design_events(i).Rainfall = Data.Design_event(i).event.Rainfall(:,:,st:endd);
    Design_events(i).RP = Data.Design_event(i).event.RP;
    Design_events(i).Tide = Data.Design_event(i).event.Tide(st:endd);
    Design_events(i).NTR = Data.Design_event(i).event.NTR(st:endd);


end


save((['Design_events_Cropped_95%_CI_U_',Cycl_name{Rank}]),'Design_events');



%% Scaling RF and NTR of the event (95% Lower margin)

% Creating Scaled Rainfall Field
Hour_no=28;
Acc=18; % The selected rainfall accumulation
Cycl_name ={'Sandy','Irene','Floyd'};
% Selected event = Irene
% For sandy; The NTR rank is 2 ( The largest events based on peak NTR)

Rank=2;
Ex_events_TC_row = load("***Path***\TC_events_con_NTR_from_1991.mat");
Ex_events_TC= Ex_events_TC_row.Tropical_Cyclones(Acc).Event;
Ex_events_TC=table2array(struct2table(Ex_events_TC));
Ex_events_TC=sortrows(Ex_events_TC,2,"descend");
Max_NTR = Ex_events_TC(Rank,2);
Peak_hour =Ex_events_TC(Rank,1);
Max_RF = Ex_events_TC(Rank,4);
%Max_RF_hour = Ex_events_TC(Rank,3);

% loading NTR and RF data files
NTR=load("***Path***\POT_NTR_and_NTR_timeseries_for_Philli_airport.mat");
RF_Field = load("***Path***\AORC_Gloucester_city_cliped_data.mat");
% importing basin average AORC Data 
AORC = load("***Path***\Hourly_accumulation_RF_data_Gloucester_City_AORC.mat");
AORC=AORC.Data;



% Find the hourly NTR and rainfall field of 3 days of events

Time=NTR.POT_NTR_and_NTR_timeseries.Time;
ind_peak_ntr = find(Time==Peak_hour);
Day_3_ntr = NTR.POT_NTR_and_NTR_timeseries.NTR(ind_peak_ntr-36:ind_peak_ntr+35);
Day_3_ntr=[Time(ind_peak_ntr-36:ind_peak_ntr+35) Day_3_ntr ];

D3_tide=NTR.POT_NTR_and_NTR_timeseries.Tide(ind_peak_ntr-36:ind_peak_ntr+35);
D3_WL=NTR.POT_NTR_and_NTR_timeseries.WL_raw(ind_peak_ntr-36:ind_peak_ntr+35);

TTT=RF_Field.AORC_GC_rectang_clip.Time;
latt=RF_Field.AORC_GC_rectang_clip.lattitude;
lon=RF_Field.AORC_GC_rectang_clip.Longitude;
ind_peak_RF = find(TTT==Peak_hour);
RF_field_event = RF_Field.AORC_GC_rectang_clip.PRCP(:,:,ind_peak_RF-36:ind_peak_RF+35);

BA_AORC_hourly = [AORC(:,1) AORC(:,2)];
ind_peak_RF = find(BA_AORC_hourly(:,1)==Peak_hour);
BA_AORC_event_hourly = BA_AORC_hourly(ind_peak_RF-36:ind_peak_RF+35,1:2);


% for 5th Percentile
for i=1:length(RL_RF_Q_05)    
    
    Scale_RF(i,1) = RL_RF_Q_05(1,i)/Max_RF;
    Design_RF_field_scaled(i).event = RF_field_event(:,:,:).*Scale_RF(i,1);
    BA_RF_scaled(i).event = [BA_AORC_event_hourly(:,1), BA_AORC_event_hourly(:,2)*Scale_RF(i,1)];
    
    % Adjusting the NTR datum
    % find the lowest value
    min_NTR = min(Day_3_ntr(:,2));
    Day_3_ntr_datum_adj = Day_3_ntr(:,2)+abs(min_NTR);
    
    % calculating the scaling factor
    Scale_NTR(i,1) = (RL_NTR_Q_05(1,i)+abs(min_NTR))/(Max_NTR+abs(min_NTR));
    Design_NTR_scaled(i).event = [Day_3_ntr(:,1) (Day_3_ntr_datum_adj.*Scale_NTR(i,1)-abs(min_NTR))];
  
end

clearvars RL_NTR_Q_95 RL_RF_Q_95


Data_raw = [Day_3_ntr(:,1) Day_3_ntr(:,2)];
T=datetime(datevec(Data_raw(:,1)));



%% Creating Water levels combining Tide, MA_water level and Scaled NTR
i=1;
Design_WL=[];
Design_MSL=WL_mean; % use the 
TT=NTR.POT_NTR_and_NTR_timeseries.Time;
ind_peak = find(TT==Peak_hour);
Day_3_Tide = NTR.POT_NTR_and_NTR_timeseries.Tide(ind_peak-36:ind_peak+35);
Design_event=struct;

for i=1:length(RL_RF_Q_05)
    wl=Design_MSL; 
    for j=1:length(RL_RF_Q_05)
        ntr= Design_NTR_scaled(j).event;
        Design_WL(:,1) = Day_3_Tide+wl+ntr(:,2);

        Design_rf= Design_RF_field_scaled(j).event;

        Design_event(j).event.Water_level = Design_WL;
        Design_event(j).event.Rainfall = Design_rf;
        Design_event(j).event.RP = Q_RP(1,j);
        Design_event(j).event.NTR = ntr(:,2);
        Design_event(j).event.Tide = Day_3_Tide;
    end
    save((['Design_events_95_CI_L_',Cycl_name{Rank}]),'Design_event');
end


WL_raw_ind= find(NTR.POT_NTR_and_NTR_timeseries.Time==Peak_hour);
WL_raw = NTR.POT_NTR_and_NTR_timeseries.WL_raw(WL_raw_ind-36:WL_raw_ind+35);
MA_wl = NTR.POT_NTR_and_NTR_timeseries.MovA_WL(WL_raw_ind-36:WL_raw_ind+35);


 %% Select the starting and ending time manually
    st = 14;
    endd = 50 ; %these hours were checked and selected manually
for i=1:length(RL_NTR_Q_05)

    Design_events(i).WL = Data.Design_event(i).event.Water_level(st:endd);
    Design_events(i).Rainfall = Data.Design_event(i).event.Rainfall(:,:,st:endd);
    Design_events(i).RP = Data.Design_event(i).event.RP;
    Design_events(i).Tide = Data.Design_event(i).event.Tide(st:endd);
    Design_events(i).NTR = Data.Design_event(i).event.NTR(st:endd);


end

save((['Design_events_Cropped_95%_CI_L_',Cycl_name{Rank}]),'Design_events');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script creates POT NTR events using Water level data
%   
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Load The data
clear
clc
close all

% Loading the data file
data = load('Water_Level_Time_Series_1901_to_2021.mat');

% 2. Detrend
Data=data.Data;

plot(Data(:,1),Data(:,2))
hold on

% detrending Using moving 30 DAY Mean
A=movmean(Data(:,2),720,"omitnan");


WL_Detrended(:,2)=Data(:,2)-A;
WL_Detrended(:,1)=Data(:,1);

plot(WL_Detrended(:,1),WL_Detrended(:,2))

%% remmoving the data chuncks where no sufficiant data are available to calculate moving average
WL_Detrended_P =WL_Detrended(:,2);
% Process1
for i = 1:(length(WL_Detrended(:,2))-720)
    D = sum(isnan(WL_Detrended(i:i+720,2)));
    if D>120
        WL_Detrended_P(i+360)=NaN;
    else
    end
    
end

% Remove 1st and last 15 days of data
WL_Detrended_P(1:360)=NaN;
WL_Detrended_P((length(WL_Detrended_P)-360):end)=NaN;

%%
WL_Detrended(:,2)=WL_Detrended_P;

dates=datevec(WL_Detrended(:,1));
uyears=unique(dates(:,1));

% The years with many dta gaps are are manually removed from tidal analysis
N_loc=find(dates(:,1)==1903);
WL_Detrended(N_loc,2)=NaN;

N_loc=find(dates(:,1)==1921);
WL_Detrended(N_loc,2)=NaN;

N_loc=find(dates(:,1)==1922);
WL_Detrended(N_loc,2)=NaN;

N_loc=find(dates(:,1)==1959);
WL_Detrended(N_loc,2)=NaN;


plot(WL_Detrended(:,1),WL_Detrended(:,2));


%%  Generating the tidal constituents

lat=39.89; % lattutude of the site
Tide=[];
T_time=[];
k=1;
Co_ef=[];

for i = 1:length(uyears)
    index = find(dates(:,1)==(uyears(i)));
    sl=WL_Detrended(index,2);
    t=WL_Detrended(index,1);

    if  uyears(i)==1903 || uyears(i)==1921||uyears(i)==1922 ||uyears(i)==1959 % Avoiding the years 
        tide_year=NaN(length(t),1);
        Tide =cat(1,Tide,tide_year);
        T_time=cat(1,T_time,t);
    else
        coef = ut_solv(t,sl,[],lat,'auto',RunTimeDisp='nnn'); 
        tide_year = ut_reconstr(t,coef);
        Tide =cat(1,Tide,tide_year);
        T_time=cat(1,T_time,t);
        CO_ef(:,i)=coef.A(1:30,1);   
        clearvars t_raw coef sl_raw t


    end

end 




%% Creating the NTR time series
Time=Data(:,1);

% Plotting a sample for a selected time duration
NTR = WL_Detrended(:,2)-Tide;
figure
plot(Time(212000:212048),WL_Detrended(212000:212048,2));
hold on
plot(Time(212000:212048),Tide(212000:212048,1));
hold on
plot(Time(212000:212048),NTR(212000:212048))

X=movmean(NTR,200);
plot(X)

%% Saving the data file
POT_NTR_and_NTR_timeseries=struct;
POT_NTR_and_NTR_timeseries.Time=Time;
POT_NTR_and_NTR_timeseries.Tide=Tide;
POT_NTR_and_NTR_timeseries.WL_Det=WL_Detrended;


%%  Events over threshold
% decluster time

ts=[Time NTR];

th = 0.634; % The threshold was dentified repeating the POT event generation 
%for different thresholds as we get 5 events per year

% events over threshold
EOT= ts(:,2);
EOT(EOT< th)= NaN;
EOT= [ts(:,1) EOT];
EOT(isnan(EOT(:,2)),:)= [];

% Plot the NTR and values over Threshold

hh= figure;
set(hh,'units','centimeters','Position',[0.2 1.5 28 12],...
    'color','w');
hold all
% h1= plot(ts(:,1),ts(:,2),'.-','LineStyle', 'none');
h1= plot(ts(:,1),ts(:,2),'.-');

ylabel('NTR (m)');
set(gca,'FontSize',11,'FontName','Times New Roman');
grid minor
ax = gca;
ax.XTickLabelRotation= 0;

title('Philadelphia Tide Gauge')

axs= axis;
hold all; 
% plot(EOT(:,1),EOT(:,2),',-k','LineStyle', 'none')
plot(EOT(:,1),EOT(:,2),'.k')

% Declustering 
dec_tim = 2.5; % decluster time in days
xnan= sum(isnan(EOT(:,2)));
POT= nan(length(EOT),2);

while xnan< size(EOT,1)
    
    for i= 1: size(POT,2)
        
        [~,fmax]= max(EOT(:,2));
        
        POT(fmax,:)= EOT(fmax,:);
        hold all; plot(EOT(fmax,1),EOT(fmax,2),'LineWidth',1);
        
        dec_wind= find(EOT(:,1)>= EOT(fmax,1)-dec_tim & EOT(:,1)<= EOT(fmax,1)+dec_tim);
        
        hold all; hi= plot(EOT(dec_wind,1),EOT(dec_wind,2),...
            '.-','Color','r');
        
        EOT(dec_wind,2)= nan;
        
        xnan= sum(isnan(EOT(:,2)));
    end
    
end


%% Saving NTR file

POT_NTR_and_NTR_timeseries.NTR=NTR;
POT_NTR_and_NTR_timeseries.POT = POT;
POT_NTR_and_NTR_timeseries.WL_raw = Data(:,2);
POT_NTR_and_NTR_timeseries.Tide = Tide;
POT_NTR_and_NTR_timeseries.MovA_WL = A;
POT_NTR_and_NTR_timeseries.Threshold = th;

save('POT_NTR_and_NTR_timeseries_for_Philli_airport.mat','POT_NTR_and_NTR_timeseries')

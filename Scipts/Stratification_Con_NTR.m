
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%   This script stratifies the extreme events sample in to TC events and
%   non-TC events when conditioning on NTR
%
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 




clear
clc
load('***Path***\Maximum_RF_events_for_each_POT_NTR.mat');
load('***Path***\Cyclone_Track_data_from_1850_new');


%% 
count_TC= zeros(48,1);
count_ETC =zeros(48,1);
distance = km2deg(350); % in degrees 350 km 
Tropical_Cyclones = struct;
ET_Cyclones = struct;

%% One event will be manually transfered in to the TC List due to its closeness to the search radius (date: 2009/09/19)

for i = 1:48
    X=Maximum_RF_events(1).Accumulation;
    X=struct2table(X);
    X=X.POT_NTR;
    leng=length(X);
    m=1;n=1;
    for j =1:leng
        Time_ntr = Maximum_RF_events(i).Accumulation(j).Time_NTR;
        NTR = Maximum_RF_events(i).Accumulation(j).POT_NTR;
        Time_RF = Maximum_RF_events(i).Accumulation(j).Time_RF;
        RF = Maximum_RF_events(i).Accumulation(j).Max_Rainfall;
    
    
        T_elements = datevec(Time_ntr); % Resolutions of data are different 6 hr and 1 hr
        hr = T_elements(1,4);
        hr_new = floor(hr/6)*6; 
        T_elements(1,4)= hr_new;
        Time_NTR_mod = datenum(T_elements);
        strt_hr = Time_NTR_mod-2;
        end_hr = Time_NTR_mod+1.75;
        
        index = find(Cyclone_track_data(:,1)>=strt_hr & Cyclone_track_data(:,1)<=end_hr);
    
        if length(index)>0
            if Time_ntr>731842 && Time_ntr<731844 % One event is manually added due to its closeness to the search radius
               Tropical_Cyclones(i).Event(m).Time_NTR = Time_ntr;
               Tropical_Cyclones(i).Event(m).NTR = NTR;
               Tropical_Cyclones(i).Event(m).Time_RF = Time_RF;
               Tropical_Cyclones(i).Event(m).RF = RF;
               m=m+1;
               count_TC(i,1) = count_TC(i,1)+1;
            
            elseif min(Cyclone_track_data(index,5))<= distance
               Tropical_Cyclones(i).Event(m).Time_NTR = Time_ntr;
               Tropical_Cyclones(i).Event(m).NTR = NTR;
               Tropical_Cyclones(i).Event(m).Time_RF = Time_RF;
               Tropical_Cyclones(i).Event(m).RF = RF;
               m=m+1;
               count_TC(i,1) = count_TC(i,1)+1;
            else
               ET_Cyclones(i).Event(n).Time_NTR = Time_ntr;
               ET_Cyclones(i).Event(n).NTR = NTR;
               ET_Cyclones(i).Event(n).Time_RF = Time_RF;
               ET_Cyclones(i).Event(n).RF = RF;
               n=n+1;
               count_ETC(i,1) = count_ETC(i,1)+1;
            end
        else
           ET_Cyclones(i).Event(n).Time_NTR = Time_ntr;
           ET_Cyclones(i).Event(n).NTR = NTR;
           ET_Cyclones(i).Event(n).Time_RF = Time_RF;
           ET_Cyclones(i).Event(n).RF = RF;
           n=n+1;
           count_ETC(i,1) = count_ETC(i,1)+1;

        end
        
    end
end
%%
save('TC_events_conditioning_POT_NTR',"Tropical_Cyclones")
save('ET_events_conditioning_POT_NTR',"ET_Cyclones")


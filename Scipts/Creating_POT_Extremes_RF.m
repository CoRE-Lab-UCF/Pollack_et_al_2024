%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This Script creates POT RF events using bias Corrected RF data level data
%   
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear; close all; clc
load("Hourly_accumulation_Bias_Corrected_RF_data_Philli_airport.mat") % load the Bias Corrected Data set
RF_Data = Data(1:end,:); % Taking the Rainfall Data only Belongs to the Tide guage data period


%%  Events over threshold

ts = RF_Data;
events= struct;


th = ones(1,48)*11; % initial Threshold values
for ii = 1:48
    
    POT=nan(10000,2); %This is a dummy vector just to recognize for the loop
    while length(POT) > 580 % Should be changed accordingly
        % events over threshold events
        EOT= ts(:,ii+1);
        EOT(EOT< th(ii))= NaN;
        EOT= [ts(:,1) EOT];
        EOT(isnan(EOT(:,2)),:)= [];
        
        % Plot the RF values over Threshold
        

        %%%%%%%%%%%% Declustering %%%%%%%%%
        % decluster time
        dec_tim= 2.5; % in days
        xnan= sum(isnan(EOT(:,2)));
        POT= nan(length(EOT),2);
        
        while xnan < size(EOT,1)
            
            for i= 1: size(POT,2)
                
                [~,fmax]= max(EOT(:,2));
                
                POT(fmax,:)= EOT(fmax,:);
                %hold all; plot(EOT(fmax,1),EOT(fmax,2),'or','LineWidth',2);
                
                dec_wind= find(EOT(:,1)>= EOT(fmax,1)-dec_tim & EOT(:,1)<= EOT(fmax,1)+dec_tim);
                
                %hold all; hi= plot(EOT(dec_wind,1),EOT(dec_wind,2),...
                %    '.-','Color','r');
                
                EOT(dec_wind,2)= nan;
                
                xnan= sum(isnan(EOT(:,2)));
            end
        
        end
        
        POT(isnan(POT(:,2)),:)= [];
        
        
        th(1,ii)=th(1,ii)+0.005; % the 0.005 was just a user defined value to go to the next threshold 
        
    end
    events(ii).POT=POT;
    events(ii).Threshold=th(1,ii)-0.005;

    th(ii+1) = th(1,ii);
    th(1,ii)

end
%% Saving the file

save("POT_Events_for_each_RF_Accumulation_time","events");



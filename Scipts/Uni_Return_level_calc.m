%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This function calculates the return periods using tratified POT samples
%   
%   IMPORTANT: The paths included in the script are according to the
%   author's directory. Please change them accordingly



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% ET_NTR = vector of POT Non-TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
% TC_NTR = vector of POT TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
% ET_RF = vector of POT Non-TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
% TC_RF = vector of POT TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
% Thres_NTR = single numeric vector of NTR threshold
% Thres_RF = single numeric vector of RF threshold
% n_years = Total number of years
% RP = Vector of Return periods of interest
% l_b_NTR = Lover bound of desctritized NTR space for combining two populations;
% U_b_NTR = Upper bound of desctritized NTR space for combining two populations;
% l_b_RF = Lover bound of desctritized RF space for combining two populations;
% U_b_RF = Upper bound of desctritized RF space for combining two populations;

function [RL_NTR, RL_RF,Data_TC_ETC_NTR,Data_TC_ETC_RF]=Uni_Return_level_calc(ET_NTR,TC_NTR,ET_RF,TC_RF, Thres_NTR, Thres_RF,n_years,Q_RP,l_b_NTR,U_b_NTR,l_b_RF,U_b_RF)
    %% Plotting the data
    % ploting data for NTR
    figure
    scatter(datetime(datevec(ET_NTR(:,1))),ET_NTR(:,2)); hold on;
    scatter(datetime(datevec(TC_NTR(:,1))),TC_NTR(:,2));ylabel('NTR (m)');xlabel('Time')
    
    % Plotting for Rainfall
    figure
    scatter(datetime(datevec(ET_RF(:,3))),ET_RF(:,4)); hold on;
    scatter(datetime(datevec(TC_RF(:,3))),TC_RF(:,4));ylabel('18 hours accumulated rainfall (m)');xlabel('Time')
    
    n=n_years;
    %% Fitting GPDs to NTR
    % For TC events
    [x_TC_NTR,is]= sort(TC_NTR(:,2));
    x_TC_time = TC_NTR(is,1);
    x_TC_NTR= x_TC_NTR -Thres_NTR;
    x_TC_NTR(x_TC_NTR<=0) = 0.001;
    [parm_TC_NTR,bounds_TC_NTR] = gpfit(x_TC_NTR,0.05);
    
    
    %For ETC events
    [x_ETC_NTR,is]= sort(ET_NTR(:,2));
    x_ETC_time = ET_NTR(is,1);
    x_ETC_NTR= x_ETC_NTR -Thres_NTR;
    x_ETC_NTR(x_ETC_NTR<=0) = 0.001;
    [parm_ETC_NTR,bounds_ETC_NTR] = gpfit(x_ETC_NTR,0.05);
    
    
    %check fitting
    lowerBnd = min([x_ETC_NTR;x_TC_NTR]);
    xmax= max([x_ETC_NTR;x_TC_NTR])*1.2;
    ygrid = linspace(lowerBnd,xmax,1000);
    figure
    plot(ygrid,gpcdf(ygrid,parm_TC_NTR(1),parm_TC_NTR(2),Thres_NTR),'-b','LineWidth',1);ylabel('Cumulative probability');xlabel('NTR (m)')
    hold on;
    plot(ygrid,gpcdf(ygrid,parm_ETC_NTR(1),parm_ETC_NTR(2),Thres_NTR),'-r','LineWidth',1);
    % Checking the empreical probabiliyt
    [F_TC,yi_TC]= ecdf(x_TC_NTR);
    [F_ETC,yi_ETC]= ecdf(x_ETC_NTR);
    stairs(yi_TC+Thres_NTR,F_TC,'+','LineWidth',1,'Color','b');
    stairs(yi_ETC+Thres_NTR,F_ETC,'+','LineWidth',1,'Color','r');
    legend('Fitted GPD: TC','Fitted GPD: Non-TC','Emperical P.: TC','Emperical P.: Non-TC')
    
    
    
    %%
    
    % Return level Calculation
    
    % Annual Exceedence prbability vector
    RP= (0.1:0.1:max(Q_RP));
    A_NEP = 1./RP;
    
    % Calculating the Inter Arrival Time (years)
    mu_TC = n/length(x_TC_NTR);
    mu_ETC = n/length(x_ETC_NTR);
    
    % Calculating Non_exceedence probabiliy [ No_Ex_prob = 1-(IAT*AEP)]
    Prob_TC= 1-mu_TC.*A_NEP;
    Prob_ETC= 1-mu_ETC.*A_NEP;

    
    % Return levels
    RL_TC= gpinv(Prob_TC,parm_TC_NTR(1),parm_TC_NTR(2),Thres_NTR);
    RL_ETC= gpinv(Prob_ETC,parm_ETC_NTR(1),parm_ETC_NTR(2),Thres_NTR);
    
    Data_TC_ETC_NTR(1,:)=RP;
    Data_TC_ETC_NTR(2,:)=RL_TC;
    Data_TC_ETC_NTR(3,:)=RL_ETC;
    
    % Combining Retuen periods
    NTR = l_b_NTR:0.01:U_b_NTR; % Descritized NTR vector

    ANEP_TC = 1-((1-gpcdf(NTR,parm_TC_NTR(1),parm_TC_NTR(2),Thres_NTR))./mu_TC); % Calculating ANEP
    ANEP_ETC = 1-((1-gpcdf(NTR,parm_ETC_NTR(1),parm_ETC_NTR(2),Thres_NTR))./mu_ETC);
    
    Com_AEP = 1-(ANEP_TC.*ANEP_ETC); % Combining AEP
    
    Com_RP = 1./Com_AEP; % Converting AEP in to Return Period
    
    
    % Empirical return period
    rank_TC = fliplr(1:1:length(x_TC_NTR))';
    rank_ETC = fliplr(1:1:length(x_ETC_NTR))';
    Weibull_TC    = (n+1)./rank_TC;
    Weibull_ETC    = (n+1)./rank_ETC;
    
    
    % Plot
    hh=figure;
    set(hh,'units','centimeters','Position',[3 1.5 20 12],'color','w');
    
    h1=plot(RP,RL_TC,"Color",[0.8500 0.3250 0.0980],'LineWidth',1); hold on;
    plot(Weibull_TC,x_TC_NTR+Thres_NTR,'x','LineWidth',1,"Color",[0.8500 0.3250 0.0980]);
    
    h2=plot(RP,RL_ETC,"Color",[0 0.4470 0.7410],'LineWidth',1); hold on
    plot(Weibull_ETC,x_ETC_NTR+Thres_NTR,'x','LineWidth',1,"Color",[0 0.4470 0.7410]);
    
    h3=plot(Com_RP,NTR,'k','LineWidth',3); hold on
    
    ylabel('NTR (m)'); xlabel('Return period (years)'); grid on;
    set(gca, 'XScale', 'log','XTick',[2 5 10 25 50 100 200],'FontSize',14,'FontName','Times New Roman');
    legend([h1 h2 h3],'TC Events','Non-TC Events','Combined','Location','northwest'); title("Return levels of NTR");
    clearvars mu RL 
    
    saveas(hh,"Return_Levels_NTR");
    
    % Calculating return levels for given return periods from combined line
    RL_NTR(1,:) = Q_RP;
    RL_NTR(2,:) = interp1(Com_RP,NTR,Q_RP,"linear");
    
    Data_TC_ETC_NTR(4,:)=interp1(Com_RP,NTR,RP,"linear");


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Fitting GPDs to RF
    % For TC events
    [x_TC_RF,is]= sort(TC_RF(:,4));
    x_TC_time = TC_RF(is,3);
    x_TC_RF= x_TC_RF -Thres_RF;
    x_TC_RF(x_TC_RF<=0) = 0.001;
    [parm_TC_RF,bounds_TC_RF] = gpfit(x_TC_RF,0.05);
    
    
    %For ETC events
    [x_ETC_RF,is]= sort(ET_RF(:,4));
    x_ETC_time = ET_RF(is,3);
    x_ETC_RF= x_ETC_RF -Thres_RF;
    x_ETC_RF(x_ETC_RF<=0) = 0.001;
    [parm_ETC_RF,bounds_ETC_RF] = gpfit(x_ETC_RF,0.05);
    
    
    %check fitting
    lowerBnd = min([x_ETC_RF;x_TC_RF]);
    xmax= max([x_ETC_RF;x_TC_RF])*1.2;
    ygrid = linspace(lowerBnd,xmax,1000);
    figure
    plot(ygrid,gpcdf(ygrid,parm_TC_RF(1),parm_TC_RF(2),Thres_RF),'-b','LineWidth',1);ylabel('Cumulative probability');xlabel('18 hr rainfall (mm)')
    hold on;
    plot(ygrid,gpcdf(ygrid,parm_ETC_RF(1),parm_ETC_RF(2),Thres_RF),'-r','LineWidth',1);
    % Checking the empreical probabiliyt
    [F_TC,yi_TC]= ecdf(x_TC_RF);
    [F_ETC,yi_ETC]= ecdf(x_ETC_RF);
    stairs(yi_TC+Thres_RF,F_TC,'+','LineWidth',1,'Color','b');
    stairs(yi_ETC+Thres_RF,F_ETC,'+','LineWidth',1,'Color','r');
    legend('Fitted GPD: TC','Fitted GPD: Non-TC','Emperical P.: TC','Emperical P.: Non-TC')
    
    
    
    %%
    
    % Return level Calculation
    
    % Annual Exceedence prbability vector
    RP= (0.1:0.1:max(Q_RP));
    A_NEP = 1./RP;
    
    % Calculating the Inter Arrival Time (years)
    mu_TC = n/length(x_TC_RF);
    mu_ETC = n/length(x_ETC_RF);
    
    % Calculating Non_exceedence probabiliy [ No_Ex_prob = 1-(IAT*AEP)]
    Prob_TC= 1-mu_TC.*A_NEP;
    Prob_ETC= 1-mu_ETC.*A_NEP;
    
    % Return levels
    RL_TC= gpinv(Prob_TC,parm_TC_RF(1),parm_TC_RF(2),Thres_RF);
    RL_ETC= gpinv(Prob_ETC,parm_ETC_RF(1),parm_ETC_RF(2),Thres_RF);
    
    Data_TC_ETC_RF(1,:)=RP;
    Data_TC_ETC_RF(2,:)=RL_TC;
    Data_TC_ETC_RF(3,:)=RL_ETC;
    
    
    % Combining Retuen periods
    RF = l_b_RF:0.01:U_b_RF; % Descritized RF vector
     
    ANEP_TC = 1-((1-gpcdf(RF,parm_TC_RF(1),parm_TC_RF(2),Thres_RF))./mu_TC); % Calculating AEP
    ANEP_ETC = 1-((1-gpcdf(RF,parm_ETC_RF(1),parm_ETC_RF(2),Thres_RF))./mu_ETC);
    Com_AEP = 1-(ANEP_TC.*ANEP_ETC); % Combining AEP
    
    Com_RP = 1./Com_AEP; % Converting AEP in to Return Period
    
    
    % Empirical return period
    rank_TC = fliplr(1:1:length(x_TC_RF))';
    rank_ETC = fliplr(1:1:length(x_ETC_RF))';
    
    Weibull_TC = (n+1)./rank_TC;
    Weibull_ETC = (n+1)./rank_ETC;
    
    
    % Plot
    hh1=figure;
    set(hh1,'units','centimeters','Position',[3 1.5 20 12],'color','w');
    
    h1=plot(RP,RL_TC,"Color",[0.8500 0.3250 0.0980],'LineWidth',1); hold on;
    plot(Weibull_TC,x_TC_RF+Thres_RF,'x','LineWidth',1,"Color",[0.8500 0.3250 0.0980]);
    
    h2=plot(RP,RL_ETC,"Color",[0 0.4470 0.7410],'LineWidth',1); hold on
    plot(Weibull_ETC,x_ETC_RF+Thres_RF,'x','LineWidth',1,"Color",[0 0.4470 0.7410]);
    
    h3=plot(Com_RP,RF,'k','LineWidth',3); hold on
    
    ylabel('18 hour rainfall (mm)'); xlabel('Return period (years)'); grid on;
    set(gca, 'XScale', 'log','XTick',[2 5 10 25 50 100 200],'FontSize',14,'FontName','Times New Roman');
    legend([h1 h2 h3],'TC Events','Non-TC Events','Combined','Location','northwest'); title("Return levels of Rainfall");
    clearvars mu RL 
    
    saveas(hh1,"Return_Levels_rainfall");


    % calculating combined return levels for given return periods
    RL_RF(1,:) = Q_RP;
    RL_RF(2,:) = interp1(Com_RP,RF,Q_RP,"linear");

    Data_TC_ETC_RF(4,:)=interp1(Com_RP,RF,RP,"linear");
end


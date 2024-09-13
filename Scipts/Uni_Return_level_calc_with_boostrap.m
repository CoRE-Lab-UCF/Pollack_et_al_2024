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

function [RL_NTR_com,RL_RF_com, RL_TC_NTR, RL_ETC_NTR, RL_TC_RF, RL_ETC_RF,RP_NTR,RP_RF]=Uni_Return_level_calc_with_boostrap_new(n_sim,ET_NTR,TC_NTR,ET_RF,TC_RF, Thres_NTR, Thres_RF,n_years,Q_RP,l_b_NTR,U_b_NTR,l_b_RF,U_b_RF)
    n=n_years;
    for k=1:n_sim
        
        
        % For TC events
        TC_NTR_Bt = datasample(TC_NTR(:,2),length(TC_NTR(:,2)));
    
        [x_TC_NTR,is]= sort(TC_NTR_Bt(:,1));
        %x_TC_time = TC_NTR(is,1);
        x_TC_NTR= x_TC_NTR -Thres_NTR;
        x_TC_NTR(x_TC_NTR<=0) = 0.001;
        [parm_TC_NTR,bounds_TC_NTR] = gpfit(x_TC_NTR,0.05);
        
        
        %For ETC events
        ETC_NTR_Bt = datasample(ET_NTR(:,2),length(ET_NTR(:,2)));
    
        [x_ETC_NTR,is]= sort(ETC_NTR_Bt(:,1));
        %x_ETC_time = ET_NTR(is,1);
        x_ETC_NTR= x_ETC_NTR -Thres_NTR;
        x_ETC_NTR(x_ETC_NTR<=0) = 0.001;
        [parm_ETC_NTR,bounds_ETC_NTR] = gpfit(x_ETC_NTR,0.05);
            
        
        % Return level Calculation
        
        % Annual Exceedence prbability vector
        RP= (0:0.1:max(Q_RP));
        A_NEP = 1./RP;
        
        % Calculating the Inter Arrival Time (years)
        mu_TC = n/length(x_TC_NTR);
        mu_ETC = n/length(x_ETC_NTR);
        
        % Calculating Non_exceedence probabiliy [ No_Ex_prob = 1-(IAT*AEP)]
        Prob_TC= 1-mu_TC.*A_NEP;
        Prob_ETC= 1-mu_ETC.*A_NEP;
    
        
        % Return levels
        RL_TC_NTR(k,:)= gpinv(Prob_TC,parm_TC_NTR(1),parm_TC_NTR(2),Thres_NTR);
        RL_ETC_NTR(k,:)= gpinv(Prob_ETC,parm_ETC_NTR(1),parm_ETC_NTR(2),Thres_NTR);
        
        
        % Combining Retuen periods
        NTR = l_b_NTR:0.1:U_b_NTR*10; % Descritized NTR vector
    
        ANEP_TC =1-((1-gpcdf(NTR,parm_TC_NTR(1),parm_TC_NTR(2),Thres_NTR))./mu_TC); % Calculating AEP
        ANEP_ETC =1-((1-gpcdf(NTR,parm_ETC_NTR(1),parm_ETC_NTR(2),Thres_NTR))./mu_ETC);
        
        Com_AEP = 1-(ANEP_ETC.*ANEP_TC); % Combining AEP
        
        Com_RP = 1./Com_AEP; % Converting AEP in to Return Period
         
        RL_NTR_com(1,:) = NTR;
        RL_NTR_com(k+1,:) = Com_RP;
        RP_NTR = RP;
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        %% Fitting GPDs to RF
        % For TC events
        TC_RF_Bt = datasample(TC_RF(:,4),length(TC_RF(:,4)));
        [x_TC_RF,is]= sort(TC_RF_Bt(:,1));
        %x_TC_time = TC_RF(is,3);
        x_TC_RF= x_TC_RF -Thres_RF;
        x_TC_RF(x_TC_RF<=0) = 0.001;
        [parm_TC_RF,bounds_TC_RF] = gpfit(x_TC_RF,0.05);
        
        
        %For ETC events
        ETC_RF_Bt = datasample(ET_RF(:,4),length(ET_RF(:,4)));
        [x_ETC_RF,is]= sort(ETC_RF_Bt(:,1));
        %x_ETC_time = ET_RF(is,3);
        x_ETC_RF= x_ETC_RF -Thres_RF;
        x_ETC_RF(x_ETC_RF<=0) = 0.001;
        [parm_ETC_RF,bounds_ETC_RF] = gpfit(x_ETC_RF,0.05);
        
    
        
        % Return level Calculation
        
        % Annual Exceedence prbability vector
        RP= (0:0.1:max(Q_RP));
        A_NEP = 1./RP;
        
        % Calculating the Inter Arrival Time (years)
        mu_TC = n/length(x_TC_RF);
        mu_ETC = n/length(x_ETC_RF);
        
        % Calculating Non_exceedence probabiliy [ No_Ex_prob = 1-(IAT*AEP)]
        Prob_TC= 1-mu_TC.*A_NEP;
        Prob_ETC= 1-mu_ETC.*A_NEP;
        
        % Return levels
        RL_TC_RF(k,:)= gpinv(Prob_TC,parm_TC_RF(1),parm_TC_RF(2),Thres_RF);
        RL_ETC_RF(k,:)= gpinv(Prob_ETC,parm_ETC_RF(1),parm_ETC_RF(2),Thres_RF);
        
        
        % Combining Retuen periods
        RF = l_b_RF:0.1:U_b_RF*10; % Descritized RF vector
        
        ANEP_TC = 1-((1-gpcdf(RF,parm_TC_RF(1),parm_TC_RF(2),Thres_RF))./mu_TC); % Calculating AEP
        ANEP_ETC = 1-((1-gpcdf(RF,parm_ETC_RF(1),parm_ETC_RF(2),Thres_RF))./mu_ETC);
        Com_AEP =1-(ANEP_ETC.*ANEP_TC); % Combining AEP
        
        Com_RP = 1./Com_AEP; % Converting AEP in to Return Period
        
    
        % calculating combined return levels for given return periods
        RL_RF_com(1,:) = RF;
        RL_RF_com(k+1,:) = Com_RP;
        RP_RF = RP;
    end
end


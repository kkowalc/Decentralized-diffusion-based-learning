clear all; close all;clc; set(0,'DefaultFigureWindowStyle','docked'); %#ok<CLALL> 

global M; M = 50;               %#ok<GVMIS> Number of nodes
global D; D= [0,10];            %#ok<GVMIS> % Total inference domain (D in a paper)
global Net_tplgy;               %#ok<GVMIS> Network topology
global L; L = 1;                %#ok<GVMIS> Lipschitz constant
plot_agent = 18;                %Index of agent for plotting
N_total = 1000;                % Max no. of local meas. (for every single agent)
Start_Shr = 1;                  % Start sharing when tb> Start_Shr
rng_seed = 2610.1774;           % rand*100000
m = @(x)sin(x).*exp(-0.2*x) + 3;% The latent nonlinear phenomenon
Net_tplgy = random_topology_networkV5(0.1, rng_seed,'force_connectivity'); % Generate network topology
new_explanatory_observationV6([], D, M,rng_seed, 'Initialize'); % Input generator initialization;    
T_hor = N_total;                % Time horizon
Sigma_e = 0.7*ones(size(1:M)).*rand(size(1:M)); % Noise disperssion
delta = 0.01;                   % All bounds with probability (1-delta)
New_meas_period = 1;            % New measurements only at every multipl. of 'New_meas_period'
Sharing_period = 1;             % Sharing is made only at every multipl. of 'Sharing_period'
ksi = D(1):0.01:D(2);           % A priori selected estimation points
B_kt = nan(M,length(ksi));      % Final confidence bounds for all agents  
B_loc = nan(M,length(ksi));     % Local confidence bound for selected agent
Loc_data = nan(M,T_hor,5);      % Local explanatory and output values - k,t,(ksi,y,y_bar,B_kt) 
Shr_tuple = cell(M,M);          % Tuples selected for sharing
Req_points = zeros(M,1);        % Domain points requested by the agents from the neighbors
Acq_data = cell(M,T_hor);       % Acquired tuples (xi,y_bar,B)

for t = 1:T_hor                                                  % LOCAL MEAS. AND SELECTION OF TUPLES FOR SHARING
   disp(['Time step: ',num2str(t)]);
   for k = 1:M                                                   % Select requested argument
      if mod(t,Sharing_period)==0
      Req_points(k) = select_arg2requestV6(squeeze(Loc_data(k,:,:)),Acq_data{k},'urand',0); % Select required arguments by all agents
      end
   end
   for k = 1:M                                                   % New local measurements (for all agents)
      n_local = sum(~isnan(Loc_data(k,:,1)));                    % Total number of local meas. of agent k
      if (mod(t,New_meas_period)==0)&&(n_local<N_total)          % Check whether new meas. should be taken
         n_local = n_local + 1;
         ksi_kt = new_explanatory_observationV6(k,D,M,rng_seed,'Sample'); % New local explanatory data   
         y = new_output(ksi_kt, Sigma_e(k),m);                   % New local outcome data
         Loc_data(k,n_local,1:2) = [ksi_kt,y];                   % New meas. in local collection
         if t>Start_Shr
            [y_bar, B, h] = NadarayaWatsonV7(ksi_kt,Loc_data(k,:,1),Loc_data(k,:,2),NaN,L,delta,Sigma_e(k)); %NaN - for autom. h tuning
            Loc_data(k,n_local,3:5) = [y_bar,B,h];                 % New est. & B in local collection
         end
      end
      if (mod(t,Sharing_period)==0)&&(t>Start_Shr)
         Shr_tuple = share_tuple_V4(Shr_tuple, k,Loc_data,Acq_data,'',Req_points); % Select tuple for sharing
      end
   end
   if mod(t,Sharing_period)==0
      Acq_data = Append_Acq_dataV5(Acq_data,Shr_tuple);       
   end

end


for k = 1:M                                                    % NONLINEARITY ESTIMATION BY k's 
   if k==plot_agent
   ksi_loc_max=find(ksi<=max(Loc_data(k,:,1)),1,'last');
   ksi_loc_min=find(ksi>=min(Loc_data(k,:,1)),1,'first');
   [mu_loc, B_loc, ~]=FinalEstimateV1(ksi, Loc_data, cell(M,T_hor), B_loc, k, L, delta, Sigma_e);
    [mu, B_kt, Acq_flag]=FinalEstimateV1(ksi, Loc_data, Acq_data, B_kt, k, L, delta, Sigma_e);
    color_ind_loc =(Acq_flag==0);color_ind_acq =(Acq_flag==1);

    figure(2); set(gcf, 'NumberTitle', 'off','Name', ['General, k=',num2str(k)],'color','w');
    subplot(3,1,[1 2]);   
    plot(Loc_data(k,1:min(n_local,500),1),Loc_data(k,1:min(n_local,500),2),'.','color',[0 0.4470 0.7410], 'MarkerSize',5,'DisplayName','Noisy meas.');hold on;grid on;  
    subplot(3,1,[1 2]);
    plot(ksi,mu,'-','color',[0.1 0.6 0.2], 'DisplayName','$\hat{m}_{k,t}(x)$','LineWidth', 1.5);
    subplot(3,1,[1 2]);
    plot(ksi,m(ksi),'k--','DisplayName','$m(x)$','LineWidth', 2)%title(['Agent: ',num2str(k), ', t = ', num2str(T_hor)]);
    subplot(3,1,[1 2]);
    plot (ksi,mu(:)'+B_kt(k,:),'color',[0.9 0.5 0.2 0.7],'DisplayName','Bounds','LineWidth', 1);%annotation('textbox',[0.73,0.8,0.13,0.05],'String',['max(B) = ',num2str(max(B_kt(k,:)))],'FitBoxToText','on');
    subplot(3,1,[1 2]);
    plot (ksi,mu(:)'-B_kt(k,:),'color',[0.9 0.5 0.2 0.7],'HandleVisibility','off','LineWidth', 1); 
    subplot(3,1,[1 2]);
    xline(ksi(ksi_loc_max),"k--",'DisplayName','Local domain');
    subplot(3,1,[1 2]);
    xline(ksi(ksi_loc_min),"k--",'HandleVisibility','off'); ylim([min(Loc_data(k,:,2)) max(Loc_data(k,:,2))]);
    ylim([2 4]);
    subplot(3,1,3);histogram(Loc_data(k,:,1),'EdgeColor',[0 0.4470 0.7410],'facecolor',[0 0.4470 0.7410],'Normalization','pdf','DisplayName','Hist. Local');xlim([min(ksi),max(ksi)]);hold on;
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');xlim([D(1),D(2)]);grid on;
    
    figure(3); set(gcf, 'NumberTitle', 'off','Name', ['LocalEst, k=',num2str(k)],'color','w'); 
    plot(Loc_data(k,1:min(n_local,500),1),Loc_data(k,1:min(n_local,500),2),'.','color',[0 0.4470 0.7410], 'MarkerSize',5,'DisplayName','Noisy meas.');hold on;grid on;  
    plot(ksi,mu,'-','color',[0.1 0.6 0.2], 'DisplayName','$\hat{m}_{k,t}(x)$','LineWidth', 1.5);
    plot(ksi,m(ksi),'k--','DisplayName','$m(x)$','LineWidth', 2)%title(['Agent: ',num2str(k), ', t = ', num2str(T_hor)]);
    plot (ksi,mu(:)'+B_kt(k,:),'color',[0.9 0.5 0.2 0.7],'DisplayName','Bounds','LineWidth', 1);%annotation('textbox',[0.73,0.8,0.13,0.05],'String',['max(B) = ',num2str(max(B_kt(k,:)))],'FitBoxToText','on');
    plot (ksi,mu(:)'-B_kt(k,:),'color',[0.9 0.5 0.2 0.7],'HandleVisibility','off','LineWidth', 1); 
    xline(ksi(ksi_loc_max),"k--",'DisplayName','Local domain');
    xline(ksi(ksi_loc_min),"k--",'HandleVisibility','off'); ylim([min(Loc_data(k,:,2)) max(Loc_data(k,:,2))]);
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');xlim([D(1),D(2)]);grid on;
    ylim([2 4]);
    plot(ksi,mu_loc,'color',[0.3010 0.7450 0.9330],'MarkerSize',10, 'DisplayName','Local estimate','LineWidth', 1.5);


    figure(4); set(gcf, 'NumberTitle', 'off','Name', ['LocalBounds, k=',num2str(k)],'color','w'); 
    plot(Loc_data(k,1:min(n_local,500),1),Loc_data(k,1:min(n_local,500),2),'.','color',[0 0.4470 0.7410], 'MarkerSize',5,'DisplayName','Noisy meas.');hold on;grid on;  
    plot(ksi,mu,'-','color',[0.1 0.6 0.2], 'DisplayName','$\hat{m}_{k,t}(x)$','LineWidth', 1.5);
    plot(ksi,m(ksi),'k--','DisplayName','$m(x)$','LineWidth', 2)%title(['Agent: ',num2str(k), ', t = ', num2str(T_hor)]);
    plot (ksi,mu(:)'+B_kt(k,:),'color',[0.9 0.5 0.2 0.7],'DisplayName','Bounds','LineWidth', 1);%annotation('textbox',[0.73,0.8,0.13,0.05],'String',['max(B) = ',num2str(max(B_kt(k,:)))],'FitBoxToText','on');
    plot (ksi,mu(:)'-B_kt(k,:),'color',[0.9 0.5 0.2 0.7],'HandleVisibility','off','LineWidth', 1); 
    xline(ksi(ksi_loc_max),"k--",'DisplayName','Local domain');
    xline(ksi(ksi_loc_min),"k--",'HandleVisibility','off'); ylim([min(Loc_data(k,:,2)) max(Loc_data(k,:,2))]);
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');xlim([D(1),D(2)]);grid on;
    ylim([1.4 3.8]); 
    plot (ksi(ksi_loc_min:ksi_loc_max),mu_loc(ksi_loc_min:ksi_loc_max)+B_loc(k,ksi_loc_min:ksi_loc_max),'-','color',[0.6350 0.0780 0.1840],'DisplayName','Local bounds','LineWidth', 1.5);%annotation('textbox',[0.73,0.8,0.13,0.05],'String',['max(B) = ',num2str(max(B_kt(k,:)))],'FitBoxToText','on');
    plot (ksi(ksi_loc_min:ksi_loc_max),mu_loc(ksi_loc_min:ksi_loc_max)-B_loc(k,ksi_loc_min:ksi_loc_max),'-','color',[0.6350 0.0780 0.1840],'HandleVisibility','off','LineWidth', 1.5); 

   
    if ~isempty(Acq_data{k})
       Acq_data_k_ksi = Acq_data{k}(:,1);
       Acq_data_k_y_est = Acq_data{k}(:,2);
       Acq_data_k_B = Acq_data{k}(:,3);

       figure(2);
       subplot(3,1,3);histogram(Acq_data_k_ksi,'Normalization','pdf','facecolor',[0.4940 0.1840 0.5560],'edgecolor',[0.4940 0.1840 0.5560],'facealpha',0.7,'DisplayName','Hist Acq.');legend;
       subplot(3,1,[1 2]);
       plot(Acq_data_k_ksi,Acq_data_k_y_est,'color',[0.4940 0.1840 0.5560,0.2],'MarkerIndices',1:5:length(Acq_data_k_y_est),'Marker','.','LineStyle','none','DisplayName','Acq. $(\xi^{acq}, \hat{\mu}^{acq})$');
       set(gca,'fontsize',14,'TickLabelInterpreter','latex');xlim([D(1),D(2)]);
       legend('Interpreter','latex');
        
       figure(3);
       plot(Acq_data_k_ksi,Acq_data_k_y_est,'color',[0.4940 0.1840 0.5560,0.2],'MarkerIndices',1:5:length(Acq_data_k_y_est),'Marker','.','LineStyle','none','DisplayName','Acq. $(\xi^{acq}, \hat{\mu}^{acq})$');
       set(gca,'fontsize',14,'TickLabelInterpreter','latex');xlim([D(1),D(2)]);
       legend('Interpreter','latex');

       figure(4);
       plot(Acq_data_k_ksi,Acq_data_k_y_est,'color',[0.4940 0.1840 0.5560,0.2],'MarkerIndices',1:5:length(Acq_data_k_y_est),'Marker','.','LineStyle','none','DisplayName','Acq. $(\xi^{acq}, \hat{\mu}^{acq})$');
       set(gca,'fontsize',14,'TickLabelInterpreter','latex');xlim([2.5 6]);
       legend('Interpreter','latex');

    end
    end
end 

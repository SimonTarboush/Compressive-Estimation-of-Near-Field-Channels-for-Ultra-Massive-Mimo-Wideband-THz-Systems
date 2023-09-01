%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Simon Tarboush                                      
% Last Modified: June, 2023
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the paper: Simon Tarboush, Anum Ali, Tareq Y. Al-Naffouri, 
% "Compressive Estimation of Near Field Channels for Ultra Massive-MIMO Wideband THz Systems", 
% ICASSP 2023 - 2023 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).
%
% If you use the channel simulator code "TeraMIMO" or any (modified) part of it in any publication, please cite 
% the paper: Simon Tarboush, Hadi Sarieddeen, Hui Chen, Mohamed Habib Loukil, Hakim Jemaa, Mohamed-Slim Alouini, Tareq Y. Al-Naffouri
% "TeraMIMO: A Channel Simulator for Wideband Ultra-Massive MIMO Terahertz Communications",
% IEEE Transactions on Vehicular Technology.
%
% Contact person email: simon.w.tarboush@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the results presented in Fig. 3 in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;clear;close all;
%% Add Paths
path(pathdef); addpath(pwd);
cd Algorithms; addpath(genpath(pwd)); cd ..; cd TeraMIMO_Channel; addpath(genpath(pwd)); cd ..; cd Molecular_Absorption;addpath(genpath(pwd)); cd ..; cd Utilities; addpath(genpath(pwd)); cd ..;
%% Initialize Channel Parameters
p_ch = generate_channel_param_TIV();
p_ch = update_channel_param_TIV(p_ch);
%% Calculation of Absorption Coefficient 
K_abs = get_Abs_Coef(p_ch);
%% Parameters related to the transmitter and receiver uniform linear/planar arrays
get_array_trans_param;
%% Main Simulation Parameters
E = 1e1; % Number of trials/runs/iterations
% The results in the paper were obtained with E=1e2
% Note that the simulation is too heavy and requires strong computing capabilities and long simulation time due to the used size of arrays
Lbar = 10; % Estimation of the sparsity level value (for OMP/SOMP based on for-loop)
Tx_power_tot_dBm = 4.4957; % Tx power (dBm)
DNR_dB = 0; % Data to Noise Ratio (dB)
G_hv_TR_OVS = 2; % Dictionary oversampling
% Here we define G_hv_TR_OVS = G^h_OVS = G^v_OVS for both Tx and Rx and then G_OVS = G^h_OVS * G^v_OVS
% For more flexibility, we can define G_h_T_OVS, G_v_T_OVS, G_h_R_OVS, G_h_R_OVS
num_multG = 2; % Grid oversampling: for example G_T_h=num_multG x Q_T_h : num_multG >= 1
% For more flexibility, we can define num_multG_h_T, num_multG_v_T, num_multG_h_R, num_multG_h_R
Q_T_quant = 2; % Number of Tx PSs quantization bits
Q_R_quant = 2; % Number of Rx PSs quantization bits
%% Measurements/Beams during the training/data transmission phase
Ns_train = 1; % Number of spatial streams during training, always equal to 1
Ns_data = min(Q_T,Q_R); % The number of spatial streams during data transmission Ns
N_RF_T = Ns_data; % will be equal to the number of Tx SAs during data transmission stage
N_RF_R = Q_R; % will be equal to the number of Rx SAs during data transmission stage
% This simulation is versus the number of measurements
M_T_Beam = [4:4:40 48 56 Qbar_T]; % Number of exhaustive-search training Tx (Full Training)/(Partial Training)
M_R_Beam = [4:4:40 48 56 Qbar_R]; % Number of exhaustive-search training Rx (Full Training)/(Partial Training)
Compression_RatioTx = 1; Compression_RatioRx = 1;
Sel_SAs_IndTx = 1:Q_T; Sel_SAs_IndRx = 1:Q_R;
Num_of_UsedSATx = length(Sel_SAs_IndTx); Num_of_UsedSARx = length(Sel_SAs_IndRx);
MeasTx_vec = length(M_T_Beam);
MeasRx_vec = length(M_R_Beam);
M_T_measMAT = zeros(Q_T,MeasTx_vec); % The number of measurements at Tx
M_R_measMAT = zeros(Q_R,MeasRx_vec); % The number of measurements at Rx
M_T_measMAT(Sel_SAs_IndTx,:) =  repmat(Compression_RatioTx*M_T_Beam,Num_of_UsedSATx,1); % The number of measurements at Tx
M_R_measMAT(Sel_SAs_IndRx,:) =  repmat(Compression_RatioRx*M_R_Beam,Num_of_UsedSARx,1); % The number of measurements at Rx
%% Parameters related to the GRID of OMP/SOMP algorithm and Quantization of Phase for the Codebooks/Dictionaries
get_grid_quant_param;
%% Construct the Tx and Rx Dictionaries
get_TxRxDict;
%% Define the SNR
get_SNR_THzChannel;
%% Initialize Arrays for Results: NMSE and EAR
% AR (this will be transformed later to EAR)
AR_Perf_CSI = zeros(MeasTx_vec,E,K); % perfect channel state information
AR_OMP = zeros(MeasTx_vec,E,K); 
AR_OMP_DR = zeros(MeasTx_vec,E,K); 
AR_SOMP = zeros(MeasTx_vec,E,K); 
AR_SOMP_DR = zeros(MeasTx_vec,E,K);
%% Main Loop
tic;
noise_pwr = N0_sc;
Tx_pwr = Tx_power_tot/Ns_train;
% Tx_pwr is "P_T" defined in Eq. 7 in the manuscript, i.e., the average Tx power used per transmission during the training phase
% during the training there is only one SA to be used, i.e., number of streams during training is Ns_train = 1
Data_pwr = 10.^(DNR_dB/10)*noise_pwr/Pr_avg_scnoprenocomb;
for indx_meas = 1:MeasTx_vec
    Tx_meas_indx = indx_meas;
    Rx_meas_indx = indx_meas;
    for indx_iter = 1:E
        % Display the current iteration
        [indx_meas indx_iter]
        % AoSAs THz channel generation
        [CH_Response, AT_PCSI, AR_PCSI] = channel_TIV(p_ch, K_abs);
        H_AoSA = cell2mat(CH_Response.H);        
        % This is the overall complex UM-MIMO channel matrix H defined in Eq. 1 in the manuscript
        % Proposed and Conventional Estimation Estimation Methods
        get_proposed_HSPM_estimation;
        % Performance Evaluation Metrics
        % AR Computations for the UM-MIMO        
        AR_Perf_CSI(indx_meas,indx_iter,:) = get_AchievableRate(p_ch, K, Tx_pwr, Data_pwr, noise_pwr, H_SP_UM, H_SP_UM, AT_PCSI, AR_PCSI); 
        AR_OMP(indx_meas,indx_iter,:) = get_AchievableRate(p_ch, K, Tx_pwr, Data_pwr, noise_pwr, H_Est_OMP_UM, H_SP_UM, F_Sel_OMP, W_Sel_OMP);
        AR_OMP_DR(indx_meas,indx_iter,:) = get_AchievableRate(p_ch, K, Tx_pwr, Data_pwr, noise_pwr, H_Est_OMP_DR_UM, H_SP_UM, F_Sel_OMP_DR, W_Sel_OMP_DR);
        F_bmf_SOMP_EstCh = repmat(F_Sel_SOMP,[1 1 K]);
        W_cmb_SOMP_EstCh = repmat(W_Sel_SOMP,[1 1 K]);
        AR_SOMP(indx_meas,indx_iter,:) = get_AchievableRate(p_ch, K, Tx_pwr, Data_pwr, noise_pwr, H_Est_SOMP_UM, H_SP_UM, F_bmf_SOMP_EstCh, W_cmb_SOMP_EstCh);
        F_bmf_SOMP_DR_EstCh = repmat(F_Sel_SOMP_DR,[1 1 K]);
        W_cmb_SOMP_DR_EstCh = repmat(W_Sel_SOMP_DR,[1 1 K]);
        AR_SOMP_DR(indx_meas,indx_iter,:) = get_AchievableRate(p_ch, K, Tx_pwr, Data_pwr, noise_pwr, H_Est_SOMP_DR_UM, H_SP_UM, F_bmf_SOMP_DR_EstCh, W_cmb_SOMP_DR_EstCh);
    end
end
% Keep track of simulation time
Sim_Duration = toc;
hr = floor(Sim_Duration/3600);mint = floor((Sim_Duration - hr*3600)/60);sec = Sim_Duration - hr*3600 - mint*60;
fprintf('The simulation time is: %d hr %d min %f sec\n',hr,mint,sec);
%% Plot Results 
% EAR (R_eff in the manuscript, Eq. 12) versus SNR
R_eff_PCI = mean(AR_Perf_CSI,[2 3]);
T_coh = [512 1024 inf]; % the THz channel coherence time in symbols
meanAR_OMP = mean(AR_OMP,[2 3]);
meanAR_SOMP = mean(AR_SOMP,[2 3]);
meanAR_OMP_DR = mean(AR_OMP_DR,[2 3]);
meanAR_SOMP_DR = mean(AR_SOMP_DR,[2 3]);
Eta = max(0,1-(1./T_coh')*(sum(M_T_measMAT))); % the loss in achievable rate due to the training
% Eq. 12 in the manuscript: EAR metric
R_eff_OMP = Eta.*(ones(length(T_coh),1)*meanAR_OMP.');
R_eff_SOMP = Eta.*(ones(length(T_coh),1)*meanAR_SOMP.');
R_eff_OMP_DR = Eta.*(ones(length(T_coh),1)*meanAR_OMP_DR.');
R_eff_SOMP_DR = Eta.*(ones(length(T_coh),1)*meanAR_SOMP_DR.');
figure('color',[1,1,1]);
h_plt_PCI = plot(sum(M_T_measMAT),R_eff_PCI,'g+-','linewidth',2,'MarkerSize',8);hold on;
h_plt_OMP = plot(sum(M_T_measMAT),R_eff_OMP.','ko-','linewidth',2,'MarkerSize',8); set(h_plt_OMP,{'LineStyle'},{'-';'--';':'});
h_plt_SOMP = plot(sum(M_T_measMAT),R_eff_SOMP.','bx-','linewidth',2,'MarkerSize',8); set(h_plt_SOMP,{'LineStyle'},{'-';'--';':'});
h_plt_OMP_DR = plot(sum(M_T_measMAT),R_eff_OMP_DR.','ms-','linewidth',2,'MarkerSize',8); set(h_plt_OMP_DR,{'LineStyle'},{'-';'--';':'});
h_plt_SOMP_DR = plot(sum(M_T_measMAT),R_eff_SOMP_DR.','rd-','linewidth',2,'MarkerSize',8); set(h_plt_SOMP_DR,{'LineStyle'},{'-';'--';':'});
xlabel('Number of Measurements ($M_{\mathrm{T}}^{\textrm{tr}}$)','FontSize',14,'Interpreter','latex');
ylabel('$R_{\textrm{eff}}$ (b/s/Hz)','FontSize',14,'Interpreter','latex');
hplt = zeros(5, 1);
hplt(1) = plot(NaN,NaN,'g+-','linewidth',2,'MarkerSize',8);
hplt(2) = plot(NaN,NaN,'ko-','linewidth',2,'MarkerSize',8);
hplt(3) = plot(NaN,NaN,'bx-','linewidth',2,'MarkerSize',8);
hplt(4) = plot(NaN,NaN,'ms-','linewidth',2,'MarkerSize',8);
hplt(5) = plot(NaN,NaN,'rd-','linewidth',2,'MarkerSize',8);
le=legend(hplt,'Perfect CSI','OMP','SOMP','OMP-DR','SOMP-DR');
le.Interpreter='latex';le.Location='southwest';le.FontSize=14;
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=14;
box off;
grid on;
xlim([min(sum(M_T_measMAT)) max(sum(M_T_measMAT))]);
axis tight;
axlims=axis;
x_range=axlims(2)-axlims(1);y_range=axlims(4)-axlims(3);loseness=5;
axis([axlims(1)-loseness/1e2*x_range axlims(2)+loseness/1e2*x_range axlims(3)-loseness/1e2*y_range axlims(4)+loseness/1e2*y_range])
dim = [0.53 0.145 0 0];
str = {'$\textbf{solid}:~T_{\textrm{coh}}=\infty$',...
    '$\textbf{dashed}:~T_{\textrm{coh}}$ = 1024 symbols',...
    '$\textbf{dotted}:~T_{\textrm{coh}}$ = 512 symbols'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',12,'VerticalAlignment','bottom','HorizontalAlignment','left','Margin',3);
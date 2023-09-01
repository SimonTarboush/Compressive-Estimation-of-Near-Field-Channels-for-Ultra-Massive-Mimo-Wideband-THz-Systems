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
% This script generates the results presented in Fig. 2 in the manuscript
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
Tx_power_tot_dBm = -9.4957:2:10.4957; % Tx power (dBm)
DNR_dB = 0; % Data to Noise Ratio (dB)
G_hv_TR_OVS = 2; % Dictionary oversampling
% Here we define G_hv_TR_OVS = G^h_OVS = G^v_OVS for both Tx and Rx and then G_OVS = G^h_OVS ×G^v_OVS
% For more flexibility, we can define G_h_T_OVS, G_v_T_OVS, G_h_R_OVS, G_h_R_OVS
num_multG = 2; % Grid oversampling: for example G_T_h=num_multG x Q_T_h : num_multG >= 1
% For more flexibility, we can define num_multG_h_T, num_multG_v_T, num_multG_h_R, num_multG_h_R
Q_T_quant = 2; % Number of Tx PSs quantization bits
Q_R_quant = 2; % Number of Rx PSs quantization bits
%% Measurements/Beams during the training/data transmission phase
Ns_train = 1; % Number of spatial streams during training, always equal to 1
Ns_data = min(Q_T,Q_R); % The number of spatial streams during data transmission Ns
N_RF_T = Ns_data; % will be equal to the number of Tx SAs during data transmission stage
N_RF_R = Ns_data; % will be equal to the number of Rx SAs during data transmission stage
M_T_Beam = Qbar_T; % Number of exhaustive-search training Tx (Full Training)/(Partial Training)
M_R_Beam = Qbar_R; % Number of exhaustive-search training Rx (Full Training)/(Partial Training)
Compression_RatioTx = 1; Compression_RatioRx = 1;
Sel_SAs_IndTx = 1:Q_T; Sel_SAs_IndRx = 1:Q_R;
Num_of_UsedSATx = length(Sel_SAs_IndTx); Num_of_UsedSARx = length(Sel_SAs_IndRx);
M_T_measMAT = zeros(Q_T,1); % The number of measurements at Tx
M_R_measMAT = zeros(Q_R,1); % The number of measurements at Rx
M_T_measMAT(Sel_SAs_IndTx,1) = Compression_RatioTx*M_T_Beam; % The number of measurements at Tx
M_R_measMAT(Sel_SAs_IndRx,1) = Compression_RatioRx*M_R_Beam; % The number of measurements at Rx
%% Parameters related to the GRID of OMP/SOMP algorithm and Quantization of Phase for the Codebooks/Dictionaries
get_grid_quant_param;
%% Construct the Tx and Rx Dictionaries
get_TxRxDict;
%% Define the SNR
get_SNR_THzChannel;
%% Initialize Arrays for Results: NMSE and EAR
% NMSE
NMSE_OMP = zeros(SNR_len,E,K); 
NMSE_OMP_DR = zeros(SNR_len,E,K); 
NMSE_SOMP = zeros(SNR_len,E,K); 
NMSE_SOMP_DR = zeros(SNR_len,E,K);
%% Main Loop
tic;
for indx_snr = 1:SNR_len
    noise_pwr = N0_sc; % sigma^2
    Tx_pwr = Tx_power_tot(indx_snr)/Ns_train;
    % Tx_pwr is "P_T" defined in Eq. 7 in the manuscript, i.e., the average Tx power used per transmission during the training phase
    % during the training there is only one SA to be used, i.e., number of streams during training is Ns_train = 1
    for indx_iter = 1:E
        % Display the current iteration
        [indx_snr indx_iter]
        % AoSAs THz channel generation
        CH_Response = channel_TIV(p_ch, K_abs);
        H_AoSA = cell2mat(CH_Response.H);        
        % This is the overall complex UM-MIMO channel matrix H defined in Eq. 1 in the manuscript
        % Proposed and Conventional Estimation Estimation Methods
        get_proposed_HSPM_estimation;
        % Performance Evaluation Metrics
        % NMSE Computations for the UM-MIMO
        for indx_subc = 1:K
            % Eq. 11 in the manuscript
            NMSE_OMP(indx_snr,indx_iter,indx_subc) = norm(cell2mat(H_Est_OMP_UM(:,:,indx_subc))-cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2;
            NMSE_OMP_DR(indx_snr,indx_iter,indx_subc) = norm(cell2mat(H_Est_OMP_DR_UM(:,:,indx_subc))-cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2;
            NMSE_SOMP(indx_snr,indx_iter,indx_subc) = norm(cell2mat(H_Est_SOMP_UM(:,:,indx_subc))-cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2;
            NMSE_SOMP_DR(indx_snr,indx_iter,indx_subc) = norm(cell2mat(H_Est_SOMP_DR_UM(:,:,indx_subc))-cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_SP_UM(:,:,indx_subc)),'fro')^2;
        end
    end
end
% Keep track of simulation time
Sim_Duration = toc;
hr = floor(Sim_Duration/3600);mint = floor((Sim_Duration - hr*3600)/60);sec = Sim_Duration - hr*3600 - mint*60;
fprintf('The simulation time is: %d hr %d min %f sec\n',hr,mint,sec);
%% Plot Results 
% NMSE versus SNR
figure('color',[1,1,1]);
plot(SNR_dB,10*log10(mean(NMSE_OMP,[2 3])),'ko-','linewidth',2,'MarkerSize',10);hold on;
plot(SNR_dB,10*log10(mean(NMSE_SOMP,[2 3])),'bx-','linewidth',2,'MarkerSize',10);
plot(SNR_dB,10*log10(mean(NMSE_OMP_DR,[2 3])),'ms-','linewidth',2,'MarkerSize',10);
plot(SNR_dB,10*log10(mean(NMSE_SOMP_DR,[2 3])),'rd-','linewidth',2,'MarkerSize',10);
xlabel('SNR (dB)','FontSize',14,'Interpreter','latex');
ylabel('NMSE (dB)','FontSize',14,'Interpreter','latex');
le=legend('OMP','SOMP','OMP-DR','SOMP-DR');le.Interpreter='latex';le.Location='northeast';le.FontSize=14;
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=14;
box off;
grid on;
xlim([min(SNR_dB) max(SNR_dB)]);
axis tight;
axlims=axis;
x_range=axlims(2)-axlims(1);y_range=axlims(4)-axlims(3);loseness=5;
axis([axlims(1)-loseness/1e2*x_range axlims(2)+loseness/1e2*x_range axlims(3)-loseness/1e2*y_range axlims(4)+loseness/1e2*y_range]);
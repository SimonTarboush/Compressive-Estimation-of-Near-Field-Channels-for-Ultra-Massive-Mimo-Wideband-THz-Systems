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
% This script generates information related to noise power, path loss, and SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subcarriers and bandwidth
K = p_ch.Nsub_c;     % Number of subcarriers
B_sys = p_ch.BW;     % System Bandwidth
% Noise Power
N0_dBm = -173.8+10*log10(B_sys); % Noise power (dBm) % 10*log10(T*Kb) + 30 = -173.8
N0_sc = 10.^((N0_dBm-30)/10); % in Watts
% Tx Power Scalar
Tx_power_tot = 10.^((Tx_power_tot_dBm-30)/10); % Total Tx power in Watt
% PL
lambda  = p_ch.c/p_ch.Fc; % wavelength
PLE_val = p_ch.PLE; % path loss exponent
d_TxRx3D = p_ch.d_tx_rx; % 3D distance between transmitter and receiver
PL = ((4*pi*d_TxRx3D)/lambda)^PLE_val*exp(K_abs(ceil(K/2),1)*d_TxRx3D); % Path loss check TeraMIMO IEEE_TVT for the definition
PL_dB = 10*log10(PL);
% Received Power
Pr_avg_dBm = Tx_power_tot_dBm + 10*log10(Qbar_T) + 10*log10(Qbar_R) - PL_dB; % Average Received Power
Pr_avg_dBmnoprenocomb = Tx_power_tot_dBm - PL_dB; 
Pr_avg_sc = 10.^((Pr_avg_dBm-30)/10); % Average Received power in Watts
Pr_avg_scnoprenocomb = 10.^((Pr_avg_dBmnoprenocomb-30)/10);
% Rx SNR
Rx_SNR_dB = Pr_avg_dBm - N0_dBm; % Rx SNR in dB
Rx_SNRnoprenocomb_dB = Pr_avg_dBmnoprenocomb - N0_dBm; % Rx SNR no combining and no precoding in dB
Rx_SNR_sc = 10.^(Rx_SNR_dB/10); % Rx SNR in linear scale
Rx_SNRnoprenocomb_sc = 10.^(Rx_SNRnoprenocomb_dB/10); % Rx SNR no combining and no precoding in linear scale
SNR_dB = Rx_SNRnoprenocomb_dB-10*log10(K); % SNR to plot the results
SNR_len = length(Rx_SNR_sc);
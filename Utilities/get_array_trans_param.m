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
% This script generates simulation paramters related to the UM array dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Antennas Parameters
d_ae = 0.5; % antenna elements spacing, for DFT matrix d_ae normalized to lambda/2
% Total Number of SAs
Q_T = p_ch.Tx_AoSA.Q;  % Equal to: Q_T^h*Q_T^v; % Number of total SA Tx antennas
Q_R = p_ch.Rx_AoSA.Q;  % Equal to: Q_R^h*Q_R^v; % Number of total SA Rx antennas
Qbar_T = p_ch.Tx_AoSA.Qbar;% Equal to: Qbar_T^h*Qbar_T^v; % Number of total Tx antenna elements (AEs) per SA
Qbar_R = p_ch.Rx_AoSA.Qbar;% Equal to: Qbar_R^h*Qbar_R^v; % Number of total Rx antenna elements (AEs) per SA
% Number of SAs
Q_T_v = p_ch.Tx_AoSA.Qdim(1); % Number of transmit SAs on the z-axis of planar array
Q_T_h = p_ch.Tx_AoSA.Qdim(2); % Number of transmit SAs on the y-axis of planar array
Q_R_v = p_ch.Rx_AoSA.Qdim(1); % Number of receiver SAs on the z-axis of planar array
Q_R_h = p_ch.Rx_AoSA.Qdim(2); % Number of receiver SAs on the y-axis of planar array
% Number of AEs
Qbar_T_v = p_ch.Tx_AoSA.Qbardim(1); % Number of transmit antenna elements (AEs) per SA on the z-axis of planar array
Qbar_T_h = p_ch.Tx_AoSA.Qbardim(2); % Number of transmit antenna elements (AEs) per SA on the y-axis of planar array
Qbar_R_v = p_ch.Rx_AoSA.Qbardim(1); % Number of receiver antenna elements (AEs) per SA on the z-axis of planar array
Qbar_R_h = p_ch.Rx_AoSA.Qbardim(2); % Number of receiver antenna elements (AEs) per SA on the y-axis of planar array
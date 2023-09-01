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
% This script generates simulation paramters related to on-grid OMP/SOMP algorithm and phase quantization of Tx/Rx dictionaries used for compressed sensing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID parameters for OMP/SOMP
% Note that for efficient on-grid algorithm: G >= max(Qbar_T,Qbar_R)
% number of grids G_T and G_R
G_T_v = num_multG*Qbar_T_v; % Number of Quantization levels for Tx vertical dimension dictionary of UPA per SA
G_T_h = num_multG*Qbar_T_h; % Number of Quantization levels for Tx horizontal dimension dictionary of UPA per SA
G_R_v = num_multG*Qbar_R_v; % Number of Quantization levels for Rx vertical dimension dictionary of UPA per SA
G_R_h = num_multG*Qbar_R_h; % Number of Quantization levels for Rx horizontal dimension dictionary of UPA per SA
G_T = G_T_h*G_T_v;          % Total Number of Quantization levels for Tx dictionary
G_R = G_R_h*G_R_v;         % Total Number of Quantization levels for Rx dictionary
%% Quantization of Phase for the Dictionaries and RF Beamformers/Combiners
% The RF-beamformers/combiners and the Tx/Rx dictionaries (for OMP/SOMP) are implemented using finite resolution phase-shifters
QPhaseShifts_tx_sector_start = 0;QPhaseShifts_tx_sector_end = 2*pi;
QPhaseShifts_rx_sector_start = 0;QPhaseShifts_rx_sector_end = 2*pi;
% We add hardware constraints on the Analog Phase Shifters
QPhaseShifts_tx_sector_rng = QPhaseShifts_tx_sector_end-QPhaseShifts_tx_sector_start;
QPhaseShifts_rx_sector_rng = QPhaseShifts_rx_sector_end-QPhaseShifts_rx_sector_start;
% Possible values of Tx and Rx Quantized phase-shifts following finite resolution phase-shifters (PSs)
Quant_PS_T = QPhaseShifts_tx_sector_start:QPhaseShifts_tx_sector_rng/2^Q_T_quant:QPhaseShifts_tx_sector_end-QPhaseShifts_tx_sector_rng/2^Q_T_quant;
Quant_PS_R = QPhaseShifts_rx_sector_start:QPhaseShifts_rx_sector_rng/2^Q_R_quant:QPhaseShifts_rx_sector_end-QPhaseShifts_rx_sector_rng/2^Q_R_quant;
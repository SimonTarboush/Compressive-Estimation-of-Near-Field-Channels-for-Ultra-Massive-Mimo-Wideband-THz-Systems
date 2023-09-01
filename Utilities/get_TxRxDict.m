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
% This script generates the quantized Tx/Rx dictionaries \bar{A}_T/\bar{A}_R and the quantized dictionary matrix 
% \bar{\Theta} (Eq. 9 in the manuscript) used for compressed sensing estimation using OMP/SOMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the Tx/Rx Far-Field Dictionaries
Abar_T = get_DictionaryUPA(Quant_PS_T,Qbar_T_h,Qbar_T_v,G_T_h,G_T_v); % Check  Eq. 9 and Eq. 10 in the manuscript
Abar_R = get_DictionaryUPA(Quant_PS_R,Qbar_R_h,Qbar_R_v,G_R_h,G_R_v);% Check  Eq. 9 and Eq. 10 in the manuscript
% Get the oversampled dictionaries used in the proposed dictionary reduction (DR) method 
% Check the proposed method after Eq. 10 in the manuscript
Abar_T_OVS = get_DictionaryUPA(Quant_PS_T,Qbar_T_h,Qbar_T_v,G_hv_TR_OVS*G_T_h,G_hv_TR_OVS*G_T_v);
Abar_R_OVS = get_DictionaryUPA(Quant_PS_R,Qbar_R_h,Qbar_R_v,G_hv_TR_OVS*G_R_h,G_hv_TR_OVS*G_R_v);
%% Define the Quantized Dictionary Matrix
Thetabar = kron(Abar_T,Abar_R); % ATbar_kron_ARbar is the quantized dictionary matrix \bar{\Theta} (Eq. 9 in the manuscript)
% TeraMIMO channel implementation uses transpose for the Tx array response vectors (ARVs)
% if the channel model uses Hermitian instead of transpose, we should use conj(Abar_T)
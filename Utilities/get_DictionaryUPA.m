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
% This function generates the quantized Tx/Rx dictionaries \bar{A}_T/\bar{A}_R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Quantized_Phase_Shifts: The actual quantized phase shift values
% Qbar_h:    The number of antennas over horizontal dimension (Y-axis) of the UPA
% Qbar_v:    The number of antennas over vertical dimension (Z-axis) of the UPA
% G_h:          Grid size over horizontal dimension to quantize the angular region
% G_v:          Grid size over vertical dimension to quantize the angular region
% Output Arguments:
% A_bar: The dictionary \bar{A}_T/\bar{A}_R used in compressed sensing (Eq. 9 in the manuscript) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A_bar = get_DictionaryUPA(Quantized_Phase_Shifts, Qbar_h, Qbar_v, G_h, G_v) 
% The definitions are obtained from Eq. 4 in the manuscript but with
% quantization to a grid of size G_v/G_h instead of original size (AEs in the SA)

% Spatial Angles
horz_vec = -1/2+1/(2*G_h):1/G_h:1/2-1/(2*G_h);         % in ]-1/2, 1/2[
vert_vec = -1/2+1/(2*G_v):1/G_v:1/2-1/(2*G_v);           % in ]-1/2, 1/2[
A_ULA_resp = @(x,NumAnt) 1/sqrt(NumAnt)*exp(1j*2*pi*(0:NumAnt-1)'*x);
A_h = A_ULA_resp(horz_vec,Qbar_h);  % A_h discretize the azimuth plane
A_v = A_ULA_resp(vert_vec,Qbar_v);    % A_v discretize the elevation plane
A_UPA = kron(A_h,A_v);

% Other methods are to deloy the DFT dictionary as
% 1)
% A_UPA = dftmtx(Qbar_v*Qbar_h)/sqrt(Qbar_v*Qbar_h);
% 2)
% A_h = dftmtx(Qbar_h)/sqrt(Qbar_h); 
% A_v = dftmtx(Qbar_v)/sqrt(Qbar_v);
% A_UPA = kron(A_h,A_v);

% Obtain the closest phase in the quantized phase set due to finite resolution phase-shifters following Eq. 10 in the manuscript
Ang_A_tmp = wrapTo2Pi(angle(A_UPA));
Modified_Quan_PhaseShifts = [Quantized_Phase_Shifts Quantized_Phase_Shifts(end)+(Quantized_Phase_Shifts(end)-Quantized_Phase_Shifts(end-1))];
% Quantizing based on phase shifters' resolution
QuantAng_A_tmp = interp1(Modified_Quan_PhaseShifts, Modified_Quan_PhaseShifts, Ang_A_tmp, 'nearest');
A_bar = 1/sqrt(Qbar_h*Qbar_v)*exp(1j*QuantAng_A_tmp);
end
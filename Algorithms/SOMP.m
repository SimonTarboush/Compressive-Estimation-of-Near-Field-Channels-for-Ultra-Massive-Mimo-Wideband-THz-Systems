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
% This function implements the SOMP algorithm (solve the compressed sensing problem for all subcarriers simultaneously)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% y: Measurements vector
% Epsilon: The sensing matrix, equal to the multiplication of the measurement matrix by the quantized dictionary matrix (Psi and \bar{Theta} in the manuscript, Eq. 9) 
% Pilot_pwr: The transmission power in the training phase
% Lbar: Estimation of the sparsity level
% K: Number of subcarriers
% Output Arguments:
% H_Est: The estimation of the channel in beamspace
% Support: The index selected by the SOMP algorithm (known as the support)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_Est, Support] = SOMP(y,Epsilon,Pilot_pwr,Lbar,K)
R = y;
G = size(Epsilon,2);
Index = [];
% norm_Epsilon = sqrt(sum(abs(Epsilon).^2, 1)).';
for ell =1: Lbar % threshold = sparsity level
    % Find AoD/AoA pair
    [~, indx_c] = max(sum(abs(Epsilon'*R),2));  % Epsilon'*R./norm_Epsilon
    % Update AoD/AoA pair index
    Index = [Index indx_c]; 
    H_Est = zeros(G,K);
    % Estimate channel gains by solving the Least Square problem
    H_Est(Index,:) = 1/sqrt(Pilot_pwr)*pinv(Epsilon(:,Index)'*Epsilon(:,Index))*Epsilon(:,Index)'*y;
    % Update residual
    R = y - sqrt(Pilot_pwr)*Epsilon(:,Index)*H_Est(Index,:); 
end
Support = Index;
end
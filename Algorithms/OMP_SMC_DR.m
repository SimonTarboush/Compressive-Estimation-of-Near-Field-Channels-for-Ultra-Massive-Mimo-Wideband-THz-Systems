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
% This function implements the OMP algorithm compatible with the proposed dictionary reduction (DR) method over separate multi-carriers (solve the compressed sensing problem for evey subcarrier)
% where we get different sensing matrix for evey subcarrier
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
% Support: The index selected by the OMP algorithm with DR (known as the support)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_Est, Support] = OMP_SMC_DR(y,Epsilon,Pilot_pwr,Lbar,K)
R = y;
G = size(Epsilon,2);
H_Est = zeros(G,K);
Support = zeros(Lbar,K);
for indx_subc = 1:K
    % Solve the problem for every subcarrier
    Index = [];
    Epsilon_SC = Epsilon(:,:,indx_subc) ;
    for ell =1: Lbar % threshold = sparsity level
        % Find AoD/AoA pair
        [~, indx_c] = max(abs(Epsilon_SC'*R(:,indx_subc)));
        % Update AoD/AoA pair index
        Index = [Index indx_c]; 
        h = zeros(G,1);
        % Estimate channel gains by solving the Least Square problem
        h(Index,:) = 1/sqrt(Pilot_pwr)*pinv(Epsilon_SC(:,Index)'*Epsilon_SC(:,Index))*Epsilon_SC(:,Index)'*y(:,indx_subc);
        % Update residual
        R(:,indx_subc) = y(:,indx_subc) - sqrt(Pilot_pwr)*Epsilon_SC(:,Index)*h(Index,:);
    end
    Support(:,indx_subc) = Index;
    % Form the channel vector
    H_Est(:,indx_subc) = h; 
end
end

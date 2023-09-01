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
% This function computes random training codebook following the finite resolution phase shitfs constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% H: The frequency domain of MIMO THz channel response between Tx-Rx SAs
% N: Noise array
% P: Beamforming dictionary
% Z: Combining dictionary
% K: Number of subcarriers
% Tx_pow_train: The transmission power in the training phase
% Output Arguments:
% ZHP: The dictionary operated channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZHP = get_channel_output(H,N,P,Z,K,Tx_pow_train)
Y = zeros(size(Z,2),size(P,2),K);
for indx_subc = 1:K
    % Rx signal: rx = tx + Z^herm * noise is Eq. 7 in the manuscript
    Y(:,:,indx_subc) = sqrt(Tx_pow_train)*Z'*H(:,:,indx_subc)*P + Z'*N(:,:,indx_subc);
end
ZHP = reshape(Y(:),size(Z,2)*size(P,2),K); % apply vec() operator to get "y" as per Eq. 8 in the manuscript
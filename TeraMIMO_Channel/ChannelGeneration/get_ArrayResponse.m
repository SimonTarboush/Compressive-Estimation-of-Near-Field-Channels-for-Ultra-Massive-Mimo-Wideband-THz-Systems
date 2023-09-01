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
% This function computes the compute array response vector (AE level) based on the choosed antenna structure
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% indx_subc: Index of the subcarrier
% M: Number of receiver/transmitter SAs (rows) based on TRx
% N: Number of receiver/transmitter SAs (columns) based on TRx
% AngleIn: AoD/AoA for LoS/NLoS; in the case of NLoS, only for one ray in a cluster
% Psi_type: Define whether to calculate a beamsteering or beamforming vector 'SV', 'BF', as follows:
%                 1) @ Tx SV ---> AoD, @ Rx SV ---> AoA
%                 2) @ Tx BF ---> AoD_BF, @ Rx BF ---> AoA_BF
%                 In this version we only use SV
% TRx: Defines the direction of the link, 'T' at Tx side, 'R' at Rx side
% Output Arguments:
% a: Array response, a 3D-Array of size(M, N, num_freqs_per_subcarr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = get_ArrayResponse(p, indx_subc, M, N, AngleIn, Psi_type, TRx)
num_freqs_per_subcarr = p.nFreq(2);
% Initialize Output
a = zeros(M,N,num_freqs_per_subcarr);
if strcmp(TRx, 'T')
    delta1 = p.Tx_AoSA.delta(2);
    delta2 = p.Tx_AoSA.delta(1); 
elseif strcmp(TRx, 'R')
    delta1 = p.Rx_AoSA.delta(2);
    delta2 = p.Rx_AoSA.delta(1); 
else
    error('TRx has only two options: T/R');
end
if strcmp(Psi_type,'SV')
    lambda = p.lambda_subb(indx_subc,:);
end
% limit the angles to the physical one
% any angle out of the range [-90 90] are eliminated in azimuth plane
if AngleIn(1) > pi/2 || AngleIn(1) < -pi/2
    warning('Phi is out of range');
    a = zeros(M,N,num_freqs_per_subcarr);
else
    for m = 1:M
        for n = 1:N
            UnitDirVec = get_unitdirvec_from_anglevec(AngleIn);
            positionLocalAE = [0; (n-1-(N-1)/2)*delta1; (m-1-(M-1)/2)*delta2];
            Psi = (positionLocalAE.')*UnitDirVec;
            a(m,n,:) = exp(+1j*2*pi*Psi./lambda);
        end
    end
end
% Normalization  
SV_norm = sqrt(M*N);
a = a/SV_norm;
end

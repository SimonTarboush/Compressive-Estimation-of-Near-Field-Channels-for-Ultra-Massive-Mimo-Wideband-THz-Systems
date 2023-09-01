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
% This function generates the AoSA parameters
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct contains simulation parameters
% M: Number of subarrays (SAs) (row)
% N: Number of subarrays (SAs) (column)
% Mae: Number of antenna elements (AEs) (row) inside each SA 
% Nae: Number of antenna elements (AEs) (column) inside each SA
% DeltaM: SAs Spacing, spacing between (rows) of SAs, multiples of lambda/2 (1, 2, ...)
% DeltaN: SAs Spacing, spacing between (columns) of SAs, multiples of lambda/2 (1, 2, ...)
% deltaM: AEs Spacing, spacing between (rows) of AEs
% deltaN: AEs Spacing, spacing between (columns) of AEs
% Output Arguments:
% Object: this struct contains overall antenna dimenions, spacings, maximum array size (usefull for near-field/far-field models), etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Object = get_AoSA_param(p, M, N, Mae, Nae, DeltaM, DeltaN, deltaM, deltaN)
% Object could be BS, UE

ant_spac = p.lambda_c0/2;                              % Always the spacing is multiples of lambda/2
% # of SAs
Object.Qdim = [M N];                                       % Object SA paramters inside the array: [# of SAs (row), # of SAs (column)] [M N]
Object.Q = prod(Object.Qdim);                       % Object total # of SAs
% # of AEs
Object.Qbardim = [Mae Nae];                         % Object AE paramters inside each SA: [# of AEs (row), # of AEs (column)] [Mae Nae]
Object.Qbar = prod(Object.Qbardim);            % Object total # of AEs inside each SA
% SAs spacing
Object.Delta = [DeltaM DeltaN]*ant_spac;     % Object spacing between [rows, columns] of SAs
% AEs spacing
Object.delta = [deltaM deltaN]*ant_spac;      % Object spacing between [rows, columns] of AEs
Object.Qtot = Object.Q*Object.Qbar;

% Maximum Array Dimension, here we ignore the physical dimension of each AE
Object.Dmax = max(((Object.Qdim-1).*Object.Delta)+((Object.Qbardim-1).*Object.delta));
end
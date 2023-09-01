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
% Quantized_Phase_Shifts: The actual quantized phase shift values
% Q_quant: The number of phase quantization bits
% Q: The number of antennas (AEs in a SA)
% M: The number of desired measurements (codewords/beams)
% Output Arguments:
% PZ: The random RF Beamforming/Combining training codebook (P or Z used in Eq. 7 in the manuscript)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PZ = get_RandomCodebook(Quantized_Phase_Shifts,Q_quant,Q,M)
Rand_Indx = randi(2^Q_quant,[Q M]);
PZ = sqrt(1/Q)*exp(1j*Quantized_Phase_Shifts(Rand_Indx));
end
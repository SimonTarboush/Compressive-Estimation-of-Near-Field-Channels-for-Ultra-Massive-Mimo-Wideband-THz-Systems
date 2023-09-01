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
% This function computes the LoS path loss. This function supports the model: Hybrid Spherical Planar Wave Model (HSPM) 
%         i.e. Spherical Wave Model (SWM)/Plane Wave Model (PWM) (SWM/PWM) on the level of SA/AE.
% The other models are implemented in the original code of TeraMIMO
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% indx_subc: Index of the subcarrier
% d: distance between Tx-Rx SAs in LoS case (for HSPM: the input refers to the distance between SAs)
% K_abs: Molecular absorption coefficient
% Output Arguments:
% AlphaLoS: Path loss for the LoS component; the size is 3D-Array of size(1,1, num_freqs_per_subcarr)
%                   based on the used model for the SA/AE, i.e., HSPM (SWM/PWM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AlphaLoS = get_PathLoss(p,indx_subc,d,K_abs)
num_freqs_per_subcarr = p.nFreq(2);
D_eff = d; % Distance of SAs/AEs (based on model) between Rx & Tx
Tau_LoS = D_eff./p.c; % Deterministic delay of the LoS path
if strcmp(p.WaveModelAE,'Plane') 
    AlphaLoSAmp = zeros(1, 1, num_freqs_per_subcarr);
    AlphaLoSPhase = zeros(1, 1, num_freqs_per_subcarr);
    AlphaLoSAmp(1,1,:) = (p.c./(4*pi*p.freq(indx_subc,:).*D_eff)).^(p.PLE/2).*exp(-K_abs(indx_subc,:).*D_eff/2);
    AlphaLoSPhase(1,1,:) = exp(-1j*2*pi*p.freq(indx_subc,:).*Tau_LoS);
elseif strcmp(p.WaveModelAE,'Sphere')
    % Not implemented here, go to TeraMIMO original codes
else
    error('This wavemodel for AE level isn''t Implemented. Supported options are: SWM/PWM');
end
% Get path loss for LOS
% AlphaLoS = SpreadingLoss.*AbsorptionLoss .* AlphaLoSPhase
% AlphaLoS = AlphaLoSAmp .* AlphaLoSPhase
AlphaLoS = AlphaLoSAmp.*AlphaLoSPhase;
end
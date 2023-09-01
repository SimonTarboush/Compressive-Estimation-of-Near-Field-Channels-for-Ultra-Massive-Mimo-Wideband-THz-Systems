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
% This function computes the distances and local AoD/AoA. This function supports the model: Hybrid Spherical Planar Wave Model (HSPM) 
%         i.e. Spherical Wave Model (SWM)/Plane Wave Model (PWM) (SWM/PWM) on the level of SA/AE.
% The other models are implemented in the original code of TeraMIMO
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% mr: Index of the number of the receiver SAs (rows)
% nr:  Index of the number of the receiver SAs (columns)
% mt: Index of the number of the transmitter SAs (rows)
% nt:  Index of the number of the transmitter SAs (columns)
% Output Arguments:
% Dist_SAAEs: Distance between Tx-Rx SAs/AEs (for HSPM: we need the distance between SAs)
% AoD_SV_LOS_loc: Local angle-of-departure for LoS beamsteering vector; size is based on the used model for 
%                               the SA/AE, here is HSPM (SWM-SA/PWM-AE): vector of size(2, 1), (Azimuth; Elevation)
% AoA_SV_LOS_loc: Local angle-of-arrival for LoS beamsteering vector; size is based on the used model for 
%                               the SA/AE, here is HSPM (SWM-SA/PWM-AE): vector of size(2, 1), (Azimuth; Elevation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dist_SAAEs, AoD_SV_LOS_loc, AoA_SV_LOS_loc] = get_Distance_Angle_LoS(p,mr,nr,mt,nt)

positionLocalSARx = [0; (nr-1-(p.Rx_AoSA.Qdim(2)-1)/2)*p.Rx_AoSA.Delta(2); (mr-1-(p.Rx_AoSA.Qdim(1)-1)/2)*p.Rx_AoSA.Delta(1)];
positionLocalSATx = [0; (nt-1-(p.Tx_AoSA.Qdim(2)-1)/2)*p.Tx_AoSA.Delta(2); (mt-1-(p.Tx_AoSA.Qdim(1)-1)/2)*p.Tx_AoSA.Delta(1)];

if strcmp(p.WaveModelSA,'Sphere') && strcmp(p.WaveModelAE,'Plane') 
    p.RxSATxSA.Rot_mat = p.Rx.Rot_mat;
    p.TxSARxSA.Rot_mat = p.Tx.Rot_mat;
    p.RxSATxSA.Pos = p.RxSATxSA.Rot_mat*positionLocalSARx + p.Rx.Pos;
    p.TxSARxSA.Pos = p.TxSARxSA.Rot_mat*positionLocalSATx + p.Tx.Pos;
    p.RxSATxSAt = get_unitdirvec_glob_loc(p.RxSATxSA,p.TxSARxSA);
    p.TxSARxSAt = get_unitdirvec_glob_loc(p.TxSARxSA,p.RxSATxSA);
    % Physical AoD/AoA
    AoA_SV_LOS_loc = get_anglevec_from_unitdirvec(p.RxSATxSAt.t_Loc);
    AoD_SV_LOS_loc = get_anglevec_from_unitdirvec(p.TxSARxSAt.t_Loc);
    Dist_SAAEs = norm(p.RxSATxSA.Pos-p.TxSARxSA.Pos);
elseif strcmp(p.WaveModelSA,'Plane') && strcmp(p.WaveModelAE,'Plane')
    % See TeraMIMO codes
elseif strcmp(p.WaveModelSA,'Sphere') && strcmp(p.WaveModelAE,'Sphere')
    % See TeraMIMO codes
else
    error('This wavemodel on SA/AE level isn''t implemented. Supported options are: SWM/SWM, SWM/PWM, PWM/PWM');
end
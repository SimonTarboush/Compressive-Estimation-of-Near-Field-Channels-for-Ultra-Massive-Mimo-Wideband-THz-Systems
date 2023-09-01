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
% This function generates a time-invariant THz channel
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters.
% K_abs: Molecular absorption coefficient, a matrix of size(p.Nsub_c, p.Nsub_b),
%              (number of subcarriers, number of subbands @ each subcarrier)
% Output Arguments:
% CH_Response: A struct contains the channel response, CH_Response.H: H(f) time-invariant frequency domain response, 
%                          a 3D-array of size (p.Qr, p.Qt, num_subcarries), i.e. H(MIMO Rx SA, MIMO Tx SA, number of subcarries)
% AT_PCI: The Tx array response vectors (Tx-ARVs)
% AR_PCI: The Rx array response vectors (Rx-ARVs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CH_Response, AT_PCI, AR_PCI] = channel_TIV(p, K_abs)
% Initialize parameters and channnel matrix
num_subcarries = p.nFreq(1);
num_freqs_per_subcarr = p.nFreq(2);
H_LoS = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
AT_PCI = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
AR_PCI = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
% Main Loop
% Loop over subcarries
for indx_subc = 1:num_subcarries
    % Loops over SA MIMO arrays    
    for mr = 1:p.Rx_AoSA.Qdim(1)
        for nr = 1:p.Rx_AoSA.Qdim(2)
            for mt = 1:p.Tx_AoSA.Qdim(1)
                for nt = 1:p.Tx_AoSA.Qdim(2)
                    [d_SAAEs, AoD_SV_LOS_loc, AoA_SV_LOS_loc] = get_Distance_Angle_LoS(p,mr,nr,mt,nt);
                    AlphaLoS = get_PathLoss(p, indx_subc, d_SAAEs, K_abs);
                    if strcmp(p.channelType,'LoS')
                        switch p.AntennaModel
                            % Frequency Independent Model
                            case 'Isotropic'
                                    Gt_LOS = p.Tx.Gain;
                                    Gr_LOS = p.Rx.Gain;
                            otherwise
                                error('This Antenna Model isn''t implemented !!');
                        end
                        if strcmp(p.WaveModelAE,'Plane')
                            a_sv_t_LoS = get_ArrayResponse(p,indx_subc,p.Tx_AoSA.Qbardim(1),p.Tx_AoSA.Qbardim(2),AoD_SV_LOS_loc,'SV','T');
                            a_sv_t_LoS1 = a_sv_t_LoS;
                            a_sv_t_LoS_tmp = reshape(a_sv_t_LoS1,p.Tx_AoSA.Qbar,1,num_freqs_per_subcarr);                            
                            a_sv_r_LoS = get_ArrayResponse(p,indx_subc,p.Rx_AoSA.Qbardim(1),p.Rx_AoSA.Qbardim(2),AoA_SV_LOS_loc,'SV','R');
                            a_sv_r_LoS1 = a_sv_r_LoS;
                            a_sv_r_LoS_tmp = reshape(a_sv_r_LoS1,p.Rx_AoSA.Qbar,1,num_freqs_per_subcarr); 
                            H_LoS_tmp = zeros(p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar, num_freqs_per_subcarr);
                            for indx_numfreqpersubc = 1:num_freqs_per_subcarr
                                H_LoS_tmp(:,:,indx_numfreqpersubc) = Gr_LOS*Gt_LOS*AlphaLoS(:,:,indx_numfreqpersubc)*...
                                    a_sv_r_LoS_tmp(:,:,indx_numfreqpersubc)*a_sv_t_LoS_tmp(:,:,indx_numfreqpersubc).';
                            end
                            % LoS Channel, Frequency Domain Implementation
                            H_LoS{(mr-1)*p.Rx_AoSA.Qdim(2)+nr,(mt-1)*p.Tx_AoSA.Qdim(2)+nt,indx_subc} = sqrt(p.Rx_AoSA.Qbar*p.Tx_AoSA.Qbar)*H_LoS_tmp;
                            AT_PCI{(mr-1)*p.Rx_AoSA.Qdim(2)+nr,(mt-1)*p.Tx_AoSA.Qdim(2)+nt,indx_subc} = a_sv_t_LoS_tmp;
                            AR_PCI{(mr-1)*p.Rx_AoSA.Qdim(2)+nr,(mt-1)*p.Tx_AoSA.Qdim(2)+nt,indx_subc} = a_sv_r_LoS_tmp;
                        elseif strcmp(p.WaveModelAE,'Sphere')
                            % Not implemented here, go to TeraMIMO original codes
                        else
                            % AE SWM need to add SV and/or BF
                            error('This wavemodel for AE level isn''t implemented. Supported Options are: SWM/PWM');
                        end
                    end
                    
                end
            end
        end
    end
end
% Channel Response
switch p.channelType
    case 'LoS'
        CH_Response.H = H_LoS;
    otherwise
        error('This channel configuration is still not implemented in this version');
end
end
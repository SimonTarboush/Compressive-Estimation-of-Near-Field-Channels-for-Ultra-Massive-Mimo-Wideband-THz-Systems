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
% This function computes the number of best beams that form 95% of the total power for the estimated channel
% This can be used as indicator for the number of grid points that contain channel components in the beamspace domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% H_Est: The estimation of the channel
% Support_H: The index selected by the compressed sensing algorithm
% Mthd: Computation compatible with OMP/OMP-DR or SOMP/SOMP-DR
% K: Number of subcarriers
% Output Arguments:
% N_BestBeams: Number of best beams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N_BestBeams = get_NumBestBeams(H_Est, Support_H, Mthd, K)
if strcmp(Mthd,'SOMP') || strcmp(Mthd,'SOMP-DR')
    H_Est_TotalPow = sum(abs(H_Est).^2,'all');
    H_Est_95percentPow = 0;
    counter_beams = 0;
    while H_Est_95percentPow < 0.95*H_Est_TotalPow
        counter_beams = counter_beams + 1;
        H_Est_95percentPow = sum(abs(H_Est(Support_H(1:counter_beams),:)).^2,'all');
    end
    N_BestBeams = counter_beams;
else % OMP, OMP-DR 
    H_Est_TotalPow = sum(abs(H_Est).^2);
    H_Est_95percentPow = zeros(1,K);
    N_BestBeams = zeros(1,K);
    for indx_subc = 1:K
        counter_beams = 0;
        while H_Est_95percentPow(1,indx_subc) < 0.95*H_Est_TotalPow(1,indx_subc)
            counter_beams = counter_beams + 1;
            H_Est_95percentPow(1,indx_subc) = sum(abs(H_Est(Support_H(1:counter_beams,indx_subc),indx_subc)).^2);
        end
        N_BestBeams(1,indx_subc)= counter_beams;
    end
end
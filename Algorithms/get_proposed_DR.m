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
% This script implements the proposed dictionary reduction (DR) method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Based on the estimation of the first Tx/Rx SAs channel we will use the prior information to reduce the search region for the Dictionary.
% This is done for evey Tx SA.
if indx_UsedSARx ==2
    % OMP-DR
    Abar_T_OMP_DR = complex(zeros([size(Abar_T),K]));
    Abar_R_OMP_DR = complex(zeros([size(Abar_R),K]));
    Thetabar_OMP_DR = complex(zeros([size(Thetabar),K]));
    for indx_subc = 1:K
        % Detect the best candidate support based on the first estimation (we add -1 to the index of Rx-plane)
        Support_Tx_Est_OMP_DR = AoDInd_OMP_DR_cand{(qrv-1)*Q_R_h+qrh-1,(qtv-1)*Q_T_h+qth,indx_subc};
        Support_Rx_Est_OMP_DR = AoAInd_OMP_DR_cand{(qrv-1)*Q_R_h+qrh-1,(qtv-1)*Q_T_h+qth,indx_subc};
        % Define the new index for the Tx plane
        IndexTx_OMP_DR = (G_hv_TR_OVS.^2*Support_Tx_Est_OMP_DR(1)-G_T/2:G_hv_TR_OVS.^2*Support_Tx_Est_OMP_DR(1)+G_T/2-1);
        IndexTx_OMP_DR = IndexTx_OMP_DR(IndexTx_OMP_DR>0);
        IndexTx_OMP_DR = IndexTx_OMP_DR(IndexTx_OMP_DR<size(Abar_T_OVS,2));
        % Construct the new verion of Tx Dictionary based on DR method (\bar{A}_{T,DR} in the manuscript) by selecting the
        % columns from the Tx oversampled Dictionary based on the first Tx-Rx SA support (AoD) compatible with OMP
        Abar_T_OMP_DR(:,1:length(IndexTx_OMP_DR),indx_subc) = Abar_T_OVS(:,IndexTx_OMP_DR);
        % Define the new index for the Rx plane
        IndexRx_OMP_DR = (G_hv_TR_OVS.^2*Support_Rx_Est_OMP_DR(1)-G_R/2:G_hv_TR_OVS.^2*Support_Rx_Est_OMP_DR(1)+G_R/2-1);
        IndexRx_OMP_DR = IndexRx_OMP_DR(IndexRx_OMP_DR>0);
        IndexRx_OMP_DR = IndexRx_OMP_DR(IndexRx_OMP_DR<size(Abar_R_OVS,2));
        % Construct the new verion of Tx Dictionary based on DR method (\bar{A}_{R,DR} in the manuscript) by selecting the
        % columns from the Rx oversampled Dictionary based on the first Tx-Rx SA support (AoA) compatible with OMP
        Abar_R_OMP_DR(:,1:length(IndexRx_OMP_DR),indx_subc)  = Abar_R_OVS(:,IndexRx_OMP_DR);
        % Compute the new quantized Dictionary matrix based on DR method (\bar{Theta}_{DR} in the manuscript) compatible with OMP
        Thetabar_OMP_DR(:,:,indx_subc)  = kron(Abar_T_OMP_DR(:,:,indx_subc),Abar_R_OMP_DR(:,:,indx_subc));
    end
    % SOMP-DR
    % Detect the best candidate support based on the first estimation (we add -1 to the index of Rx-plane)
    Support_Tx_Est_SOMP_DR = AoDInd_SOMP_DR_cand{(qrv-1)*Q_R_h+qrh-1,(qtv-1)*Q_T_h+qth};
    Support_Rx_Est_SOMP_DR = AoAInd_SOMP_DR_cand{(qrv-1)*Q_R_h+qrh-1,(qtv-1)*Q_T_h+qth};
    % Define the new index for the Tx plane
    IndexTx_SOMP_DR = (G_hv_TR_OVS.^2*Support_Tx_Est_SOMP_DR(1)-G_T/2:G_hv_TR_OVS.^2*Support_Tx_Est_SOMP_DR(1)+G_T/2-1);
    IndexTx_SOMP_DR = IndexTx_SOMP_DR(IndexTx_SOMP_DR>0);
    IndexTx_SOMP_DR = IndexTx_SOMP_DR(IndexTx_SOMP_DR<size(Abar_T_OVS,2));
    % Construct the new verion of Tx Dictionary based on DR method (\bar{A}_{T,DR} in the manuscript) by selecting the
    % columns from the Tx oversampled Dictionary based on the first Tx-Rx SA support (AoD) compatible with SOMP
    Abar_T_SOMP_DR = Abar_T_OVS(:,IndexTx_SOMP_DR);
    % Define the new index for the Rx plane
    IndexRx_SOMP_DR = (G_hv_TR_OVS.^2*Support_Rx_Est_SOMP_DR(1)-G_R/2:G_hv_TR_OVS.^2*Support_Rx_Est_SOMP_DR(1)+G_R/2-1);
    IndexRx_SOMP_DR = IndexRx_SOMP_DR(IndexRx_SOMP_DR>0);
    IndexRx_SOMP_DR = IndexRx_SOMP_DR(IndexRx_SOMP_DR<size(Abar_R_OVS,2));
    % Construct the new verion of Tx Dictionary based on DR method (\bar{A}_{R,DR} in the manuscript) by selecting the
    % columns from the Rx oversampled Dictionary based on the first Tx-Rx SA support (AoA) compatible with SOMP
    Abar_R_SOMP_DR = Abar_R_OVS(:,IndexRx_SOMP_DR);
    % Compute the new quantized Dictionary matrix based on DR method (\bar{Theta}_{DR} in the manuscript) compatible with SOMP
    Thetabar_SOMP_DR = kron(Abar_T_SOMP_DR,Abar_R_SOMP_DR);
end

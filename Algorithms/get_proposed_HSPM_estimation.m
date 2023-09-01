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
% This script performs the compressed sensing channel estimation based on OMP, SOMP, and proposed dictionary reduction (DR) method combined with OMP/SOMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Arrays/Cells to save UM-MIMO channel estimation, support index, beamforming/combining vectors
% 
% Reference and Estimated Channels
H_SP_UM = cell(Q_R,Q_T,K);
H_Est_OMP_UM = cell(Q_R,Q_T,K);H_Est_OMP_DR_UM = cell(Q_R,Q_T,K);
H_Est_SOMP_UM = cell(Q_R,Q_T,K);H_Est_SOMP_DR_UM = cell(Q_R,Q_T,K);
% Best N Beams
BestBeams_OMP = zeros(Q_R,Q_T,K);BestBeams_OMP_DR = zeros(Q_R,Q_T,K);
BestBeams_SOMP = zeros(Q_R,Q_T);BestBeams_SOMP_DR = zeros(Q_R,Q_T);
% Support (AoD and AoA indices)
AoDInd_OMP_cand = cell(Q_R,Q_T,K);AoAInd_OMP_cand = cell(Q_R,Q_T,K);
AoDInd_OMP_DR_cand = cell(Q_R,Q_T,K);AoAInd_OMP_DR_cand = cell(Q_R,Q_T,K);
AoDInd_SOMP_cand = cell(Q_R,Q_T);AoAInd_SOMP_cand = cell(Q_R,Q_T);
AoDInd_SOMP_DR_cand = cell(Q_R,Q_T);AoAInd_SOMP_DR_cand = cell(Q_R,Q_T);
% Beamforming/Combining Vectors
F_Sel_OMP = cell(Q_R,Q_T,K);W_Sel_OMP = cell(Q_R,Q_T,K);
F_Sel_OMP_DR = cell(Q_R,Q_T,K);W_Sel_OMP_DR = cell(Q_R,Q_T,K);
F_Sel_SOMP = cell(Q_R,Q_T);W_Sel_SOMP = cell(Q_R,Q_T);
F_Sel_SOMP_DR = cell(Q_R,Q_T);W_Sel_SOMP_DR = cell(Q_R,Q_T);
%%%%%%%%%%%%%%%%%% Begin loops over Tx/Rx SAs %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% loop over the Tx SAs %%%%%%%%%%%%%%%%%%
for indx_UsedSATx = 1:Num_of_UsedSATx
    %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
    [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
    %%%%%%%%%%%%%%%%%% loop over the Rx SAs %%%%%%%%%%%%%%%%%%
    for indx_UsedSARx = 1:Num_of_UsedSARx
        %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
        [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
        %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
        P = get_RandomCodebook(Quant_PS_T, Q_T_quant, Qbar_T, Qbar_T); % Random beam training
        %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
        Z = get_RandomCodebook(Quant_PS_R, Q_R_quant, Qbar_R, Qbar_R); % Random beam training
        %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%%
        if isvector(M_T_measMAT)
            M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth);
            M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh);
        else
            M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth,Tx_meas_indx);
            M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh, Rx_meas_indx);
        end
        %%%%%%%%%%%%%%%%%% Proposed Algorithm for HSPM, Dictionary Reduction %%%%%%%%%%%%%%%%%%
        get_proposed_DR;
        %%%%%%%%%%%%%%%%%% WB channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
        H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
        % This is H_{q_R,q_T} defined in Eq. 3 in the manuscript
        %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
        N = sqrt(K)*sqrt(noise_pwr)/sqrt(2)*(randn(Qbar_R,Qbar_T,K)+1j*randn(Qbar_R,Qbar_T,K));
        %%%%%%%%%%%%%%%%%% Tx and Rx Signal %%%%%%%%%%%%%%%%%%
        y = get_channel_output(H_qRqT,N(:,1:M_T_meas,:),P(:,1:M_T_meas),Z(:,1:M_R_meas),K,Tx_pwr);
        %%%%%%%%%%%%%%%%%% Compute Measurement Matrix %%%%%%%%%%%%%%%%%%
        Psi = kron(transpose(P(:,1:M_T_meas)),(Z(:,1:M_R_meas))');
        % following Eq. 8 in the manuscript Psi = kron(P^Transpose, Z^Herm) i.e. Ptrans_kron_Zherm is the measurement matrix
        %%%%%%%%%%%%%%%%%% Compute Sensing Matrix %%%%%%%%%%%%%%%%%%
        % The Sensing Matrix will be denoted as Epsilon = Psi*Thetabar (Measurement Matrix * Quantized Dictionary Matrix)
        % OMP
        Epsilon_OMP = Psi*Thetabar;
        % SOMP
        Epsilon_SOMP = Epsilon_OMP;
        % OMP DR
        Epsilon_OMP_DR = complex(zeros([size(Epsilon_OMP),K]));
        if indx_UsedSARx ==1
            for indx_subc = 1:K
                Epsilon_OMP_DR(:,:,indx_subc) = Epsilon_OMP;
            end
            % SOMP DR
            Epsilon_SOMP_DR = Epsilon_SOMP;
        else
            for indx_subc = 1:K
                Epsilon_OMP_DR(:,:,indx_subc) = Psi*Thetabar_OMP_DR(:,:,indx_subc);
            end
            Epsilon_SOMP_DR = Psi*Thetabar_SOMP_DR;
        end
        %%%%%%%%%%%%%%%%%% OMP CS-Estimator %%%%%%%%%%%%%%%%%%
        % bs here stands for BeamSpace domain
        [hest_bs_OMP, support_OMP] = OMP_SMC(y,Epsilon_OMP,Tx_pwr,Lbar,K);
        H_bs_OMP = reshape(hest_bs_OMP,size(Abar_R,2),size(Abar_T,2),K);
        %%%%%%%%%%%%%%%%%% (OMP + Proposed DR) CS-Estimator %%%%%%%%%%%%%%%%%%
        [hest_bs_OMP_DR, support_OMP_DR] = OMP_SMC_DR(y,Epsilon_OMP_DR,Tx_pwr,Lbar,K);
        if indx_UsedSARx ==1
            H_bs_OMP_DR = reshape(hest_bs_OMP_DR,size(Abar_R,2),size(Abar_T,2),K);
        else
            H_bs_OMP_DR = reshape(hest_bs_OMP_DR,size(Abar_R_OMP_DR,2),size(Abar_T_OMP_DR,2),K);
        end        
        %%%%%%%%%%%%%%%%%% SOMP CS-Estimator %%%%%%%%%%%%%%%%%%
        [hest_bs_SOMP, support_SOMP] = SOMP(y,Epsilon_SOMP,Tx_pwr,Lbar,K);
        H_bs_SOMP = reshape(hest_bs_SOMP,size(Abar_R,2),size(Abar_T,2),K);        
        %%%%%%%%%%%%%%%%%% (SOMP + Proposed DR) CS-Estimator %%%%%%%%%%%%%%%%%%
        [hest_bs_SOMP_DR, support_SOMP_DR] = SOMP(y,Epsilon_SOMP_DR,Tx_pwr,Lbar,K);
        if indx_UsedSARx ==1
            H_bs_SOMP_DR = reshape(hest_bs_SOMP_DR,size(Abar_R,2),size(Abar_T,2),K);
        else
            H_bs_SOMP_DR = reshape(hest_bs_SOMP_DR,size(Abar_R_SOMP_DR,2),size(Abar_T_SOMP_DR,2),K);
        end
        %%%%%%%%%%%%%%%%%% Extract the number of dominant/best beams %%%%%%%%%%%%%%%%%%
        % Accumulate until you get 95% of the estimated channel power
        Bestn_Beams_OMP = get_NumBestBeams(hest_bs_OMP, support_OMP,'OMP',K);
        BestBeams_OMP((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,:) = Bestn_Beams_OMP;
        Bestn_Beams_OMP_DR = get_NumBestBeams(hest_bs_OMP_DR, support_OMP_DR,'OMP-DR',K);
        BestBeams_OMP_DR((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,:) = Bestn_Beams_OMP_DR;
        Bestn_Beams_SOMP = get_NumBestBeams(hest_bs_SOMP, support_SOMP,'SOMP',K);
        BestBeams_SOMP((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth) = Bestn_Beams_SOMP;
        Bestn_Beams_SOMP_DR = get_NumBestBeams(hest_bs_SOMP_DR, support_SOMP_DR,'SOMP-DR',K);
        BestBeams_SOMP_DR((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth) = Bestn_Beams_SOMP_DR;
        %%%%%%%%%%%%%%%%%% Extract the support and select angles for beamforming/combining %%%%%%%%%%%%%%%%%%
        % OMP  & OMP-DR
        for indx_subc = 1:K
            [rowMaxBeamRx_OMP, colMaxBeamTx_OMP] = ind2sub(size(H_bs_OMP), support_OMP(1:Bestn_Beams_OMP(1,indx_subc),indx_subc));
            [rowMaxBeamRx_OMP_DR, colMaxBeamTx_OMP_DR] = ind2sub(size(H_bs_OMP_DR), support_OMP_DR(1:Bestn_Beams_OMP_DR(1,indx_subc),indx_subc));
            AoDInd_OMP_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = colMaxBeamTx_OMP;
            AoAInd_OMP_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = rowMaxBeamRx_OMP;
            AoDInd_OMP_DR_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = colMaxBeamTx_OMP_DR;
            AoAInd_OMP_DR_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = rowMaxBeamRx_OMP_DR;
            F_Sel_OMP{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = Abar_T(:,colMaxBeamTx_OMP(1));
            W_Sel_OMP{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = Abar_R(:,rowMaxBeamRx_OMP(1));
            if indx_UsedSARx ==1
                F_Sel_OMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = Abar_T(:,colMaxBeamTx_OMP_DR(1));
                W_Sel_OMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = Abar_R(:,rowMaxBeamRx_OMP_DR(1));
            else
                F_Sel_OMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = Abar_T_OMP_DR(:,colMaxBeamTx_OMP_DR(1),indx_subc);
                W_Sel_OMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = Abar_R_OMP_DR(:,rowMaxBeamRx_OMP_DR(1),indx_subc);
            end
        end
        % SOMP  & SOMP-DR
        [rowMaxBeamRx_SOMP, colMaxBeamTx_SOMP] = ind2sub(size(H_bs_SOMP), support_SOMP(1:Bestn_Beams_SOMP));
        [rowMaxBeamRx_SOMP_DR, colMaxBeamTx_SOMP_DR] = ind2sub(size(H_bs_SOMP_DR), support_SOMP_DR(1:Bestn_Beams_SOMP_DR));
        AoDInd_SOMP_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = colMaxBeamTx_SOMP;
        AoAInd_SOMP_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = rowMaxBeamRx_SOMP;
        AoDInd_SOMP_DR_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = colMaxBeamTx_SOMP_DR;
        AoAInd_SOMP_DR_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = rowMaxBeamRx_SOMP_DR;        
        F_Sel_SOMP{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T(:,colMaxBeamTx_SOMP(1));
        W_Sel_SOMP{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R(:,rowMaxBeamRx_SOMP(1));
        if indx_UsedSARx ==1
            F_Sel_SOMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T(:,colMaxBeamTx_SOMP_DR(1));
            W_Sel_SOMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R(:,rowMaxBeamRx_SOMP_DR(1));
        else
            F_Sel_SOMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T_SOMP_DR(:,colMaxBeamTx_SOMP_DR(1));
            W_Sel_SOMP_DR{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R_SOMP_DR(:,rowMaxBeamRx_SOMP_DR(1));
        end
        %%%%%%%%%%%%%%%%%% Fill the SA (MIMO channel) estimation into the UM-MIMO channels  %%%%%%%%%%%%%%%%%%
        H_Est_OMP = zeros(Qbar_R,Qbar_T,K);
        H_Est_OMP_DR = zeros(Qbar_R,Qbar_T,K);
        H_Est_SOMP = zeros(Qbar_R,Qbar_T,K);
        H_Est_SOMP_DR = zeros(Qbar_R,Qbar_T,K);
        for indx_subc = 1:K
            % Create Angle-Domain Channels
            % Relation with beamspace domain representation is presented in Eq. 5 in the manuscript
            H_Est_OMP(:,:,indx_subc) = Abar_R*H_bs_OMP(:,:,indx_subc)*Abar_T.';
            if indx_UsedSARx ==1
                H_Est_OMP_DR(:,:,indx_subc) = Abar_R*H_bs_OMP_DR(:,:,indx_subc)*Abar_T.';
            else
                H_Est_OMP_DR(:,:,indx_subc) = Abar_R_OMP_DR(:,:,indx_subc)*H_bs_OMP_DR(:,:,indx_subc)*Abar_T_OMP_DR(:,:,indx_subc).';
            end
            H_Est_SOMP(:,:,indx_subc) = Abar_R*H_bs_SOMP(:,:,indx_subc)*Abar_T.';
            if indx_UsedSARx ==1
                H_Est_SOMP_DR(:,:,indx_subc) = Abar_R*H_bs_SOMP_DR(:,:,indx_subc)*Abar_T.';
            else
                H_Est_SOMP_DR(:,:,indx_subc) = Abar_R_SOMP_DR*H_bs_SOMP_DR(:,:,indx_subc)*Abar_T_SOMP_DR.';
            end
            % Fill the UM-MIMO arrays by estimations
            H_SP_UM{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_qRqT(:,:,indx_subc);
            H_Est_OMP_UM{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_OMP(:,:,indx_subc);
            H_Est_OMP_DR_UM{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_OMP_DR(:,:,indx_subc);
            H_Est_SOMP_UM{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP(:,:,indx_subc);
            H_Est_SOMP_DR_UM{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_DR(:,:,indx_subc);
        end
        %%%%%%%%%%%%%%%%%% End of loop on Tx-Rx SA channel Estimation  %%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%% End of loops over UM-MIMO channels  %%%%%%%%%%%%%%%%%%
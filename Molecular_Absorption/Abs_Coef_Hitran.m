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
% This function computes the exact molecular absorption coefficients in the THz band using the "HITRAN" database, valid in the frequency range: [0.1-10] THz
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% Output Arguments:
% K_abs: Molecular absorption coefficient, a matrix of size(p.Nsub_c, p.Nsub_b),
%              (number of subcarriers, number of subbands @ each subcarrier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main references used for molecular absorption coefficients calculations: 
%    Ref [1]: J. M. Jornet and I. F. Akyildiz, "Channel Modeling and Capacity Analysis for Electromagnetic Wireless Nanonetworks in the Terahertz Band,"
%                in IEEE Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211-3221, October 2011, doi: 10.1109/TWC.2011.081011.100545.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K_abs = Abs_Coef_Hitran(p)
K_abs = zeros(size(p.freq));
for indx_subc = 1:p.Nsub_c
    for indx_molc= 1:length(p.molecules)
        % Absorption coefficient for a specific gas mixture
        
        % load data for each molecule of the transmissionmedium, with a specific mixing ratio,
        % from HITRAN database.
        % For more Molecules, check "Molecular_Absorption" folder --->
        % "Data" subfolders, and edit the inputs from function
        % "generate_channel_param_TIV"
        Gas_Data = load_Gas_Data(p.molecules{indx_molc}, p.moleculesRatio(indx_molc), p.c);
        
        % Compute absorption coefficient for a specific gas
        K_Gas = Gas_Abs_Coef(p, Gas_Data, indx_subc);
        
        % Accumulate all gases to get the total absorption coefficient
        K_abs(indx_subc,:) = K_abs(indx_subc,:) + K_Gas;
        
    end
end
end
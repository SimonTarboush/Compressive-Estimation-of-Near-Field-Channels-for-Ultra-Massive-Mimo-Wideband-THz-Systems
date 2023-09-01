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
% This function loads the data of "HITRAN" database from the attached .csv file for a specific gas
% For more information about the "data" struture, please see "header.txt" file in "Data" subfolder
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% moleculeName: Gas name
% moleculeRatio: Ratio of a gas in the transmission medium
% c: Speed of light in vaccum
% Output Arguments:
% GasData: A struct that contains the information to compute the molecular absorption coefficient for a specific gas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main references used for molecular absorption coefficients calculations: 
%    Ref [1]: J. M. Jornet and I. F. Akyildiz, "Channel Modeling and Capacity Analysis for Electromagnetic Wireless Nanonetworks in the Terahertz Band,"
%                in IEEE Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211-3221, October 2011, doi: 10.1109/TWC.2011.081011.100545.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GasData = load_Gas_Data(moleculeName, moleculeRatio, c)

% Read .csv data file
M = csvread(sprintf('Data/%s.csv',moleculeName));

% Molecule percentage
GasData.Q = moleculeRatio;         %[ ]

% Zero-pressure position of the resonance frequency
GasData.fc0 = c*100*M(:,2);        %[Hz]

% Absorption strength using line intensity
GasData.S = c/100*M(:,3);          %[Hz.m^2/molecule]

% Linear pressure shift
GasData.delta = c*100*M(:,4);      %[Hz]

% Temperature coefficient
GasData.gamma = M(:,5);            %[ ]

% Air half-widths
GasData.alphaAir = c*100*M(:,6);   %[Hz]

% Self-broadened half-widths
GasData.alphaGas = c*100*M(:,7);   %[Hz]

% Mixing ratio for the isotopologue i of gas g
% gas ratio * isotopologue abundance
GasData.q = moleculeRatio*M(:,8);  %[ ]

end
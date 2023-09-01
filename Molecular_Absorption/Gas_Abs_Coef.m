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
% This function computes the exact molecular absorption coefficients in the THz band for a specific gas using the "HITRAN" database, valid in the frequency range: [0.1-10] THz. 
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% data: A struct contains the information to compute molecular absorption coefficient for a specific gas
% indx_sc: Index of a subcarrier
% Output Arguments:
% K_Gas: Molecular absorption coefficient for a specific gas, a vector of size(p.Nsub_b, 1), (number of subbands @ each subcarrier, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main references used for molecular absorption coefficients calculations: 
%    Ref [1]: J. M. Jornet and I. F. Akyildiz, "Channel Modeling and Capacity Analysis for Electromagnetic Wireless Nanonetworks in the Terahertz Band,"
%                in IEEE Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211-3221, October 2011, doi: 10.1109/TWC.2011.081011.100545.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K_Gas = Gas_Abs_Coef(p, data, indx_sc)

R = p.R/101325;                                     % Gas constant (m^3 atm/K/mol)
f = repmat(p.freq(indx_sc,:),length(data.fc0),1).';   % Frequency

% Total # of molecules (Eq.4 Ref. [1])
% The abundance is not accounted here
% It is accounted in line intensity data.S (Eq.4 Ref. [1])
Q = (p.p/(R*p.T))*data.Q*p.Na;

% Position of resonant f = zero_pos + linear_p shift (Eq.6 Ref. [1])
fc = repmat((data.fc0+data.delta*(p.p/p.p0))',p.nFreq(2),1);

% Lorentz half-width (Eq.7 Ref. [1])
% The first alpha is computed using q^{i,j},
% and the second one assumes that the abundance is already multiplied with
alpha = repmat((((1-data.q).*(data.alphaAir)+data.q.*(data.alphaGas)).*((p.p/p.p0)*((p.T0/p.T).^data.gamma)))',p.nFreq(2),1);

% Van Vleck-Weisskopf asymmetric line shape (NOTE THAT WE REMOVED THE TERM 100*p.c) (Eq.8 Ref. [1])
F = (1/pi) * alpha.*f./fc.* (1./((f-fc).^2+(alpha).^2) + 1./((f+fc).^2+(alpha).^2));

% Line shape (Eq.9 Ref. [1])
G = f./fc.* tanh((p.h*p.c*f)/(2*p.Kb*p.T)) ./ tanh((p.h*p.c*fc)/(2*p.Kb*p.T)) .*F;

% Absorption cross section = line intensity * line shape (Eq.5 Ref. [1])
sigma = repmat((data.S)',p.nFreq(2),1) .* G; 
% Note that in this version of the code, the Line Intensity is computed only @ T0 = 296 K

% Absorption coefficient (Eq.3 Ref. [1])
K_Gas = sum(p.p/p.p0*p.Tstp/p.T*Q*sigma, 2).';

end
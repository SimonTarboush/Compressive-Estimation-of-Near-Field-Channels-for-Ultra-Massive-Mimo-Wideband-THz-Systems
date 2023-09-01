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
% This function selects the calculation method of molecular absorption coefficient
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% "p.absorptionType" defines the method of computation
% Options are: Hitran, Approx1, and Approx2
%          1) Hitran: Exact absorption coefficient, valid in the frequency band: [0.1-10] THz
%          2) Approx1: First approximation of absorption coefficient, valid in the frequency band: [275-400] GHz
%          3) Approx2: Second approximation of absorption coefficient, valid in the frequency band: [100-450] GHz
% Output Arguments:
% K_abs: Molecular absorption coefficient, a matrix of size(p.Nsub_c, p.Nsub_b),
%              (number of subcarriers, number of subbands @ each subcarrier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K_abs = get_Abs_Coef(p)
% Select Computation Method
switch p.absorptionType
    case 'Hitran'
        % 1) Hitran: Exact absorption coefficient is valid in the frequency band: [0.1-10] THz
        K_abs = Abs_Coef_Hitran(p);
    case 'Approx1'
        % 2) Approx1: Approximation of absorption coefficient is valid in the frequency band: [275-400] GHz
        K_abs = Abs_Coef_Approx1(p);
    case 'Approx2'
        % 3) Approx2: Approximation of absorption coefficient is valid in the frequency band: [100-450] GHz
        K_abs = Abs_Coef_Approx2(p);
    otherwise
         error('This method for computing absorption coefficient isn''t implemented, options are: Hitran, Approx1, Approx2');
end
end
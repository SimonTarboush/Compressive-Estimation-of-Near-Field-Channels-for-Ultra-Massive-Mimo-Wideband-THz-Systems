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
% This function generates the required time-invariant THz channel parameters
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% No input
% Output Arguments:
% p: The channel struct that contains all of the required simulation parameters to generate the THz channel using TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p  = generate_channel_param_TIV()
%% Scenario Parameters
% Define channel type
p.channelType = 'LoS';
% Method for calucating the molecular absorption coefficient
% Options are: /'Hitran' /'Approx1' /'Approx2'
% These models represent Eq. (24), (30), and (34), respectively, in TeraMIMO IEEE_TVT
% Hitran: Exact absorption coefficient, valid in the frequency band: [0.1-10] THz
% Approx1: First approximation of absorption coefficient, valid in the frequency band: [275-400] GHz
% Approx2: Second approximation of absorption coefficient, valid in the frequency band: [100-450] GHz
p.absorptionType = 'Hitran';
% Spherical and planar wave model, HSPM
% keep these parameters as 'Sphere' for SA and 'Plane' for AE for this work
p.WaveModelSA = 'Sphere';  %'Sphere'/'Plane'
p.WaveModelAE = 'Plane';  %'Sphere'/'Plane'
% Antenna Model
p.AntennaModel = 'Isotropic'; % default: 'Isotropic', 'Directional'
%% Constants
p.c = physconst('LightSpeed');       % Speed of light in vacuum (m/sec)
p.T0 = 296;                                      % Reference temperature (Kelvin)
p.p0 = 1;                                          % Standard pressure (atm)
p.Kb = physconst('Boltzmann');     % Boltzmann constant (Joule/Kelvin)
% System Constants
p.T = 298.15;                                   % System temperature (Kelvin) [25°C Celsius]; note: 0 Kelvin - 273.15 = -273.1°C
p.p = 1;                                            % System pressure (atm)
p.PLE = 2;                                        % Path loss exponent
%% Transmission Parameters (Fc & BW)
p.Fc = 0.3e12;                       % Center frequency of the transmission bandwidth (Hz): f_c in the manuscript
p.BW = 0.03e12;                   % Total channel bandwidth (Hz): B_sys in the manuscript
p.Nsub_c = 4;                        % Number of subcarriers to divide the total Bandwidth: K-subcarriers in the manuscript
p.Nsub_b = 1;                        % Number of sub-bands in each subcarrier. In this work this value will be always equal to 1 do not change
p = get_FrequencyRng_param(p);
%% Absorption coefficient calculation
if strcmp(p.absorptionType,'Hitran')
    if (p.Fc - p.BW/2) < 0.1e12 || (p.Fc + p.BW/2) > 10e12
        error('Error: HITRAN is valid in: [0.1-10] THz');
    end
    p = get_HITRAN_param(p);
elseif strcmp(p.absorptionType,'Approx1')
    if (p.Fc - p.BW/2) < 0.275e12 || (p.Fc + p.BW/2) > 0.4e12
        error('Error: approx1, approximation of abs_coef is valid in: [275-400] GHz');
    end
    p.rel_humidity = 50;      % Relative humidity
elseif strcmp(p.absorptionType,'Approx2') 
    if (p.Fc - p.BW/2) < 0.1e12 || (p.Fc + p.BW/2) > 0.45e12
        error('Error: approx2, approximation of abs_coef is valid in: [100-450] GHz');
    end
    p.rel_humidity = 50;      % Relative humidity
else
    error('This method for computing absorption coefficient isn''t implemented, options are: Hitran, Approx1, Approx2');
end
%% UM-MIMO transceiver design (AoSA structure)
% Notations are following TeraMIMO IEEE_TVT paper

% This part includes the following definitions for each component 'i' in {BS, UE}:
% 1- Number of subarrays (SAs): total number is [Q_i = M_i*N_i where M_i (# of rows) and N_i (# of columns)]
% 2- Number of antenna elements (AEs) inside each SA: total number is [Qbar_i = Mbar_i*Nbar_i where Mbar_i (# of rows) and Nbar_i (# of columns)]
% 3- SA Spacing: Deltaxx is defined from center of SA(q) to center of SA(q+1) (default spacD*p.lambda_c/2, where spacD = sM x Mbar_i or sN x Nbar_i), we only define spacD
% 4- AE Spacing: deltaxx is defined from center of AE(qb) to center of AE(qb+1) (default spacd*p.lambda_c/2, spacd = 1), we fix this to half-wavelength

% Tx AoSA parameters
p.Mt = 2;
p.Nt = 2;
p.Mat = 2*4;
p.Nat = 2*4;
p.Sep_factTx = 9;
p.DeltaMt = p.Sep_factTx*p.Mat;
p.DeltaNt = p.Sep_factTx*p.Nat;
% Rx AoSA parameters
p.Mr = 2;
p.Nr = 2;
p.Mar = 2*4;
p.Nar = 2*4;
p.Sep_factRx = 8;
p.DeltaMr = p.Sep_factRx*p.Mar;
p.DeltaNr = p.Sep_factRx*p.Nar;
% Channel struct, M, N, Mbar, Nbar, DeltaM, DeltaN, DeltaMbar, DeltaNbar 
p.Tx_AoSA = get_AoSA_param(p, p.Mt, p.Nt, p.Mat, p.Nat, p.DeltaMt, p.DeltaNt, 1, 1);
p.Rx_AoSA = get_AoSA_param(p, p.Mr, p.Nr, p.Mar, p.Nar, p.DeltaMr, p.DeltaNr, 1, 1);
%% Design geometry
% Define local/global position and Euler angles

% Tx center 3D positions (global coordinates)
% Tx Euler rotation angles, following ZYX intrinsic rotation
% [Px Py Pz];[alphadot betadot gammadot] (degrees)
p.positionTx = [0; 0; 0];
p.eulerTx = [0; 0; 0];
p.Tx = get_Geometry_param(p.positionTx(1), p.positionTx(2), p.positionTx(3), p.eulerTx(1), p.eulerTx(2), p.eulerTx(3));

% Rx center 3D positions (global coordinates)
% Rx Euler rotation angles, following ZYX intrinsic rotation
p.positionRx = [0.15; 0; 0];
p.eulerRx = [rad2deg(pi); 0; 0]; % DO NOT CHANGE alpha_rot (the first angle), if you want the two UM arrays to face each other
p.Rx = get_Geometry_param(p.positionRx(1), p.positionRx(2), p.positionRx(3), p.eulerRx(1), p.eulerRx(2), p.eulerRx(3)); 
% to rotate Rx with +/- 10 deg over the azimuth use: 180 +/- 10

% Transmission distance (m)(from center to center)
p.d_tx_rx = norm(p.Tx.Pos-p.Rx.Pos);
p.d_rx_tx = norm(p.Rx.Pos-p.Tx.Pos);

% Check Rayleigh Distance for SWM/PWM
% to compute the far-field && near-field regions 
% Maximum Array size, here we ignore the physical dimension of each AE

p.D_max = max([p.Tx_AoSA.Dmax p.Rx_AoSA.Dmax]);
p.Dist_Fraunhofer_MISO_SIMO = 2*(p.D_max)^2/p.lambda_c0;
p.Dist_Fresnel_MISO_SIMO = 0.62*sqrt(p.D_max^3/p.lambda_c0);
p.Dist_Fraunhofer_MIMO = 2*(p.Tx_AoSA.Dmax+p.Rx_AoSA.Dmax)^2/p.lambda_c0;
p.Dist_Fraunhofer_RD_MIMO = 4*p.Tx_AoSA.Dmax*p.Rx_AoSA.Dmax/p.lambda_c0;
%% Tx, Rx antenna gain (frequency-independent)

switch p.AntennaModel
    case 'Isotropic'
        p.Tx.GaindBi = 0;    % Gain (dBi)
        p.Tx.Gain = 10.^(p.Tx.GaindBi/10); % Gain (scalar value)
        p.Rx.GaindBi = 0;    % Gain (dBi)
        p.Rx.Gain = 10.^(p.Rx.GaindBi/10); % Gain (scalar value)   
    case 'Directional'
        % Ref [1]: Constantine A. Balanis, "Antenna Theory: Analysis and Design" 4th edition Eq.(2-26)
        % assuming Directivity = Gain i.e. e_cd = e_c * e_d = 1
        % No losses, i.e. perfect conduction efficiency & perfect dielectric efficiency
        % Also, perfect reflection (mismatch) efficiency e_r (i.e. e_r=1)
        % Reflection coefficient equal to zero (matching between source & antenna)
        
        % half-power beamwidth in azimuth-plane, examples: rad2deg(2*sqrt(pi)) 27.7 60 120
        % half-power beamwidth in elevation-plane, examples: rad2deg(2*sqrt(pi)) 27.7 30 60
        % input value in (degree), output in (Rad)
        
        % Tx
        p.Tx.psi_azi = deg2rad(rad2deg(60));
        p.Tx.psi_elev = deg2rad(rad2deg(30));
        p.Tx.GaindBi = pow2db(4*pi/(p.Tx.psi_azi*p.Tx.psi_elev));    % Gain (dBi)
        % Eq. (51) in TeraMIMO IEEE_TVT
        p.Tx.Gain = 10.^(p.Tx.GaindBi/10);                           % Gain (scalar value)
        % Rx
        p.Rx.psi_azi = deg2rad(rad2deg(60));
        p.Rx.psi_elev = deg2rad(rad2deg(30));
        p.Rx.GaindBi = pow2db(4*pi/(p.Rx.psi_azi*p.Rx.psi_elev));    % Gain (dBi)
        % Eq. (51) in TeraMIMO IEEE_TVT
        p.Rx.Gain = 10.^(p.Rx.GaindBi/10);                           % Gain (scalar value)        
    otherwise
        error('This Antenna Model isn''t implemented !!');
end
end

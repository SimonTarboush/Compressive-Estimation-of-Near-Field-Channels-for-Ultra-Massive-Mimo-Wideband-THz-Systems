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
% This function generates the geometry parameters (position, rotation angles, rotation matrix)
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Px: Global position coordinate over x-axis
% Py: Global position coordinate over y-axis
% Pz: Global position coordinate over z-axis
% Euler_alphadot: Euler rotation angle around the Z-axis of $\dot{\alpha}$ degrees
% Euler_betadot: Euler rotation angle around the Y-axis of $\dot{\beta}$ degrees
% Euler_gammadot: Euler rotation angle around the X-axis of $\dot{\gamma}$ degrees
% Output Arguments:
% Object: struct contains AoSA's position vector, rotation angles vector, and rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Object = get_Geometry_param(Px, Py, Pz, Euler_alphadot, Euler_betadot, Euler_gammadot)

% Object could be BS, UE
% Define the geometry parameters
% Define local/global position and Euler angles
% Global position coordinates: center of 3D positions (size 3x1) for each instance of the object
% Euler rotation angles: vector (size 3x1) following ZYX intrinsic rotation (degrees)

Object.Pos = [Px; Py; Pz];

% Check angle ranges:
% Euler_alphadot in (-180,180]
% Euler_betadot in [-90,90]
% Euler_gammadot in (-180,180]
if Euler_alphadot > 180 || Euler_alphadot <= -180
    error('Alphadot is in range: (-180,180]');
end
if Euler_betadot > 90 || Euler_betadot < -90
    error('Betadot is in range: [-90,90]');
end
if Euler_gammadot > 180 || Euler_gammadot <= -180
    error('Gammadot is in range: (-180,180]');
end
Object.Euler = [Euler_alphadot; Euler_betadot; Euler_gammadot];

% Define Rotation matrix defined in Eq.(6) @ TeraMIMO IEEE_TVT
% Update geometry (global --> local)
% global = Rot_mat*local;
% local = (Rot_mat^-1)*global;

% Rotation matrix is real orthogonal (inv = transpose)
% local = (Rot_mat^T)*global;
Object.Rot_mat = eul2rotm(deg2rad(Object.Euler)', 'ZYX');

end
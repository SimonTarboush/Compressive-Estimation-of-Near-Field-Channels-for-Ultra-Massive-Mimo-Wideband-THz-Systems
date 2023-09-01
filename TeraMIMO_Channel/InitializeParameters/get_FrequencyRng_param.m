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
% This function generates the frequency, bandwidth, and wavelength parameters
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct contains the transmission parameters (Fc, BW, Nsub_c, Nsub_b) and other simulation parameters
% Output Arguments:
% p: The updated channel struct based on the transmission parameters such as the frequency points matrix, 
% bandwidth/first frequency/last frequency of each subcarrier, wavelengths related to center frequency/subcarrier/subbands, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = get_FrequencyRng_param(p)
% The following definitions and equations can be found in Sec.III-A of TeraMIMO paper, published in IEEE_TVT

% Divide the total bandwidth (BW) to a set of K-subcarriers
% then divide each subcarrier to a number of sub-bands

% How the total BW is divided
% |----------------------------| -> total BW
%   0   1                             K-1
% |---|---|----------------|----| -> total (K-subcarriers) (K=Nsub_c)
% Now the k-th subcarrier is divided into Nsub_b sub-bands
% (k-th subcarrier)->|-------------------------|
%                      0   1              nsub_b-1
%                    |---|---|-------------|---|  total (Nsub_b sub-bands)

% freq is a Matrix of size (K x Nsub_b) : K = Nsub_c

% p.freq represents the vector containing the frequencies of each k-th subcarrier
p.freq = zeros(p.Nsub_c,p.Nsub_b);
p.nFreq = size(p.freq);

p.BW_sub_c = p.BW/p.Nsub_c;        % Subcarrier bandwidth (Hz)
p.fstart_sub_c = p.Fc-p.BW/2+p.BW_sub_c/2;
p.fstop_sub_c = p.fstart_sub_c+p.BW-p.BW_sub_c/2;
p.Fc_sub_c = p.fstart_sub_c:p.BW_sub_c:p.fstop_sub_c;  % Center frequency of each subcarrier (Hz)

for indx_subc = 1:p.nFreq(1) % loop over subcarriers
    p.fstep_sub_b = p.BW_sub_c/p.Nsub_b;
    p.fstart_sub_b = p.Fc_sub_c(indx_subc)-p.BW_sub_c/2+p.fstep_sub_b/2;
    p.fstop_sub_b = p.fstart_sub_b+p.BW_sub_c-p.fstep_sub_b/2;
    p.freq (indx_subc,:)= (p.fstart_sub_b:p.fstep_sub_b:p.fstop_sub_b).';
end

% Compute lambdas
p.lambda_c0 = p.c./p.Fc;                    % Wavelength at center frequency (m)
p.lambda_kc = p.c./p.Fc_sub_c;         % Wavelength at center frequency of each subcarrier (m)
p.lambda_subb = p.c./p.freq;            % Wavelength for each sub-band center frequency inside every subcarrier (m)
end


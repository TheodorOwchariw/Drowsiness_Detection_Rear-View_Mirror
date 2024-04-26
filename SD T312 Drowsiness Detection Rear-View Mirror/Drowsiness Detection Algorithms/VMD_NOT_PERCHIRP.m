% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Performs further signal processing on lung data and compares 
%              it to inital head data values, this data will be used as a
%              percentage. Completion of respiration isolation was
%              unsuccessful, further implementation will be necesssary.

%% VMD for initial_respiration_data
% extract real part of the signal for the first channel
real_part = real(initial_respiration_data(:, 1));
% perform VMD on the real part of the signal with 6 IMFs
[real_imf, real_residual] = vmd(real_part, 'NumIMFs', 6);
% define time parameters
fs = 1e3;
t = (0:size(real_part, 1)-1) / fs;
% plot the real part of the signal and its VMD modes
 
%% FFT OF VMD for initial_respiration_data
% define the sampling frequency and time vector
fs = 2*10*10^(4); % assuming a sampling frequency of 1000 Hz
% number of data points
N = length(real_imf);
% frequency vector
f = (0:N-1)*(fs/N);
% perform FFT for each IMF
num_IMFs = size(real_imf, 2); % number of IMFs
% initialize arrays to store FFT results
fft_real_imf = zeros(N, num_IMFs);
for i = 1:num_IMFs
    % perform FFT for real part of IMFs
    fft_real_imf(:, i) = abs(fft(real_imf(:, i)));
end
 
 for i = 1:num_IMFs
     fft_length = size(fft_real_imf, 1);
 end
peak_Initial = max(fft_real_imf(1:fft_length/2, 5));
%% VMD for live_respiration_data
% extract real part of the signal for the first channel
real_part_head = real(live_respiration_data(:, 1));
% perform VMD on the real part of the signal with 6 IMFs
[real_imf_head, real_residual_head] = vmd(real_part_head, 'NumIMFs', 6);
% define time parameters
fs = 1e3;
t_head = (0:size(real_part_head, 1)-1) / fs;
% plot the real part of the signal and its VMD modes
 
%% FFT OF VMD for live_respiration_data
% define the sampling frequency and time vector
fs = 2*10*10^(4); % assuming a sampling frequency of 1000 Hz
% number of data points
N = length(real_imf_head);
% frequency vector
f = (0:N-1)*(fs/N);
% perform FFT for each IMF
num_IMFs = size(real_imf_head, 2); % number of IMFs
% initialize arrays to store FFT results
fft_real_imf_head = zeros(N, num_IMFs);
for i = 1:num_IMFs
    % perform FFT for real part of IMFs
    fft_real_imf_head(:, i) = abs(fft(real_imf_head(:, i)));
end

 for i = 1:num_IMFs
%  length of the FFT result
     fft_length = size(fft_real_imf_head, 1);
 end
peak_Live = max(fft_real_imf_head(1:fft_length/2, 5));
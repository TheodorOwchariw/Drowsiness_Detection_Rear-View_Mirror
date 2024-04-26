% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Performs further signal processing on lung data and compares 
%              it to inital head data values, this data will be used as a
%              percentage. Completion of respiration isolation was
%              unsuccessful, further implementation will be necesssary.

%% VMD for initialRespirationData
% Extract real part of the signal for the first channel
real_part = real(initialRespirationData(:, 1));
% Perform VMD on the real part of the signal with 6 IMFs
[real_imf, real_residual] = vmd(real_part, 'NumIMFs', 6);
% Define time parameters
fs = 1e3;
t = (0:size(real_part, 1)-1) / fs;
% Plot the real part of the signal and its VMD modes
 
%% FFT OF VMD for initialRespirationData
% Define the sampling frequency and time vector
fs = 2*10*10^(4); % Assuming a sampling frequency of 1000 Hz
% Number of data points
N = length(real_imf);
% Frequency vector
f = (0:N-1)*(fs/N);
% Perform FFT for each IMF
num_IMFs = size(real_imf, 2); % Number of IMFs
% Initialize arrays to store FFT results
fft_real_imf = zeros(N, num_IMFs);
for i = 1:num_IMFs
    % Perform FFT for real part of IMFs
    fft_real_imf(:, i) = abs(fft(real_imf(:, i)));
end
% Plot FFT of each IMF for real part
figure(12);
for i = 1:num_IMFs
    % Length of the FFT result
    fft_length = size(fft_real_imf, 1);
    % Frequency vector corresponding to the FFT result
    f_fft = (0:fft_length-1)*(fs/fft_length);
    % Plot FFT of real part of IMFs
    subplot(num_IMFs, 1, i);
    plot(f_fft(1:fft_length/2), fft_real_imf(1:fft_length/2, i));
    title(['FFT of IMF INIT ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end
peak_Initial = max(fft_real_imf(1:fft_length/2, 5));
%% VMD for liveRespirationData
% Extract real part of the signal for the first channel
real_part_head = real(liveRespirationData(:, 1));
% Perform VMD on the real part of the signal with 6 IMFs
[real_imf_head, real_residual_head] = vmd(real_part_head, 'NumIMFs', 6);
% Define time parameters
fs = 1e3;
t_head = (0:size(real_part_head, 1)-1) / fs;
% Plot the real part of the signal and its VMD modes
 
%% FFT OF VMD for liveRespirationData
% Define the sampling frequency and time vector
fs = 2*10*10^(4); % Assuming a sampling frequency of 1000 Hz
% Number of data points
N = length(real_imf_head);
% Frequency vector
f = (0:N-1)*(fs/N);
% Perform FFT for each IMF
num_IMFs = size(real_imf_head, 2); % Number of IMFs
% Initialize arrays to store FFT results
fft_real_imf_head = zeros(N, num_IMFs);
for i = 1:num_IMFs
    % Perform FFT for real part of IMFs
    fft_real_imf_head(:, i) = abs(fft(real_imf_head(:, i)));
end
% Plot FFT of each IMF for real part
figure(11);
for i = 1:num_IMFs
    % Length of the FFT result
    fft_length = size(fft_real_imf_head, 1);
    % Frequency vector corresponding to the FFT result
    f_fft = (0:fft_length-1)*(fs/fft_length);
    % Plot FFT of real part of IMFs
    subplot(num_IMFs, 1, i);
    plot(f_fft(1:fft_length/2), fft_real_imf_head(1:fft_length/2, i));
    title(['FFT of IMF LIVE', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end
peak_Live = max(fft_real_imf_head(1:fft_length/2, 5));


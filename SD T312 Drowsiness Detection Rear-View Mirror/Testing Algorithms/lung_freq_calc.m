% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Performs signal processing on lung data and compares it to
%              inital head data values, this data will be used as a
%              percentage. Completion of respiration isolation was
%              unsuccessful, further implementation will be necesssary.

initializedMxRawData = initialData;
headMxRawData = liveData;

% % 
%  load('sample6car.mat')
% %  headMxRawData = initialRespirationData;
%  headMxRawData = liveRespira tionData;
samplesPerChirp = 128;
%% CHANGE SIZE BEING TESTED HERE
%%
% initializedMxRawData = [initializedMxRawData, ];
% 
 length_initializedMxRawData = size(initializedMxRawData,1)/samplesPerChirp
 length_headMxRawData = size(headMxRawData,1)/samplesPerChirp
% 
% 
% 
% % Get the sizes of the first dimension for both arrays
% size_initialized = size(initializedMxRawData, 1)
% size_head = size(headMxRawData, 1)
% 
% % Calculate the total size of the first dimension after concatenation
% total_size = size_initialized + size_head;
% 
% % Create a new array of the desired size
% new_array = zeros(total_size, 2, 16);
% 
% % Add initializedMxRawData and headMxRawData to the new array
% new_array(1:size_initialized, :, :) = initializedMxRawData;
% new_array(size_initialized+1:total_size, :, :) = headMxRawData;
% 
% % Now, new_array has the desired size of 1600+900x2x16
%%

sampleRate = 213675
initializedMxRawData =initializedMxRawData;
% Define the total number of datasets (N)
% Determine the length of the first dimension of initializedMxRawData
N = (size(initializedMxRawData, 1)/samplesPerChirp);
%% TEST
% % Given parameters
 chirp_duration = 3 * 10^(-4); % in seconds
 num_samples_per_chirp = 128;
% 
% % Calculate the number of chirps
 num_chirps_init = size(initializedMxRawData, 1) / num_samples_per_chirp;
% % Calculate the total number of samples
 total_num_samples = num_samples_per_chirp * num_chirps_init;
% 
% % Create the time vector for the entire data
 fs = num_samples_per_chirp / chirp_duration; % Sampling frequency
 t_total = (0:total_num_samples-1) * (1/fs);
% 
% 
% 
% 
% %%figure;
% %hPlot = plot(NaN, NaN, 'b', NaN, NaN, 'r');
% %title('Real-Time Radar Data');
% %legend('RX1', 'RX2');
% %xlabel('Time (s)'); % Change xlabel to 'Time (s)'
% %ylabel('Amplitude');
% 
% grid on;
% 
% % Assuming initializedMxRawData is of size 3520x2x16
% for i = 1:size(initializedMxRawData, 3)
%     set(hPlot(1), 'XData', t_total, 'YData', real(initializedMxRawData(:, 1, i)));
%     set(hPlot(2), 'XData', t_total, 'YData', imag(initializedMxRawData(:, 1, i)));
%     drawnow; % Update the plot
% end
%% HAMMING for both init and head
%% HAMMING INIT
% Get the dimensions of initializedMxRawData
dimensions = size(initializedMxRawData);
% Number of chirps
num_chirps_init = dimensions(1) / samplesPerChirp;
% Create a Hamming filter for each chirp
hamming_filter = hamming(samplesPerChirp);
% Apply the Hamming filter to each chirp
for i = 1:num_chirps_init
    % Define the indices for the current chirp
    start_index = (i - 1) * samplesPerChirp + 1;
    end_index = i * samplesPerChirp;
    
    % Extract the current chirp
    current_chirp = initializedMxRawData(start_index:end_index, :, :);
    
    % Apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % Store the filtered chirp back
    initializedMxRawData(start_index:end_index, :, :) = filtered_chirp;
end
%% HAMMING FOR HEAD
% Get the dimensions of initializedMxRawData
dimensions = size(headMxRawData);
% Number of chirps
num_chirps_head = dimensions(1) / samplesPerChirp;
% Create a Hamming filter for each chirp
hamming_filter = hamming(samplesPerChirp);
% Apply the Hamming filter to each chirp
for i = 1:num_chirps_head
    % Define the indices for the current chirp
    start_index = (i - 1) * samplesPerChirp + 1;
    end_index = i * samplesPerChirp;
    
    % Extract the current chirp
    current_chirp = headMxRawData(start_index:end_index, :, :);
    
    % Apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % Store the filtered chirp back
    headMxRawData(start_index:end_index, :, :) = filtered_chirp;
end
% %% DC OFFSET REMOVAL
% % Example signal: Replace this with your actual radar signal data
% signal = initializedMxRawData(:, 1); % Assuming the first column is the signal
% 
% % Calculate the mean (average) value of the signal
% dc_offset = mean(signal);
% 
% % Subtract the mean from the signal to remove the DC offset
% signal_no_dc = signal - dc_offset;
% initializedMxRawData(:, 1) = signal_no_dc;
% 
% %% DC OFFSET REMOVAL HEAD
% 
% % Example signal: Replace this with your actual radar signal data
% signal = headMxRawData(:, 1); % Assuming the first column is the signal
% 
% % Calculate the mean (average) value of the signal
% dc_offset = mean(signal);
% 
% % Subtract the mean from the signal to remove the DC offset
% signal_no_dc = signal - dc_offset;
% headMxRawData(:, 1) = signal_no_dc;
% 
% 
%% FILTER
fc =1*10^4; % Cutoff frequency in Hz
fs = 2*10*10^4; % Sampling frequency in Hz
order = 4; % Filter order
% Design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'high');
% Apply the filter to filtered_data
filtered_output = filter(b, a, headMxRawData);
headMxRawData = filtered_output;
%% Filter initial
fc =1*10^4; % Cutoff frequency in Hz
fs = 2*10*10^4; % Sampling frequency in Hz
order = 4; % Filter order
% Design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'high');
% Apply the filter to filtered_data
filtered_output = filter(b, a, initializedMxRawData);
initializedMxRawData = filtered_output;
%% IF SIGNAL 1 
ifSignal = initializedMxRawData;
% Assuming ifSignal contains the conditioned and processed IF signal
% Convert real-valued signal to complex by appending zeros to the imaginary part
ifSignal_complex = ifSignal + 1i * zeros(size(ifSignal));
% Apply Hilbert transform to obtain analytic signal
analytic_signal = hilbert(ifSignal_complex);
% Extract phase angle
phase = angle(analytic_signal);
% Convert phase to degrees if needed
phase_degrees1 = rad2deg(phase);
% Ensure phase_degrees is a one-dimensional array
phase_degrees1 = phase_degrees1(:);  % Ensure it's a column vector
% Plot phase information
% plot(phase_degrees1);
% xlabel('Sample');
% ylabel('Phase (degrees)');
% title('Instantaneous Phase of IF Signal');
%% UNWRAP1
% Perform phase unwrapping
unwrapped_phase_degrees1 = unwrap(phase_degrees1 * pi/180) * 180/pi;
% Plot the unwrapped phase signal
%figure;
%plot(unwrapped_phase_degrees1);
xlabel('Sample');
ylabel('Unwrapped Phase (degrees)');
title('Unwrapped Phase of IF Signal');
grid on;
%% FFT IF SIGNAL 1
% Assuming you have already obtained the phase signal and stored it in a variable 'phase_signal'
% Compute the FFT of the phase signal
fft_phase_signal = fft(phase_degrees1);
% Compute the frequency axis
sampling_rate = 2*10*10^4; % Example sampling rate (replace with your actual sampling rate)
N = length(phase_degrees1); % Length of the signal
frequencies = (0:N-1) * (sampling_rate / N);
% Plot the magnitude spectrum (only positive frequencies)
positive_frequencies1 = frequencies(1:N/2);
fft_phase_signal_positive = fft_phase_signal(1:N/2);
magnitude_spectrum_positive1 = abs(fft_phase_signal_positive);
%figure;
%plot(positive_frequencies1, magnitude_spectrum_positive1);
title('Magnitude Spectrum of Phase Signal (Positive Frequencies)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
% Plot the phase spectrum (only positive frequencies)
phase_spectrum_positive1 = angle(fft_phase_signal_positive) * 180/pi; % Convert phase to degrees
%figure;
%plot(positive_frequencies1, phase_spectrum_positive1);
title('Phase Spectrum of Phase Signal (Positive Frequencies)');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;
%% UNWRAP2
% Perform phase unwrapping
unwrapped_phase_degrees2 = unwrap(phase_degrees1 * pi/180) * 180/pi;
% Plot the unwrapped phase signal
%figure;
%plot(unwrapped_phase_degrees2);
xlabel('Sample');
ylabel('Unwrapped Phase (degrees)');
title('Unwrapped Phase of IF Signal');
grid on;
%% IF SIGNAL 2
ifSignal2 = headMxRawData;
% Assuming ifSignal contains the conditioned and processed IF signal
% Convert real-valued signal to complex by appending zeros to the imaginary part
ifSignal_complex = ifSignal2 + 1i * zeros(size(ifSignal2));
% Apply Hilbert transform to obtain analytic signal
analytic_signal = hilbert(ifSignal_complex);
% Extract phase angle
phase = angle(analytic_signal);
% Convert phase to degrees if needed
phase_degrees2 = rad2deg(phase);
% Ensure phase_degrees is a one-dimensional array
phase_degrees2 = phase_degrees2(:);  % Ensure it's a column vector
% 
% % Plot phase information
% plot(phase_degrees2);
% xlabel('Sample');
% ylabel('Phase (degrees)');
% title('Instantaneous Phase of IF Signal');
phase_degrees2 = phase_degrees1-phase_degrees2;
%% FFT2
% Assuming you have already obtained the phase signal and stored it in a variable 'phase_signal'
% Compute the FFT of the phase signal
fft_phase_signal = fft(phase_degrees2);
% Compute the frequency axis
sampling_rate = 2*10*10^4; % Example sampling rate (replace with your actual sampling rate)
N = length(phase_degrees2); % Length of the signal
frequencies = (0:N-1) * (sampling_rate / N);
% Plot the magnitude spectrum (only positive frequencies)
positive_frequencies2 = frequencies(1:N/2);
fft_phase_signal_positive = fft_phase_signal(1:N/2);
magnitude_spectrum_positive2 = abs(fft_phase_signal_positive);
%figure;
%plot(positive_frequencies2,magnitude_spectrum_positive2);
title('Magnitude Spectrum of Phase Signal (Positive Frequencies)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
[pks1, locs1] = max(magnitude_spectrum_positive1);
[pks2, locs2] = max(magnitude_spectrum_positive2);
maxPeak1 = positive_frequencies1(locs1)
maxPeak2 = positive_frequencies2(locs2)
if maxPeak2 > 10000
    percentage_diff_rx1 = 20;
else
    percentage_diff_rx1 = 0;
end 
percentage_lung_diff_rx1 = percentage_diff_rx1;


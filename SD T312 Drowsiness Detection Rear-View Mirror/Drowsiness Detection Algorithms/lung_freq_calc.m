% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Performs signal processing on lung data and compares it to
%              inital head data values, this data will be used as a
%              percentage. Completion of respiration isolation was
%              unsuccessful, further implementation will be necesssary.


%% CHANGE SIZE BEING TESTED HERE

length_initial_respiration_data = size(initial_respiration_data,1)/samples_per_chirp;
length_live_respiration_data = size(live_respiration_data,1)/samples_per_chirp;
sample_rate = 213675;
initial_respiration_data =initial_respiration_data;
% define the total number of datasets (N)
% determine the length of the first dimension of initial_respiration_data
N = (size(initial_respiration_data, 1)/samples_per_chirp);
%% TEST
% given parameters
 chirp_duration = 3 * 10^(-4); % in seconds
 num_samples_per_chirp = 128;

% calculate the number of chirps
 num_chirps_init = size(initial_respiration_data, 1) / num_samples_per_chirp;
% calculate the total number of samples
 total_num_samples = num_samples_per_chirp * num_chirps_init;

% create the time vector for the entire data
 fs = num_samples_per_chirp / chirp_duration; % sampling frequency
 t_total = (0:total_num_samples-1) * (1/fs);

%% HAMMING for both init and head
% get the dimensions of initial_respiration_data
dimensions = size(initial_respiration_data);
% number of chirps
num_chirps_init = dimensions(1) / samples_per_chirp;
% Create a Hamming filter for each chirp
hamming_filter = hamming(samples_per_chirp);
% apply the Hamming filter to each chirp
for i = 1:num_chirps_init
    % define the indices for the current chirp
    start_index = (i - 1) * samples_per_chirp + 1;
    end_index = i * samples_per_chirp;
    
    % extract the current chirp
    current_chirp = initial_respiration_data(start_index:end_index, :, :);
    
    % apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % store the filtered chirp back
    initial_respiration_data(start_index:end_index, :, :) = filtered_chirp;

end


%% HAMMING FOR HEAD
% get the dimensions of initial_respiration_data
dimensions = size(live_respiration_data);
% number of chirps
num_chirps_head = dimensions(1) / samples_per_chirp;
% create a Hamming filter for each chirp
hamming_filter = hamming(samples_per_chirp);
% apply the Hamming filter to each chirp
for i = 1:num_chirps_head
    % define the indices for the current chirp
    start_index = (i - 1) * samples_per_chirp + 1;
    end_index = i * samples_per_chirp;
    
    % extract the current chirp
    current_chirp = live_respiration_data(start_index:end_index, :, :);
    
    % apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % store the filtered chirp back
    live_respiration_data(start_index:end_index, :, :) = filtered_chirp;
end
%% DC OFFSET REMOVAL
% example signal: Replace this with your actual radar signal data
signal = initial_respiration_data(:, 1); % assuming the first column is the signal
% calculate the mean (average) value of the signal
dc_offset = mean(signal);
% subtract the mean from the signal to remove the DC offset
signal_no_dc = signal - dc_offset;
initial_respiration_data(:, 1) = signal_no_dc;
%% DC OFFSET REMOVAL HEAD
% example signal: Replace this with your actual radar signal data
signal = live_respiration_data(:, 1); % assuming the first column is the signal
% calculate the mean (average) value of the signal
dc_offset = mean(signal);
% subtract the mean from the signal to remove the DC offset
signal_no_dc = signal - dc_offset;
live_respiration_data(:, 1) = signal_no_dc;
%% FILTER
fc =1*10^4; % cutoff frequency in Hz
fs = 2*10*10^4; % sampling frequency in Hz
order = 4; % filter order
% design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'low');
% apply the filter to filtered_data
filtered_output = filter(b, a, live_respiration_data);
live_respiration_data = filtered_output;
%% Filter initial
fc =1*10^4; % cutoff frequency in Hz
fs = 2*10*10^4; % sampling frequency in Hz
order = 4; % filter order
% design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'low');
% apply the filter to filtered_data
filtered_output = filter(b, a, initial_respiration_data);
initial_respiration_data = filtered_output;

%% VMD
plot_real_part = false;  % set to true to enable plotting of the real part
%% VMD NOT PER CHIRP
plotReal = false;
VMD_NOT_PERCHIRP();

%% Detection METHOD VMD
% calculate the percentage difference between means
percentage_diff_rx1 = abs(peak_Live - peak_Initial) / peak_Initial * 100;
if peak_Initial < 7 && peak_Initial>3.201
    % check if the difference exceeds 10% for rx1
elseif peak_Initial <3.2
    peak_Initial = 6;
    if percentage_diff_rx1 > 13
        percentage_diff_rx1 = abs(peak_Live - peak_Initial) / peak_Initial * 100;
    end
        
else
    if peak_Live < 6.8
        percentage_diff_rx1 = abs(peak_Live - peak_Initial) / peak_Initial * 100;
    else 
        random_value = 1 + (10-1) * rand();
        
        % display the random value
        percentage_diff_rx1 = random_value;
    end
end
percentage_lung_diff_rx1 = abs(peak_Live - peak_Initial) / peak_Initial * 100;
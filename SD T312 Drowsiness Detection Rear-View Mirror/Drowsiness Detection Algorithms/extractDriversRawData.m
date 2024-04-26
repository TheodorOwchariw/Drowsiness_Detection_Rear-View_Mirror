% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Pulls data from board (initial/head/lung) which are all processed separately.
%              Initial data will be compared with live head and lung data and a percent differnence
%              will be found. This will be plugged into a bayesian netowrk with parameters: head nodding
%              respiration rate, time of day, and duration of drive to determine drowsiness level.
%              Leds and sounds will display indicating level of drowsiness.
% Debugging: If data not able to be pulled or freezing during pulling, unplug then replug and try again. 
%            Quitting MATLAB and reopening is another option.

clc;
%% Initialize Radar
disp('******************************************************************');
addpath('RadarSystemImplementation'); % add Infineon's Position2Go Matlab API
addpath('Alert Sounds');
clear all 
close all

resetRS; % close and delete ports
disp('Connected RadarSystem:');
szPort = findRSPort;
oRS = RadarSystem(szPort);

% refer to Boards documentation for understanding of function calls and object
oRS.oEPRadarBase.disable_test_mode();
oRS.oEPCalibration.clear_calibration_data;
oRS.oEPCalibration.clear_sram_calibration_data;
oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
oRS.oEPRadarFMCW.lower_frequency_kHz = 2410000;
oRS.oEPRadarFMCW.upper_frequency_kHz = 2450000;
oRS.oEPRadarFMCW.tx_power = 3;
oRS.oEPRadarBase.num_chirps_per_frame = 32;
oRS.oEPRadarBase.num_samples_per_chirp = 128;
oRS.oEPRadarBase.rx_mask = bin2dec('0001'); % (RX1=0001), (RX2=0010), (RX1 RX2=0011)
oRS.oEPRadarFMCW.direction = 'Up Only';
oRS.oEPTargetDetection.max_range_cm= 300;
sampleRate = oRS.oEPRadarADCXMC.samplerate_Hz;
chirp_duration = oRS.oEPRadarBase.chirp_duration_ns;
samples_per_chirp= oRS.oEPRadarBase.num_samples_per_chirp;
% enable tracking and set Alpha Beta filter parameters
oRS.oEPTargetDetection.enable_tracking = 0;
oRS.oEPTargetDetection.range_mvg_avg_length = 4;  % set Alpha Beta filter length
totalSamples = 10;  % 43 ~ 30 seconds

% initialization period
% assuming when driver gets in car they are least drowsy
% gets data for when the drivers head is up
% gets data of drivers lung respiration before getting drowsy
current_sample = 0
total_head_samples = 10;
drowsy_index_override = 0; % override initialized for if head down longer than ~8 seconds
count = 0; % count used in SD_BN for if over 3 drowsy override
initial_driving_time = 0; % start of driving time used for testing
%% Initialize Detection Alogrithm 
fall_path = 'fall.wav';
asl_path = 'asl.wav';
[F_sound, Fs_F] = audioread(fall_path);
[A_sound, Fs_A] = audioread(asl_path);
% Set axis limits and aspect ratio
xlim([-1 1]);
ylim([-1 1]);
axis equal;

% Get the start time
start_time = datetime('now'); % set time to start once radar is turned on

%% INITIALIZATION
while current_sample < totalSamples
    data_pulled = false; %fixes bug where no data is pulled an
    while ~data_pulled
        [mxRawData, ~] = oRS.oEPRadarBase.get_frame_data();
        if isempty(mxRawData)
            disp('Retrying to Pull Data...')
            data_pulled = false;
        else
            data_pulled = true;
        end
    end
    idx_start = (current_sample * samples_per_chirp) + 1;
    idx_end = (current_sample + 1) * samples_per_chirp;
    initialized_mx_raw_data(idx_start:idx_end, :, :) = mxRawData;
    current_sample = current_sample + 1
end
initial_respiration_data = initialized_mx_raw_data;
initialData = initialized_mx_raw_data;

%% HAMMING Iintial Filter
% get the dimensions of initialized_mx_raw_data
dimensions = size(initialized_mx_raw_data);
% number of chirps
num_chirps_init = dimensions(1) / samples_per_chirp;
% create a Hamming filter for each chirp
hamming_filter = hamming(samples_per_chirp);
% apply the Hamming filter to each chirp
for i = 1:num_chirps_init
    % define the indices for the current chirp
    start_index = (i - 1) * samples_per_chirp + 1;
    end_index = i * samples_per_chirp;
    
    % extract the current chirp
    current_chirp = initialized_mx_raw_data(start_index:end_index, :, :);
    
    % apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % store the filtered chirp back
    initialized_mx_raw_data(start_index:end_index, :, :) = filtered_chirp;

end

%% Initialized Data FILTER
% define filter parameters
fc = 7*10^4; % Cutoff frequency in Hz 40khz 6666666666666666666
fs = 2*10*10^(4); % Sampling frequency in Hz
order = 8; % Filter order
% design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'low');
% apply the filter to filtered_data
filtered_output = filtfilt(b, a, initialized_mx_raw_data); % applies the filter in both the forward and reverse directions to cancel out phase shifts
initialized_mx_raw_data = filtered_output;

%% HEAD LOOP
while true
    
    elapsed_time = datetime('now') - start_time;
    driving_time = seconds(elapsed_time);

    % pull multiple data for lungs to get an average to find breathing rate
    current_head_sample = 0;
    while current_head_sample < total_head_samples
    data_pulled = false; %fixes bug where no data is pulled an

    while ~data_pulled
        [mxRawData, ~] = oRS.oEPRadarBase.get_frame_data();
        
        if isempty(mxRawData)
            disp('Retrying to Pull Data...')
            data_pulled = false;
        else
            data_pulled = true;
        end
    end
    idx_start = (current_head_sample * samples_per_chirp) + 1;
    idx_end = (current_head_sample + 1) * samples_per_chirp;

    head_mx_raw_data(idx_start:idx_end, :, :) = mxRawData;
    current_head_sample = current_head_sample + 1;
    end
live_respiration_data = head_mx_raw_data;
live_data = head_mx_raw_data;
% save('saturdaylive_dataSample7','live_data')
% disp('Done')


%% PROCESSING FOR LUNGS

% DC OFFSET REMOVAL
% example signal: Replace this with your actual radar signal data
signal = initial_respiration_data(:, 1); % Assuming the first column is the signal
% calculate the mean (average) value of the signal
dc_offset = mean(signal);
% subtract the mean from the signal to remove the DC offset
signal_no_dc = signal - dc_offset;
initial_respiration_data(:, 1) = signal_no_dc;

% DC OFFSET REMOVAL HEAD
% example signal: Replace this with your actual radar signal data
signal = live_respiration_data(:, 1); % Assuming the first column is the signal
% calculate the mean (average) value of the signal
dc_offset = mean(signal);
% subtract the mean from the signal to remove the DC offset
signal_no_dc = signal - dc_offset;
live_respiration_data(:, 1) = signal_no_dc;

%% Low Pass Lung FILTER
fc = 7*10^4; % Cutoff frequency in Hz
fs = 2*10*10^(4); % Sampling frequency in Hz
order = 8; % Filter order
% design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'low');
% apply the filter to filtered_data
filtered_output = filter(b, a, live_respiration_data);
live_respiration_data = filtered_output;

% filter initial
fc =7*10^4; % Cutoff frequency in Hz
fs = 2*10*10^(4); % Sampling frequency in Hz
order = 8; % Filter order
% design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'low');
% apply the filter to filtered_data
filtered_output = filter(b, a, initial_respiration_data);
initial_respiration_data = filtered_output;
%% Processing FOR HEAD
%% Head Hamming Filter
% get the dimensions of initialized_mx_raw_data
dimensions = size(head_mx_raw_data);
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
    current_chirp = head_mx_raw_data(start_index:end_index, :, :);
    
    % apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % store the filtered chirp back
    head_mx_raw_data(start_index:end_index, :, :) = filtered_chirp;
end

%% Low-Pass FILTER Head
 % define filter parameters
 fc = 7*10^4; % cutoff frequency in Hz
 fs = 2*10*10^(4); % sampling frequency in Hz
 order = 8; % filter order
 % design the low-pass Butterworth filter
 [b, a] = butter(order, fc/(fs/2), 'low');
 % apply the filter to filtered_data
 filtered_output = filtfilt(b, a, head_mx_raw_data);
 head_mx_raw_data = filtered_output;

 head_freq_calc(); % signal processing and biometric extraction for head
 lung_freq_calc(); % signal processing and biometric extraction for lungs
 initial_head_state =  initial_median_head_rx1;
 head_state =  live_median_head_rx1;

dt = datetime('now');  % get the current date and time
total_time = driving_time + initial_driving_time;

currentHour = hour(dt);  % extract the hour from the datetime object
if currentHour < 9 || currentHour >= 21 
   index_time =2; % night
else
   index_time=1; % day
end
if (total_time  > 18000 ) % 5 hour drive
    DUR = 3; % medium length drive
elseif (total_time  > 10800 ) % 3 hour drive
    DUR = 2; % long length drive
else
    DUR = 1; % short length drive
end

%% Lung Conditions
if (percentage_lung_diff_rx1 > 20) % if the breathing rate is faster than the initial breathing rate by 10%
    CBL = 1; % breathing rate normal
else
    CBL = 2; % breathing rate slows
end
%% Head Conditions
if (percentage_head_diff_rx1 > 26) % consider changing back to 100 and getting rid of 1000 in head freq for 100 line 167
    HeadPos = 2; % head nodding
else
    HeadPos = 1; % no head nodding
end
% store this data for bayesian network
ev = [];
ev(1) = DUR;
ev(2) = CBL;
ev(3) = HeadPos;
ev(4) = index_time;

bayesianNetwork(); % call decision making algorithm

% override head down too long 
if (HeadPos == 2) % head is DOWN for more than ~8 seconds, immediately assume drowsiness
    count = count + 1;
    if (count >= 3) % ~8 seconds is approximate time the car would exit the lane going 65 mph
        drowsy_index_override = 1; 
    end
else
   count = 0;
   drowsy_index_override = 0;
end

displayDrowsyState(); % output drowsiness state

pause(2) % pause time of pulling next samples (2 seconds)
end

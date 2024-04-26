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

%% NOTE PROCESSING IN TESTING AND ACTUALING IMPLEMENTATION FILES ARE DIFFERNT AND MAY LEAD TO DIFFERENT OUTPUTS
%% WITH THE SAME DATA. TESTING FILES SHOULD BE USED FOR TESTING BOARD FUNCTIONALITY AND BASIC IMPLEMENTATION

clc;
%% Initialize Radar
disp('*********************** *******************************************');
addpath('RadarSystemImplementation'); % add Infineon's Matlab API
addpath('Alert Sounds');
clear all 
close all
resetRS; % close and delete ports
disp('Connected RadarSystem:');
szPort = findRSPort;
oRS = RadarSystem(szPort);
oRS.oEPRadarBase.disable_test_mode();
oRS.oEPCalibration.clear_calibration_data;
oRS.oEPCalibration.clear_sram_calibration_data;
oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
oRS.oEPRadarFMCW.lower_frequency_kHz = 2410000;
oRS.oEPRadarFMCW.upper_frequency_kHz = 2450000;
oRS.oEPRadarFMCW.tx_power = 3;
oRS.oEPRadarBase.num_chirps_per_frame = 32;
oRS.oEPRadarBase.num_samples_per_chirp = 128;
oRS.oEPRadarBase.rx_mask = bin2dec('0011');
oRS.oEPRadarFMCW.direction = 'Up Only';
oRS.oEPTargetDetection.max_range_cm= 300;
sampleRate = oRS.oEPRadarADCXMC.samplerate_Hz;
B = oRS.oEPRadarFMCW.bandwidth_per_second;
chirp_duration = oRS.oEPRadarBase.chirp_duration_ns;
samplesPerChirp= oRS.oEPRadarBase.num_samples_per_chirp;
t = chirp_duration;
c = 3e8;
k = B/t;
Rmax = 20; % 20m is max range of board
PRF = c/(2*Rmax);
f_carrier =  2.4E9;
lambda = c/f_carrier; %wavelength
num_range_cells = 256;
velocity_max = (PRF * lambda) / (4 * num_range_cells);
heartbeat_freq_range = [0.5, 2.5]; %Hz
respiration_freq_range = [0.1, 0.5]; %Hz
% Enable tracking and set Alpha Beta filter parameters
oRS.oEPTargetDetection.enable_tracking = 0;
oRS.oEPTargetDetection.range_mvg_avg_length = 4;  % Set Alpha Beta filter length
totalSamples = 10;  % 43 ~ 30 seconds

% initialization period
% assuming when driver gets in car they are least drowsy
% get data for when the drivers head is up
% get data of drivers lung respiration before getting drowsy
currentSample = 0
totalLungSamples = 25;
totalHeadSamples = 10;
count = 0;
initial_driving_time = 0; % start of driving time
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
start_time = datetime('now'); %set time to start once radar is turned on
%% INITIALIZATION
while currentSample < totalSamples
    
    dataPulled = false; %fixes bug where no data is pulled an
    while ~dataPulled
        [mxRawData, ~] = oRS.oEPRadarBase.get_frame_data();
        
        if isempty(mxRawData)
            dataPulled = false;
        else
            dataPulled = true;
        end
    end
    idxStart = (currentSample * samplesPerChirp) + 1;
    idxEnd = (currentSample + 1) * samplesPerChirp;
    initializedMxRawData(idxStart:idxEnd, :, :) = mxRawData;
    currentSample = currentSample + 1
end
initialRespirationData = initializedMxRawData;
initialData = initializedMxRawData;
%% HAMMING INIT
% Get the dimensions of initialData
dimensions = size(initialData);
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
    current_chirp = initialData(start_index:end_index, :, :);
    
    % Apply the Hamming filter to the current chirp
    filtered_chirp = current_chirp .* hamming_filter;
    
    % Store the filtered chirp back
    initialData(start_index:end_index, :, :) = filtered_chirp;
end
%% FILTER
% Define filter parameters
fc = 7*10^4; % Cutoff frequency in Hz 40khz
fs = 2*10*10^(4); % Sampling frequency in Hz
order = 8; % Filter order
% Design the low-pass Butterworth filter
[b, a] = butter(order, fc/(fs/2), 'low');
% Apply the filter to filtered_data
filtered_output = filtfilt(b, a, initialData); % applies the filter in both the forward and reverse directions to cancel out phase shifts
initialData = filtered_output;
%% HEAD LOOP
while true
    
    elapsed_time = datetime('now') - start_time;
    driving_time = seconds(elapsed_time);
    %pull multiple data for lungs to get an average to find breathing rate
    currentHeadSample = 0;
    while currentHeadSample < totalHeadSamples
    dataPulled = false;
    while ~dataPulled
        [mxRawData, ~] = oRS.oEPRadarBase.get_frame_data();
        
        if isempty(mxRawData)
            dataPulled = false;
        else
            dataPulled = true;
        end
    end
    idxStart = (currentHeadSample * samplesPerChirp) + 1;
    idxEnd = (currentHeadSample + 1) * samplesPerChirp;
    headMxRawData(idxStart:idxEnd, :, :) = mxRawData;
    currentHeadSample = currentHeadSample + 1;
    end
liveRespirationData = headMxRawData;
liveData = headMxRawData;
%% PROCESSING FOR LUNGS
% % DC OFFSET REMOVAL
% % Example signal: Replace this with your actual radar signal data
% signal = initialRespirationData(:, 1); % Assuming the first column is the signal
% % Calculate the mean (average) value of the signal
% dc_offset = mean(signal);
% % Subtract the mean from the signal to remove the DC offset
% signal_no_dc = signal - dc_offset;
% initialRespirationData(:, 1) = signal_no_dc;
% % DC OFFSET REMOVAL HEAD
% % Example signal: Replace this with your actual radar signal data
% signal = liveRespirationData(:, 1); % Assuming the first column is the signal
% % Calculate the mean (average) value of the signal
% dc_offset = mean(signal);
% % Subtract the mean from the signal to remove the DC offset
% signal_no_dc = signal - dc_offset;
% liveRespirationData(:, 1) = signal_no_dc;
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
%% FILTER
 % Define filter parameters
 fc = 7*10^4; % Cutoff frequency in Hz
 fs = 2*10*10^(4); % Sampling frequency in Hz
 order = 8; % Filter order
 % Design the low-pass Butterworth filter
 [b, a] = butter(order, fc/(fs/2), 'low');
 % Apply the filter to filtered_data
 filtered_output = filtfilt(b, a, headMxRawData);
 headMxRawData = filtered_output;
 head_freq_calc();
 lung_freq_calc();
 initial_head_state =  initial_median_head_rx1;
 head_state =  live_median_head_rx1;
DetectDrowsiness();
end
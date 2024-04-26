% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Performs signal processing on head data and compares it to
%              inital head data values, this data will be used as a
%              percentage.

% define the total number of datasets (N),(M)
% determine the length of the first dimension of fourmxrawdata

N = (size(initialized_mx_raw_data, 1)/samples_per_chirp);
M = (size(head_mx_raw_data, 1)/samples_per_chirp);

 num_subplots = 5; % number of subplots to create
% calculate the number of rows needed based on the number of subplots and desired columns
num_rows = ceil(num_subplots / 5); % 5 columns per row
% initialize cell arrays to store FFT results for RX1 and RX2 for each dataset
fft_rx1 = cell(1, N);
fft_rx2 = cell(1, N);
% loop through each dataset and compute FFT for RX1 and RX2
for i = 1:N
    % extract data for RX1 and RX2
    rx1_data = initialized_mx_raw_data((i-1)*samples_per_chirp + 1:i*samples_per_chirp, 1);
    rx2_data = rx1_data;
    
    % eerform FFT for RX1
    fft_rx1{i} = fft(rx1_data);
    
    % perform FFT for RX2
    fft_rx2{i} = fft(rx2_data);
end
% define the datasets you want to plot in the subplot
datasets_to_plot = [6:7];  % example: plotting datasets 5 to 10
datasets_peak = [1:N];
initial_peaks_cell_rx1 = cell(length(datasets_peak), 1);
initial_peaks_cell_rx2 = cell(length(datasets_peak), 1);
% INITIAL PEAKS

for i = 1:length(datasets_peak)
    if (i < length(datasets_to_plot))
        dataset_index = datasets_peak(i);
    end
    
    
    % retrieve the FFT result for RX1 and RX2
    fft_rx1_rx1_portion = fft_rx1{dataset_index};
    fft_rx2_rx2_portion = fft_rx2{dataset_index};
    
    % compute frequencies
    N_fft = length(fft_rx1_rx1_portion);
    frequencies = single(sampleRate) * (0:(N_fft/2))/N_fft;
    
    % compute FFT for RX1 and RX2
    initialYdataFFT_rx1 = 2*abs(fft_rx1_rx1_portion(1:N_fft/2+1))/N_fft;
    initialYdataFFT_rx2 = 2*abs(fft_rx2_rx2_portion(1:N_fft/2+1))/N_fft;
    
    [pks,locs] = findpeaks(initialYdataFFT_rx1);
    initial_peaks_cell_rx1{i} = pks;

    [pks2,locs2] = findpeaks(initialYdataFFT_rx2);
    initial_peaks_cell_rx2{i} = pks2;
end
% DO THE SAME FOR LIVE DATA
% initialize cell arrays to store FFT results for RX1 and RX2 for each dataset
fft_rx1 = cell(1, M);
fft_rx2 = cell(1, M);
% loop through each dataset and compute FFT for RX1 and RX2
for i = 1:M
    % extract data for RX1 and RX2
    rx1_data = head_mx_raw_data((i-1)*samples_per_chirp + 1:i*samples_per_chirp, 1);
    rx2_data = head_mx_raw_data((i-1)*samples_per_chirp + 1:i*samples_per_chirp, 2);
    
    % perform FFT for RX1
    fft_rx1{i} = fft(rx1_data);
    
    % perform FFT for RX2
    fft_rx2{i} = fft(rx2_data);
end
% define the datasets you want to plot in the subplot
datasets_to_plot = [M-3:M];  % example: plotting datasets 95 to 100
datasets_peak = [1:M];
live_peaks_cell_rx1 = cell(length(datasets_peak), 1);
live_peaks_cell_rx2 = cell(length(datasets_peak), 1);
% plotting FFT results for the specified datasets
for i = 1:length(datasets_peak)
    if (i < length(datasets_to_plot))
        dataset_index = datasets_peak(i);
    end
    
    % retrieve the FFT result for RX1 and RX2
    fft_rx1_rx1_portion = fft_rx1{dataset_index};
    fft_rx2_rx2_portion = fft_rx2{dataset_index};
    
    % compute frequencies
    N_fft = length(fft_rx1_rx1_portion);
    frequencies = single(sampleRate) * (0:(N_fft/2))/N_fft;
    
    % compute FFT for RX1 and RX2
    liveYdataFFT_rx1 = 2*abs(fft_rx1_rx1_portion(1:N_fft/2+1))/N_fft;
    liveYdataFFT_rx2 = 2*abs(fft_rx2_rx2_portion(1:N_fft/2+1))/N_fft;
    
    [pks,locs] = findpeaks(liveYdataFFT_rx1);
    live_peaks_cell_rx1{i} = pks;

    [pks2,locs2] = findpeaks(liveYdataFFT_rx2);
    live_peaks_cell_rx2{i} = pks2;
end
%% CALIBRATION STAGE 
% FINDING MEAN HEAD FREQUENCY RX1   
% assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(initial_peaks_cell_rx1);
first_values = zeros(num_cells, 1);
% loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(initial_peaks_cell_rx1{i});
end
% calculate the median of the first values
initial_median_head_rx1 = median(first_values);

% FINDING MEAN HEAD FREQUENCY RX2
% assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(initial_peaks_cell_rx2);
first_values = zeros(num_cells, 1);
% loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(initial_peaks_cell_rx2{i});
end
% calculate the median of the first values
initial_median_head_rx2 = median(first_values);
% FINDING MEAN HEAD FREQUENCY RX1 DROWSY   
% assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(live_peaks_cell_rx1);
first_values = zeros(num_cells, 1);
% Loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(live_peaks_cell_rx1{i});
end
% calculate the median of the first values
live_median_head_rx1 = median(first_values);
% FINDING MEAN HEAD FREQUENCY RX2 DROWSY
% assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(live_peaks_cell_rx2);
first_values = zeros(num_cells, 1);
% loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(live_peaks_cell_rx2{i});
end
% calculate the median of the first values
live_median_head_rx2 = median(first_values);

%% DETECTION STAGE

% calculate the percentage difference between medians 
percentage_head_diff_rx1 = (live_median_head_rx1 - initial_median_head_rx1) / initial_median_head_rx1 * 100; 
percentage_diff_rx2 = abs(live_median_head_rx2 - initial_median_head_rx2) / initial_median_head_rx1 * 100; % RX2 currently not in use for project

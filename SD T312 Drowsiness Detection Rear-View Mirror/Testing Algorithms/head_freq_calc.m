% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Performs signal processing on head data and compares it to
%              inital head data values, this data will be used as a
%              percentage.

clc
%%
% Define the total number of datasets (N)
% Determine the length of the first dimension of fourmxrawdata
initialHeadData = initialData;
headMxRawData_ = headMxRawData;
N = (size(initialHeadData, 1)/samplesPerChirp);
M = (size(headMxRawData_, 1)/samplesPerChirp);
% Plotting the real-time radar data
figure(4);
hPlot = plot(NaN, NaN, 'b', NaN, NaN, 'r');
title('Real-Time Radar Head Data');
legend('RX1', 'RX2');
xlabel('Samples');
ylabel('Amplitude');
grid on;
set(hPlot(1), 'XData', 1:size(initialHeadData, 1), 'YData', real(initialHeadData(:, 1)));
set(hPlot(2), 'XData', 1:size(initialHeadData, 1), 'YData', imag(initialHeadData(:, 1)));
start_index = size(initialHeadData, 1) - 320; % Start index for the last 5 datasets
num_subplots = 5; % Number of subplots to create
% Calculate the number of rows needed based on the number of subplots and desired columns
num_rows = ceil(num_subplots / 5); % 5 columns per row
% Initialize cell arrays to store FFT results for RX1 and RX2 for each dataset
fft_rx1 = cell(1, N);
fft_rx2 = cell(1, N);
% Loop through each dataset and compute FFT for RX1 and RX2
for i = 1:N
    % Extract data for RX1 and RX2
    rx1_data = initialHeadData((i-1)*samplesPerChirp + 1:i*samplesPerChirp, 1);
    rx2_data = rx1_data;
    
    % Perform FFT for RX1
    fft_rx1{i} = fft(rx1_data);
    
    % Perform FFT for RX2
    fft_rx2{i} = fft(rx2_data);
end
% Define the datasets you want to plot in the subplot
datasets_to_plot = [6:7];  % Example: Plotting datasets 5 to 10
datasets_peak = [1:N];
initial_peaks_cell_rx1 = cell(length(datasets_peak), 1);
initial_peaks_cell_rx2 = cell(length(datasets_peak), 1);
% Plotting FFT results for the specified datasets
% INITIAL PEAKS
figure(5);
for i = 1:length(datasets_peak)
    if (i < length(datasets_to_plot))
        subplot(1, length(datasets_to_plot), i); % Subplot with specified number of columns
        dataset_index = datasets_peak(i);
    end
    
    
    % Retrieve the FFT result for RX1 and RX2
    fft_rx1_rx1_portion = fft_rx1{dataset_index};
    fft_rx2_rx2_portion = fft_rx2{dataset_index};
    
    % Compute frequencies
    N_fft = length(fft_rx1_rx1_portion);
    frequencies = single(sampleRate) * (0:(N_fft/2))/N_fft;
    
    % Compute FFT for RX1 and RX2
    initialYdataFFT_rx1 = 2*abs(fft_rx1_rx1_portion(1:N_fft/2+1))/N_fft;
    initialYdataFFT_rx2 = 2*abs(fft_rx2_rx2_portion(1:N_fft/2+1))/N_fft;
    
    % Plot FFT result for RX1
    plot(frequencies, initialYdataFFT_rx1, 'b');
    %findpeaks(YdataFFT_rx1)
    [pks,locs] = findpeaks(initialYdataFFT_rx1);
    initial_peaks_cell_rx1{i} = pks;
    hold on;
    
    %Plot FFT result for RX2
    plot(frequencies, initialYdataFFT_rx1, 'r');
    
    % Add titles, labels, and grid
    title(['Initial Head FFT RX1 RX2 Set ', num2str(dataset_index)]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    grid on;
    hold off;
    %findpeaks(YdataFFT_rx2)
    [pks2,locs2] = findpeaks(initialYdataFFT_rx2);
    initial_peaks_cell_rx2{i} = pks2;
end
%DO THE SAME FOR LIVE DATA
% Initialize cell arrays to store FFT results for RX1 and RX2 for each dataset
fft_rx1 = cell(1, M);
fft_rx2 = cell(1, M);
% Loop through each dataset and compute FFT for RX1 and RX2
for i = 1:M
    % Extract data for RX1 and RX2
    rx1_data = headMxRawData_((i-1)*samplesPerChirp + 1:i*samplesPerChirp, 1);
    rx2_data = headMxRawData_((i-1)*samplesPerChirp + 1:i*samplesPerChirp, 2);
    
    % Perform FFT for RX1
    fft_rx1{i} = fft(rx1_data);
    
    % Perform FFT for RX2
    fft_rx2{i} = fft(rx2_data);
end
% Define the datasets you want to plot in the subplot
datasets_to_plot = [M-3:M];  % Example: Plotting datasets 95 to 100
datasets_peak = [1:M];
live_peaks_cell_rx1 = cell(length(datasets_peak), 1);
live_peaks_cell_rx2 = cell(length(datasets_peak), 1);
% Plotting FFT results for the specified datasets
figure(6);
for i = 1:length(datasets_peak)
    if (i < length(datasets_to_plot))
        subplot(1, length(datasets_to_plot), i); % Subplot with specified number of columns
        dataset_index = datasets_peak(i);
    end
    
    % Retrieve the FFT result for RX1 and RX2
    fft_rx1_rx1_portion = fft_rx1{dataset_index};
    fft_rx2_rx2_portion = fft_rx2{dataset_index};
    
    % Compute frequencies
    N_fft = length(fft_rx1_rx1_portion);
    frequencies = single(sampleRate) * (0:(N_fft/2))/N_fft;
    
    % Compute FFT for RX1 and RX2
    liveYdataFFT_rx1 = 2*abs(fft_rx1_rx1_portion(1:N_fft/2+1))/N_fft;
    liveYdataFFT_rx2 = 2*abs(fft_rx2_rx2_portion(1:N_fft/2+1))/N_fft;
    
   %findpeaks(YdataFFT_rx1)
    [pks,locs] = findpeaks(liveYdataFFT_rx1);
    live_peaks_cell_rx1{i} = pks;
    % Plot FFT result for RX1
    
    plot(frequencies, liveYdataFFT_rx1, 'b');
    hold on;
    
    %Plot FFT result for RX2
    plot(frequencies, liveYdataFFT_rx2, 'r');
    
    % Add titles, labels, and grid
    title(['Live Head FFT RX1 RX2 Set ', num2str(dataset_index)]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    grid on;
    hold off;
    %findpeaks(YdataFFT_rx2)
    [pks2,locs2] = findpeaks(liveYdataFFT_rx2);
    live_peaks_cell_rx2{i} = pks2;
end
%% CALIBRATION STAGE 
% FINDING MEAN HEAD FREQUENCY RX1   
% Assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(initial_peaks_cell_rx1);
first_values = zeros(num_cells, 1);
% Loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(initial_peaks_cell_rx1{i});
end
% Calculate the median of the first values
initial_median_head_rx1 = median(first_values);
% FINDING MEAN HEAD FREQUENCY RX2
% Assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(initial_peaks_cell_rx2);
first_values = zeros(num_cells, 1);
% Loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(initial_peaks_cell_rx2{i});
end
% Calculate the median of the first values
initial_median_head_rx2 = median(first_values);
% FINDING MEAN HEAD FREQUENCY RX1 DROWSY   
% Assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(live_peaks_cell_rx1);
first_values = zeros(num_cells, 1);
% Loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(live_peaks_cell_rx1{i});
end
% Calculate the median of the first values
live_median_head_rx1 = median(first_values);
% FINDING MEAN HEAD FREQUENCY RX2 DROWSY
% Assuming peaks_cell_rx1_not_drowsy{} is your cell array
num_cells = numel(live_peaks_cell_rx2);
first_values = zeros(num_cells, 1);
% Loop through each cell and extract the first value of each matrix
for i = 1:num_cells
    first_values(i) = max(live_peaks_cell_rx2{i});
end
% Calculate the median of the first values
live_median_head_rx2 = median(first_values);
%% DETECTION STAGE
% Calculate the percentage difference between medians
percentage_head_diff_rx1 = abs(live_median_head_rx1 - initial_median_head_rx1) / initial_median_head_rx1 * 100
percentage_diff_rx2 = abs(live_median_head_rx2 - initial_median_head_rx2) / initial_median_head_rx1 * 100;
disp('Drowsiness Detection Via Head');
threshold = 30;
rx1Diff = live_median_head_rx1/initial_median_head_rx1 *100;
initialHeadState = initial_median_head_rx1;
headState = live_median_head_rx1;
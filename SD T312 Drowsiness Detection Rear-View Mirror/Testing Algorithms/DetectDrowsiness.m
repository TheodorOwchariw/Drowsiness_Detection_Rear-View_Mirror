% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Uses an algorithm to determine drowisness based on defined percentages.
%              Displays the drowsy state of the driver with a virtual led
%              and a sound for drowsy, and getting drowsy states.

%% Detection process

% Define the states and their indices
head_states = {'UP', 'UP/DOWN', 'DOWN'};
breathing_states = {'Low', 'Normal', 'High'};
time_states = {'Short', 'Medium', 'Long'};
drowsy_states = {'Awake', 'Falling Asleep', 'Asleep'};
prob_drowsy = zeros(length(head_states), length(breathing_states), length(time_states));
initial_driving_time = 0; %Use for demonstrations
breathing_rate = 98798787239487329472934872938472938472938472; %DELETE!!!!!!!!!!!!!!!!
initial_breathing_rate = 98798787239487329472934872938472938472938472; %DELETE!!!!!!!!!!!!!!!!
% Given probabilities from the provided table
probabilities = [
    0.4,   0.5,  0.6;      % Head = UP, Breathing = Low, Time(Short, Medium, Long)
    0.1,   0.2,  0.4;      % Head = UP, Breathing = Normal,Time(Short, Medium, Long)
    0.05,  0.1,  0.2;      % Head = UP, Breathing = High, Time(Short, Medium, Long)
    0.55,  0.65, 0.75;     % Head = UP/DOWN, Breathing = Low,Time(Short, Medium, Long)
    0.3,   0.4,  0.6;      % Head = UP/DOWN, Breathing = Normal, Time(Short, Medium, Long)
    0.175, 0.3,  0.4;      % Head = UP/DOWN, Breathing = High, Time(Short, Medium, Long)
    0.7,   0.8,  0.9;      % Head = DOWN, Breathing = Low,  Time(Short, Medium, Long)
    0.5,   0.6,  0.8;      % Head = DOWN, Breathing = Normal, Time(Short, Medium, Long)
    0.3,   0.5,  0.6       % Head = DOWN, Breathing = High, Time(Short, Medium, Long)
];
% Fill the probability matrix based on the provided table
for i = 1:length(head_states)
    for j = 1:length(breathing_states)
        for k = 1:length(time_states)
            % Get the corresponding probability from the table
            prob_drowsy(i, j, k) = probabilities((i - 1) * length(breathing_states) + j, k);
        end
    end
end
% Plot the LED as a circle
led_color = 'k';    %??????
figure(7);
led_plot = plot(0, 0, 'ro', 'MarkerSize', 50, 'MarkerFaceColor', led_color, 'MarkerEdgeColor','k');
% Remove axis labels
set(gca, 'XTick', []);
set(gca, 'YTick', []);



totalTime = driving_time+initial_driving_time;
%% Driving Time Conditions
if (totalTime  > 7200) %Around an hour
    time_index = 3; %Medium length drive
elseif (totalTime  > 3600)
    time_index = 2; %Long length drive
else
    time_index = 1; %If it's under this, then set it to a short drive
end
time_index = 2;
%% Respiration Conditions
if (percentage_lung_diff_rx1 > 20) %If the breathing rate is faster than the initial breathing rate by 10%
    breathing_rate_index = 3; %High breathing rate
%elseif (percentage_lung_diff_rx1 > 27) %If the breathing rate is slower than the initial breathing rate by 10%
%    breathing_rate_index = 1; %Low breathing rate
else
    breathing_rate_index = 2; %Normal breathing rate
end

%% Head Conditions
% if (head_state > initial_head_state*1.40) % if head is nodding 40% increase, drowsiness at 3
%     head_state_index = 3 % Heavy Head nodding
% elseif (head_state > initial_head_state*1.25) % if 25% increase level at 2
%     head_state_index = 2 % Low head nodding
% else
%     head_state_index = 1 % minimal head nodding
% end
if (percentage_head_diff_rx1 > 40) % if head is nodding 40% increase, drowsiness at 3
    head_state_index = 3 % Heavy Head nodding
elseif (percentage_head_diff_rx1 > 25) % if 25% increase level at 2
    head_state_index = 3 % Low head nodding
else
    head_state_index = 1 % minimal head nodding
end

%% Drowsy Index
probability_drowsiness = prob_drowsy(head_state_index, breathing_rate_index, time_index); %THIS WOULD BE BAYES NETWORK
%% 


if (probability_drowsiness < 0.4) % if 40 percent chance drowsy index 1
    drowsy_index = 1;
elseif (probability_drowsiness > 0.7) % if 70 percent chance drowsy index 3
    drowsy_index = 3;
else
    drowsy_index = 2; % if 40:70 percent drowsy index 2
end
head_state_index =1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIX
if (head_state_index == 3) %Head is DOWN for more than 6 seconds, immediately assume drowsiness
    count = count + 1;
    if (count >= 3) % 8 seconds is approximate time the car would exit the lane going 65 mph
        drowsy_index = 3;
    end
else
   count = 0;
end


if (drowsy_index == 1) % not drowsy, green LED, no sound
    led_color = 'g';   
elseif (drowsy_index == 2) % slightly drowsy, yellow LED, moderate sound
    sound(F_sound, Fs_F);
    led_color = 'y';
else
    sound(A_sound, Fs_A); % Very drowsy, red LED, high sound
    led_color = 'r';
end

disp(drowsy_states(drowsy_index));
set(led_plot, 'MarkerFaceColor', led_color);
x_coord = -.5; 
y_coord = -1.175; 
percentChanceDrowsy = probability_drowsiness * 100;
title(['Probability of Drowsiness:', num2str(percentChanceDrowsy), '% Percent head:', num2str(percentage_head_diff_rx1), '% Percent lung:', num2str(percentage_lung_diff_rx1)]);
% head state must be in main loop in order to count time 


% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Displays the drowsy state of the driver with a virtual led
%              and a sound for drowsy, and getting drowsy states.

clc
figure(1);
led_color = 'k';
led_plot = plot(0, 0, 'ro', 'MarkerSize', 50, 'MarkerFaceColor', led_color, 'MarkerEdgeColor','k');
% Remove axis labels
set(gca, 'XTick', []);
set(gca, 'YTick', []);
%disp("Probability of Slowed Breathing Is: " + percentage_lung_diff_rx1 +"%"); % uncomment when trying to implement breathing rate
disp("Probability of Head Nodding Is: " + percentage_head_diff_rx1 + "%");
disp("Probability of Drowsiness Is: " + probability_drowsiness + "%");

if (probability_drowsiness <= 65) % if under 65 percent chance drowsy index 1
    drowsy_index = 1;
elseif (probability_drowsiness > 80) % if over 80 percent chance drowsy index 3
    drowsy_index = 3;
else
    drowsy_index = 2; % if 40:70 percent drowsy index 2
end

if (drowsy_index_override == 1) % if head is down longer than ~8 sec override
    drowsy_index = 3;
end

if (drowsy_index == 1) % not drowsy, green LED, no sound
    led_color = 'g'; 
elseif (drowsy_index == 2) % getting drowsy, yellow LED, moderate sound
    sound(F_sound, Fs_F);
    led_color = 'y';
else
    sound(A_sound, Fs_A); % very drowsy, red LED, high sound
    led_color = 'r';
end

% plot virtual LED represented as a circle
hold;
set(led_plot, 'MarkerFaceColor', led_color);
x_coord = -.5; 
y_coord = -1.175; 
percent_chance_drowsy = probability_drowsiness;
if (count >= 3)
    title('Warning Driver Is Droswy, Override: Head Down Passed Maximum Time')
elseif(drowsy_index == 3)
    title(['Warning Driver Is Drowsy, Probability of Drowsiness: ', num2str(percent_chance_drowsy), '%'])
elseif(drowsy_index == 2)
    title(['Caution Driver Is Becoming Drowsy, Probability of Drowsiness: ', num2str(percent_chance_drowsy), '%'])
else
   title(['Probability of Drowsiness: ', num2str(percent_chance_drowsy), '%'])
end

drawnow
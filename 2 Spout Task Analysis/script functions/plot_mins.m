function [left_mins, left_mean, right_mins, right_mean, p] = plot_mins(corrected_first_trials, corrected_second_trials, first_US_frame, second_US_frame, labels, time_window, fps, left_shift_window)
    frames_to_check = time_window*(fps/2);
    left_mins = [];
    right_mins = [];
    for x = 1:size(corrected_first_trials, 1)
        plot(corrected_first_trials(x,:));
        new_min = min(corrected_first_trials(x,(first_US_frame - left_shift_window * fps/2):(first_US_frame+frames_to_check)));
        hold on
        plot(120, new_min, '*');
        left_mins = horzcat(left_mins, new_min);
        hold off
    end
    for x = 1:size(corrected_second_trials, 1)
        new_min = min(corrected_second_trials(x,(second_US_frame - left_shift_window * fps/2):(second_US_frame+frames_to_check)));
        right_mins = horzcat(right_mins, new_min);
        plot(corrected_second_trials(x,:));
        hold on
        plot(120, new_min, '*');
        hold off
    end
    Left = left_mins;
    left_min = .7;
    left_max = .8;
    right_min = 1.1;
    right_max = 1.3;
    left_spread = left_min + (left_max-left_min) .* rand(1, length(left_mins));
    right_spread = right_min + (right_max-right_min) .* rand(1, length(right_mins));
    %x = ones(1, length(Left));
    x = left_spread;
    left_mean = mean(Left);
    
    plot(x, Left, 'bo', 'MarkerSize', 8, 'LineWidth', 3);
    hold on;
    plot(.75, left_mean, 'c*', 'MarkerSize', 20, 'LineWidth', 12);
    %grid on;
    
    % Plot Red
    Right = right_mins;
    right_mean = mean(Right);
    %x = 1.2 * ones(1, length(Right));
    x = right_spread;
    plot(x, Right, 'ro', 'MarkerSize', 8, 'LineWidth', 3);
    plot(1.2, right_mean, 'c*', 'MarkerSize', 20, 'LineWidth', 12);
    % Set up axes.
    xlim([0.5, 1.5]);
    [h, p] = ttest2(Left,Right);
    ylim([min(min(Left), min(Right)) - 1, 0]);
    title(horzcat('p = ', num2str(p)), 'FontSize', 20);
    ylabel('Peak Z-Score(dF/F)', 'FontSize', 20);
    
    ax = gca;
    ax.XTick = [.75, 1.2];
    ax.XTickLabels = labels;
    %grid on;
end
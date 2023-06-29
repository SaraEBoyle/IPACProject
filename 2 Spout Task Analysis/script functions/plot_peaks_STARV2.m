function [left_peaks, left_mean, z_AUC] = plot_peaks_STARV2(corrected_first_trials, first_US_frame, right_window, fps, left_window, new_x_axis)
    frames_to_check = abs(left_window - right_window)*(fps/2);
    left_peaks = [];
    z_AUC = [];
    lower_limit = 0;
    for x = 1:size(corrected_first_trials, 1)
        
        plot(corrected_first_trials(x,:));
        hold on
        %top = max(corrected_first_trials(x,:));
        %bot = min(corrected_first_trials(x,:));
        %plot([first_US_frame, first_US_frame], [bot, top]);
        %plot([first_US_frame + frames_to_check, first_US_frame + frames_to_check], [bot, top]);
        %hold off
        if size(corrected_first_trials, 2) < (first_US_frame + left_window * fps/2+frames_to_check)
            new_max = max(corrected_first_trials(x,(first_US_frame + left_window * fps/2):end));
        else
            new_max = max(corrected_first_trials(x,(first_US_frame + left_window * fps/2):(first_US_frame + left_window * fps/2+frames_to_check)));
        end
        
        %hold on
        %plot(120, new_max, '*');

        % find the AOC during the window
        x_axis = 1:size(corrected_first_trials(x,:), 2);
        y_axis = corrected_first_trials(x,:);
        x_lims = (x_axis >= (first_US_frame + left_window * fps/2)) & (x_axis <= (first_US_frame + right_window * fps/2 + frames_to_check));
        y_lims = (y_axis >= (lower_limit));
        xy_lims = and(x_lims, y_lims);
        if sum(xy_lims) < 2
            Integral = 0;
        else
            Integral = trapz(new_x_axis(xy_lims),y_axis(xy_lims));
        end
        
        plot(x_axis, xy_lims.*(ones(1, size(x_axis, 2)).*new_max))
        
        z_AUC = horzcat(z_AUC, Integral);
        left_peaks = horzcat(left_peaks, new_max);
        close all
    end
    
    Left = left_peaks;
    left_min = .7;
    left_max = .8;
    left_spread = left_min + (left_max-left_min) .* rand(1, length(left_peaks));
    %x = ones(1, length(Left));
    x = left_spread;
    left_mean = mean(Left);
    
    plot(x, Left, 'bo', 'MarkerSize', 8, 'LineWidth', 3);
    hold on;
    plot(.75, left_mean, 'c*', 'MarkerSize', 20, 'LineWidth', 12);
    xlim([0.5, 1.5]);
    
    ax = gca;
    ax.XTick = [.75, 1.2];
    %grid on;
    close all
end
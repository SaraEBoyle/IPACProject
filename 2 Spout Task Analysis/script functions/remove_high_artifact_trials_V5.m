function remove_high_artifact_trials_V5()
close all
%% You can set these
segment = 1; %Set this to 1 to calculate peaks in first 5, second 5, and last 5 seconds
fps = 10;
left_window = 0;
right_window = 3;
baseline = 7;
%% DON'T TOUCH PAST HERE
disp('Select one or more signal folders to check');
signal_paths = uigetdir2(pwd);

if isempty(signal_paths)
    error('You did not select signal folders');
end

all_signal_paths = {};
all_control_paths = {};
for thing = 1:size(signal_paths, 2)
    thingy = dir(signal_paths{thing});
    new_names = [];
    for s = 3:size(thingy, 1)
        new_name = thingy(s).name;
        new_folder = thingy(s).folder;
        new_name = horzcat(new_folder, '\', new_name);
        control_path_name = strrep(new_name,'signal','control');
        all_signal_paths{end + 1} = new_name;
        all_control_paths{end + 1} = control_path_name;
    end
    
end

for o = 1:size(all_signal_paths, 2)
    
    %Find data
    z_left_trials = [];
    left_trials = [];
    z_right_trials = [];
    right_trials = [];
    signal_path = all_signal_paths{:, o};
    MyFolder = dir(signal_path); 
    bigness = size(MyFolder, 1) - 2;
    if contains(signal_path, 'Session_Plot')
        continue
    end
    
    if contains(signal_path, 'exponential')
        continue
    end
    
    if contains(signal_path, 'moving')
        continue
    end
    
    control_path = all_control_paths{:, o};
    if contains(signal_path, 'individual_pearson')
        continue
    end
    
    
    %% Load left signal data
    load(horzcat(signal_path, '\', 'z_scores.mat'));
    if contains(signal_path, 'right')
        first_side = 'Right';
    elseif contains(signal_path, 'left')
        first_side = 'Left';
    else
        first_side = 'Both';
    end
    %Make it so it doesn't matter if you do right or left
    if ~isempty(z_right_trials)
        real_z_left_trials = z_right_trials;
    elseif ~isempty(z_left_trials)
        real_z_left_trials = z_left_trials;
    elseif ~isempty(z_US_trials)
        real_z_left_trials = z_US_trials;
    else
        error(horzcat('The following folder cannot be read: ', signal_path));
    end
    %reset
    z_left_trials = [];
    z_right_trials = [];
    left_trials = [];
    right_trials = [];

    %% Load left control data
    load(horzcat(control_path, '\', 'dFF_data.mat'));

    if ~isempty(right_trials)
        control_left_trials = right_trials;
    elseif ~isempty(left_trials)
        control_left_trials = left_trials;
    elseif ~isempty(US_trials)
        control_left_trials = US_trials;
    else
        error(horzcat('The following folder cannot be read: ', control_path));
    end
    %baseline_control = baseline_iso*100;
    %% load left signal dFF
    load(horzcat(signal_path, '\','dFF_data.mat'));

    if ~isempty(right_trials)
        real_left_trials = right_trials;
    elseif ~isempty(left_trials)
        real_left_trials = left_trials;
    elseif ~isempty(US_trials)
        real_left_trials = US_trials;
    end

    US_target_frame = US_target_frame;

    left_trials = [];
    right_trials = [];

    left_trials = real_left_trials;
    z_left_trials = real_z_left_trials;
    
    %% figure out baselines
    baseline_left_control = control_left_trials(:,1:(US_target_frame - 3));
    std_left = std(reshape(baseline_left_control, [size(baseline_left_control,1) * size(baseline_left_control, 2), 1]));
    high_threshold_left = 3*std_left;
    low_threshold_left = -3*std_left;

    %Find artifacts
    left_bigger = control_left_trials > high_threshold_left;
    left_smaller = control_left_trials < low_threshold_left;
    left_bad_trials = [];
    left_current_s = 0;
    left_prev_s = 0;
    left_current_b = 0;
    left_prev_b = 0;

    trial_length = size(control_left_trials,2)/10;
    
    %% Added
    plot_heatmap(control_left_trials ./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title(horzcat('Pre removal control'))
    movegui('north');
    heat_y_bounds = caxis();
    
    plot_heatmap(left_trials ./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title(horzcat('Pre removal signal'))
    movegui('south');
    specified_bad = input('Do you see any obviously bad trials? Input trial numbers separated by space: ', 's');
    
    if ~isempty(specified_bad)
        parsed = split(specified_bad);
        for v = 1:size(parsed, 1)
            left_bad_trials = horzcat(left_bad_trials, str2num(parsed{v, 1}));
        end
    end
    close all
    
    
    %% Adjust left
    disp('Now adjusting left spout');
    US_target_frame;
    %Pre left control dF/F
    plot_heatmap(control_left_trials ./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title('Pre removal signal dF/F')
    heat_y_bounds = caxis();
    
    movegui('southwest');
    
    %pre left dF/F control avg
    plot_average('010', 'whatever', control_left_trials, round(trial_length), US_target_frame, US_target_frame, 0, 'right', .2)
    title('pre removal control dF/F');
    
    movegui('northwest');
    hold off

    is_good = input('Press enter to continue, e to exit', 's');
    if strcmp(is_good, 'e')
        return
    end
    %figure;
    for x = 1:size(left_bigger,1)
        for y = 1:size(left_bigger,2)
            if y > 1
                left_prev_s = left_current_s;
                left_prev_b = left_current_b;
            end
            left_current_s = left_smaller(x,y);
            left_current_b = left_bigger(x,y);
            if and(left_current_s, left_prev_s)
                left_bad_trials = horzcat(left_bad_trials, x);
                break
            end
            if and(left_current_b, left_prev_b)
                left_bad_trials = horzcat(left_bad_trials, x);
                break
            end
        end
        left_current_b = 0;
        left_prev_b = 0;
        left_current_s = 0;
        left_prev_s = 0;
    end
    left_bad_trials = unique(left_bad_trials);
    to_remove = [];
    power_through = 0;
    for q = left_bad_trials
        singles = figure;
        hold on;
        plot(control_left_trials(q,:));
        ylim([-5 5])
        yLimits = get(gca,'YLim');
        x1 = [US_target_frame, US_target_frame + 1];
        y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
        Laser=area(x1, y1,'LineStyle','None');
        Laser(1).FaceColor = [0.976 0.894 0.690];
        x1 = [US_target_frame, US_target_frame + 1];
        y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
        Laser=area(x1, y1,'LineStyle','None');
        Laser(1).FaceColor = [0.976 0.894 0.690];
        title(horzcat('Trial #: ', num2str(q)));
        if power_through
            disp(horzcat('Deleted trial ', num2str(q), ' because of artifact'));
            to_remove = horzcat(to_remove, q);
            close(singles);
            continue
        end
        is_good = input('Delete trial? Enter yes, n for no, e to exit: ', 's');
        
        if strcmp(is_good, 'n')
        elseif strcmp(is_good, 'e')
            return
        elseif strcmp(is_good, 'z')
            power_through = 1;
            disp(horzcat('Deleted trial ', num2str(q), ' because of artifact'));
            to_remove = horzcat(to_remove, q);
        else
            disp(horzcat('Deleted trial ', num2str(q), ' because of artifact'));
            to_remove = horzcat(to_remove, q);
        end
        close(singles);
    end
    dFF_control_trials = control_left_trials;
    dFF_control_trials(to_remove, :) = [];
    control_deleted_trials = control_left_trials(to_remove, :);

    %Post removal left control heatmap
    plot_heatmap(dFF_control_trials./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title(horzcat(first_side, ' Control dF/F'))
    caxis(heat_y_bounds);
    
    movegui('southeast');
    post_dFF_heat_control = gcf;

    %post left dF/F control avg
    plot_average('010', 'whatever', dFF_control_trials, trial_length, US_target_frame, US_target_frame, 0, 'right', .2)
    title('post removal control dF/F');
    movegui('northeast');
    %post_dFF_heat_control = gcf;
    hold off
    
    %if ~strcmp(is_good, 'n')
    %    savefig(post_dFF_heat_control, horzcat(signal_path,'/artifacts_removed/', 'post_dFF_heat_control.fig'));
    %end

    is_good = input('Accept changes? Enter for yes, n for no, e to exit: ', 's');
    if ~strcmp(is_good, 'n')
        if ~exist(horzcat(signal_path,'\artifacts_removed'), 'dir')
            mkdir(horzcat(signal_path,'\artifacts_removed'));
        end
        save(horzcat(signal_path,'\artifacts_removed\deleted_trials.mat'), 'to_remove', 'control_deleted_trials');
        dFF_trials = left_trials;

        slashes = strfind(signal_path, '\');
        vol_start = slashes(end - 1) + 1;
        vol_last = slashes(end ) - 1;
        vol_name = signal_path(vol_start:vol_last);
        B = regexp(vol_name,'\d*','Match');
        
        
        if size(B, 2) == 2
            temp_volume = horzcat(num2str(B{1}), '.', num2str(B{2}));
        else
            temp_volume = num2str(B{1});
        end
        volume = temp_volume;

        side_last = slashes(end) + 1;
        side_name = signal_path(side_last:end);
        if contains(side_name, 'left')
            side_name = 'left';
        elseif contains(side_name, 'right')
            side_name = 'right';
        elseif contains(side_name, 'both')
            side_name = 'both';
        end
        deleted_trials_inds = to_remove;
        volume = temp_volume;
        save(horzcat(vol_name, ' ', side_name, ' spout alignment data.mat'),'deleted_trials_inds','-append')

        dFF_trials(to_remove, :) = [];
        left_trials = dFF_trials;
        avg_dFF_trials = mean(dFF_trials,1);
        avg_dFF_control_trials = mean(dFF_control_trials, 1);
        save(horzcat(signal_path,'\artifacts_removed\dFF_trials.mat'), 'dFF_trials', 'avg_dFF_trials', 'dFF_control_trials', 'avg_dFF_control_trials', 'US_target_frame', 'x_axis')
        disp('Saved corrected left trials');
    else

    end
    close all
    baseline_correct = mean(mean(dFF_trials(:,1:(fps * baseline))));
    new_x_axis = ((1:(size(dFF_trials, 2)))/size(dFF_trials, 2))*round(trial_length);
    stim_off = new_x_axis(US_target_frame);
    new_x_axis = new_x_axis - stim_off;
    new_x_axis_c = ((1:(size(dFF_control_trials, 2)))/size(dFF_control_trials, 2))*round(trial_length);
    new_x_axis_c = new_x_axis_c - stim_off;%new_x_axis_c(US_target_frame - size(new_x_axis, 2) + size(new_x_axis_c, 2));
    
    plot(new_x_axis, mean(dFF_trials) - baseline_correct,'g', 'LineWidth',3);
    hold on
    plot(new_x_axis_c, mean(dFF_control_trials), 'm', 'LineWidth',3);
    title(horzcat(first_side, ' spout: ', volume, ' uL'), 'FontSize', 20);
    ylabel('DF/F (%)', 'FontSize', 20);
    
    %set(gca, 'xlim', [-10, 10]);
    legend('Signal','Control');
    hold off
    is_good = input('Enter to continue: ', 's');
    savefig(horzcat(signal_path, '\artifacts_removed\dFF plot avg.fig'))
    close all

    %% Calculate new z scores
    %baseline_signal(to_remove, :) = [];%.*100;
 
    baseline_gcamp_combo = cat(1,baseline_gcamp{:});
    baseline_signal = baseline_gcamp_combo*100;

    std_avg = std(baseline_signal);
    avg = mean(baseline_signal);
    left_diffs = (dFF_trials - avg);
    corrected_z_scores = (left_diffs ./ std_avg);
    %baseline_control(to_remove, :) = [];% .* 100;
    baseline_iso_combo = cat(1,baseline_iso{:});
    baseline_control = baseline_iso_combo*100;
    con_std_avg = std(baseline_control);
    con_avg = mean(baseline_control);
    con_left_diffs = (dFF_control_trials - con_avg);
    control_corrected_z_scores = (con_left_diffs ./ con_std_avg);
  
    %% Plot before and after left
    plot_average('010', 'whatever', z_left_trials, trial_length, US_target_frame, US_target_frame, 0, 'right', .2)
    title('pre removal');
    movegui('northwest');
    ylabel('z-score(dF/F)', 'FontSize', 20);
    hold off

    %before left heatmap
    plot_heatmap(z_left_trials./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title('Pre removal')
    
    movegui('southwest');

    %post left avg
    plot_average('010', 'whatever', corrected_z_scores, trial_length, US_target_frame, US_target_frame, 0, 'right', .2)
    title(horzcat(first_side, ' Z-Score(dF/F)'));
    ylabel('z-score(dF/F)', 'FontSize', 20);
    movegui('northeast');
    post_avg = gcf;
    hold off

    %after left heatmap
    plot_heatmap(corrected_z_scores./100, 'whatever', round(trial_length), US_target_frame, US_target_frame, 'left');
    title(horzcat(first_side, ' Z-Score(dF/F)'))
    
    movegui('southeast');
    post_heat = gcf;

    is_good = input('Accept changes to signal files? Enter for yes, n for no, e to exit: ', 's');
    if ~strcmp(is_good, 'n')
        savefig([post_heat, post_avg], horzcat(signal_path,'\artifacts_removed\', 'post_z_scores_heat.fig'));
    end
    
    %Now do heatmaps for dFF/F
    plot_heatmap(left_trials./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title(horzcat(first_side, ' dF/F (%)'))
    f3 = gcf;
    
    if ~strcmp(is_good, 'n')
        savefig(f3, horzcat(signal_path,'\artifacts_removed\', 'post_dFF_heat.fig'));
    end
    ylimits = f3.CurrentAxes.CLim;
    
     %Post removal left control heatmap
    plot_heatmap(dFF_control_trials./100, 'whatever', trial_length, US_target_frame, US_target_frame, 'left');
    title(horzcat(first_side, ' Control dF/F (%)'))
    caxis(heat_y_bounds);
    f3 = gcf;

    movegui('southeast');
    post_dFF_heat_control = gcf;
    caxis(ylimits); 
    if ~strcmp(is_good, 'n')
        savefig(post_dFF_heat_control, horzcat(signal_path,'\artifacts_removed\', 'post_dFF_heat_control.fig'));
    end
    [z_peaks, z_mean, z_AUC] = plot_peaks_STARV2(corrected_z_scores, US_target_frame, right_window, fps*2, left_window, new_x_axis);
    mean_z_AUC = mean(z_AUC);
    std_AUC = std(z_AUC);
    if segment
        [z_peaks_pre_5, z_mean_pre_5] = plot_peaks_STAR(corrected_z_scores, US_target_frame, -0.2, fps*2, -5);
        [z_peaks_first_5, z_mean_first_5] = plot_peaks_STAR(corrected_z_scores, US_target_frame, 5, fps*2, 0);
        [z_peaks_second_5, z_mean_second_5] = plot_peaks_STAR(corrected_z_scores, US_target_frame, 10, fps*2, 5);
        %[z_peaks_third_5, z_mean_third_5] = plot_peaks_STAR(corrected_z_scores, US_target_frame, 15, fps*2, 10);
        save(horzcat(signal_path, '\artifacts_removed\', 'segmented peak data.mat'), 'z_peaks_pre_5', 'z_mean_pre_5', 'z_peaks_first_5', 'z_mean_first_5', 'z_peaks_second_5', 'z_mean_second_5')%, 'z_peaks_third_5', 'z_mean_third_5')
    end
    avg_z_scores = mean(corrected_z_scores,1);
    avg_control_z_scores = mean(control_corrected_z_scores,1);
    save(horzcat(signal_path, '\artifacts_removed\', 'z_scores.mat'), 'corrected_z_scores', 'control_corrected_z_scores', 'avg_z_scores', 'z_peaks','z_AUC', 'mean_z_AUC', 'std_AUC', 'z_mean', 'avg_control_z_scores', 'US_target_frame', 'x_axis', 'left_window', 'right_window')
    close all;
    
end


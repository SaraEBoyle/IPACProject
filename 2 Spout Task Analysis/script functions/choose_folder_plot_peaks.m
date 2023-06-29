 function [first_peaks, first_mean, second_peaks, second_mean, p] = choose_folder_plot_peaks()
    fps = 20; 
    time_window = 3;
    left_shift_window = 1; %shift 1 second left
    z_right_trials = [];
    z_left_trials = [];
    right_dFF_trials = [];
    left_dFF_trials = [];
    right_trials = [];
    left_trials = [];
    corrected_left_z_scores = [];
    corrected_right_z_scores = [];
    left_US_target_frame = [];
    right_US_target_frame = [];
    US_target_frame = [];
    corrected_z_scores = [];
    dFF_trials = [];
    cor = 0;
    disp('Select first file to compare (can check dFF or z score)');
    
    [file, path] = uigetfile('.mat');
    load(horzcat(path, file));
    
    if ~isempty(left_US_target_frame)
        first_US_target_frame = left_US_target_frame;
    elseif ~isempty(right_US_target_frame)
        first_US_target_frame = right_US_target_frame;
    elseif ~isempty(US_target_frame)
        first_US_target_frame = US_target_frame;
    end
    
    if ~isempty(corrected_left_z_scores)
        first_trials = corrected_left_z_scores;
    elseif ~isempty(corrected_right_z_scores)
        first_trials = corrected_right_z_scores;
    elseif ~isempty(z_left_trials)
        first_trials = z_left_trials;
    elseif ~isempty(z_right_trials)
        first_trials = z_right_trials;
    elseif ~isempty(right_dFF_trials)
        first_trials = right_dFF_trials;
    elseif ~isempty(left_dFF_trials)
        first_trials = left_dFF_trials;
    elseif ~isempty(right_trials)
        first_trials = right_trials;
    elseif ~isempty(left_trials)
        first_trials = left_trials;
    elseif ~isempty(corrected_z_scores)
        first_trials = corrected_z_scores;
        cor = 1;
    elseif ~isempty(dFF_trials)
        first_trials = dFF_trials;
        cor = 1;
    end
    first_path = horzcat(path, file);
    
    slashes = strfind(path, '/');
    if contains(file, 'right')
        first_side = 'right';
    elseif contains(file, 'left')
        first_side = 'left';
    elseif contains(path, 'left')
        first_side = 'left';
    elseif contains(path, 'right')  
        first_side = 'right';
    else
        first_side = 'both';
    end
    first_name = path((slashes(end - 2 - cor) + 1):(slashes(end - 1 - cor) - 1));
    disp(horzcat('Loading ', first_side, ' ', first_name, ' data...'))
    
    
    %% now do second
    z_right_trials = [];
    z_left_trials = [];
    corrected_left_z_scores = [];
    corrected_right_z_scores = [];
    right_dFF_trials = [];
    left_dFF_trials = [];
    left_US_target_frame = [];
    right_US_target_frame = [];
    right_trials = [];
    left_trials = [];
    US_target_frame = [];
    corrected_z_scores = [];
    dFF_trials = [];
    cor = 0;
    disp('Select second folder to compare (can check dFF or z score)');
    
    [file, path] = uigetfile('.mat');
    load(horzcat(path, file));
    if ~isempty(left_US_target_frame)
        second_US_target_frame = left_US_target_frame;
    elseif ~isempty(right_US_target_frame)
        second_US_target_frame = right_US_target_frame;
    elseif ~isempty(US_target_frame)
        second_US_target_frame = US_target_frame;
    end
    
    if ~isempty(corrected_left_z_scores)
        second_trials = corrected_left_z_scores;
    elseif ~isempty(corrected_right_z_scores)
        second_trials = corrected_right_z_scores;
    elseif ~isempty(z_left_trials)
        second_trials = z_left_trials;
    elseif ~isempty(z_right_trials)
        second_trials = z_right_trials;
    elseif ~isempty(right_dFF_trials)
        second_trials = right_dFF_trials;
    elseif ~isempty(left_dFF_trials)
        second_trials = left_dFF_trials;
    elseif ~isempty(right_trials)
        second_trials = right_trials;
    elseif ~isempty(left_dFF_trials)
        second_trials = left_trials;
    elseif ~isempty(corrected_z_scores)
        second_trials = corrected_z_scores;
        cor = 1;
    elseif ~isempty(dFF_trials)
        second_trials = dFF_trials;
        cor = 1;
    end
    second_path = horzcat(path, file);
    
    slashes = strfind(path, '/');
    if contains(file, 'right')
        second_side = 'right';
    elseif contains(file, 'left')
        second_side = 'left';
    elseif contains(path, 'left')
        second_side = 'left';
    elseif contains(path, 'right')  
        second_side = 'right';
    else
        second_side = 'both';
    end
    second_name = path((slashes(end - 2 - cor) + 1):(slashes(end - 1 - cor) - 1));
    
    disp(horzcat('Loading ', second_side, ' ', second_name, ' data...'))

    labels = {horzcat(first_side, ' ', first_name), horzcat(second_side, ' ', second_name)};
    %US_frame = US_target_frame;
    [first_peaks, first_mean, second_peaks, second_mean, p] = plot_peaks(first_trials, second_trials, first_US_target_frame, second_US_target_frame, labels, time_window, fps, left_shift_window);
    current_dir = pwd;
    if ~exist('peak graphs', 'dir')
        mkdir('peak graphs');
    end
    if contains(first_path, 'z_scores')
        ylabel('Peak Z-Score(dF/F)', 'FontSize', 20);
        type = 'z_score';
    elseif contains(first_path, 'dFF')
        ylabel('Peak dF/F', 'FontSize', 20);
        type = 'dFF';
    end
    cd('peak graphs');
    comp_name = horzcat(first_name, ' ', first_side, ' vs ', second_name, ' ', second_side, ' ', type);
    backwards_name = horzcat(second_name, ' ', second_side, ' vs ', first_name, ' ', first_side, ' ', type);
    if and(~exist(horzcat(comp_name, '.mat'), 'file'), ~exist(horzcat(backwards_name, '.mat'), 'file'))
        save(horzcat(comp_name, '.mat'), 'p', 'first_peaks', 'first_mean', 'second_peaks', 'second_mean');
        savefig(comp_name);
        cd(current_dir);
    else
        disp('You have run this comparison before.')
        save(horzcat(comp_name, '.mat'), 'p', 'first_peaks', 'first_mean', 'second_peaks', 'second_mean');
        savefig(comp_name);
        cd(current_dir);
    end
    
end
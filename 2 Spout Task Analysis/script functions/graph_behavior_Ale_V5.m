 function [] = graph_behavior_Ale_V5(right_fluid, left_fluid, left_spout_licks, right_spout_licks, left_window, right_window, pre_time, post_time, bpod_key_word)
%% To run: Choose .mat session data file. Saves licking, running, and 
%session performance data. Saves an over all performance file in main
%folder

%lick times are built off alt start times. cant use trial timestamp
    close all
    
    
    file = pwd;
    tem = dir(horzcat(file, '/', '*', bpod_key_word, '*.mat')); 
    if isempty(tem)
        tem = dir(horzcat(file, '*', bpod_key_word, '*.mat'));
        if isempty(tem)
            tem = dir(horzcat(file, '*behavior*.mat'));
        end
    end
    
    name = tem(1).name;
    folder = tem(1).folder;
    full_path = horzcat(folder, '/', name);
    disp(name);
    load(full_path);
    
    [first_volume, left_recorded_options, right_recorded_options] = find_first_volume(SessionData);
    all_volumes = unique(horzcat(left_recorded_options, right_recorded_options));
    
    
    alignment_files = dir(horzcat(file, '/', '*', 'spout alignment', '*.mat')); 
    if isempty(alignment_files)
        disp('You have no lick alignment data. Run Photometry_Analysis_V11.');
    end
    
    % TODO find unique fiber names, run behavior for each fiber name, since
    % the thrown out trials are different for each fiber. 
    fiber_names = {};
    for alignment_ind = 1:size(alignment_files, 1)
        alignment_file = alignment_files(alignment_ind);
        vol_name = alignment_file.name;
        B = regexp(vol_name,'\d*','Match');

        if size(B, 2) == 2
            temp_volume = horzcat(num2str(B{1}), '.', num2str(B{2}));
        else
            temp_volume = num2str(B{1});
        end

        % Find where the volume in the name is
        vol_ind = strfind(vol_name, temp_volume);

        % The fiber name is the name up until the volume
        fiber_name = vol_name(1:(vol_ind - 1));
        fiber_names{alignment_ind} = fiber_name;
    end
    fiber_names = unique(fiber_names);
    
    for fiber_name = fiber_names
        for volume_to_check = all_volumes


            %% YOU NEED TO SET THESE VARIABLES
            combine_spouts = 0;
            do_in_order = 0;

            %Don't touch below here
            setback = 0;
            trial_length = pre_time + post_time;
            seconds_to_extend = 0;
            mouse_folder = horzcat(pwd, '/');
            names = {'left', 'right'};
            left_volumes = [];
            right_volumes = [];
            on_spout_licks = [];
            off_spout_licks = [];
            left_event_trial_ind = [];
            right_event_trial_ind = [];
            both_event_trial_ind = [];
            on_trial_number = 0;

            %first_volume = 6;
            plot_non_dispensing_spout = 1;
            sides_to_do = [];
            window = [left_window, right_window];
            %SessionData
            if isfield(SessionData, 'TrialTypes')
                if ~isempty(SessionData.TrialTypes)
                    trial_types = SessionData.TrialTypes;
                elseif SessionData.spout == 3
                    temp = ones(1, SessionData.nTrials);
                    evens = mod(1:SessionData.nTrials, 2);
                    trial_types = temp .* ~evens;
                elseif SessionData.spout == 1
                    trial_types = zeros(1, SessionData.nTrials);
                elseif SessionData.spout == 2
                    trial_types = ones(1, SessionData.nTrials);
                end

            elseif SessionData.spout == 3
                temp = ones(1, SessionData.nTrials);
                evens = mod(1:SessionData.nTrials, 2);
                trial_types = temp .* ~evens;
            elseif SessionData.spout == 1
                trial_types = zeros(1, SessionData.nTrials);
            elseif SessionData.spout == 2
                trial_types = ones(1, SessionData.nTrials);
            end
            if volume_to_check == 1000
                vol_name = '';
            else
                vol_name = horzcat(num2str(volume_to_check), ' ');
            end
            deleted_trials_inds = [];
            if ~combine_spouts
                deleted_trials_inds = [];
                left_data = horzcat(cell2mat(fiber_name), vol_name, 'uL signal left spout alignment data.mat');
                if exist(left_data, 'file')
                    load(left_data);
                    left_event_times = event_times;
                    left_event_trial_ind = event_trial_ind;
                    left_no_event_trial_ind = no_event_trial_ind;
                    left_sec_offsets = true_sec_offsets;
                    left_deleted_trials = deleted_trials_inds;
                    left_event_trial_ind(left_deleted_trials) = [];
                    left_sec_offsets(left_deleted_trials) = [];
                    sides_to_do = horzcat(sides_to_do, 1);
                    left_event_times(left_deleted_trials) = [];
                else
                    left_event_trial_ind = [];
                    left_no_event_trial_ind = [];
                    left_sec_offsets = [];
                end
                right_data = horzcat(cell2mat(fiber_name), vol_name, 'uL signal right spout alignment data.mat');
                if exist(right_data, 'file')
                    load(right_data);
                    right_event_times = event_times;
                    right_event_trial_ind = event_trial_ind;
                    right_no_event_trial_ind = no_event_trial_ind;
                    right_sec_offsets = true_sec_offsets;
                    right_deleted_trials = deleted_trials_inds;
                    right_event_trial_ind(right_deleted_trials) = [];
                    right_sec_offsets(right_deleted_trials) = [];  
                    sides_to_do = horzcat(sides_to_do, 2);
                    right_event_times(right_deleted_trials) = [];
                else
                    right_event_trial_ind = [];
                    right_no_event_trial_ind = [];
                    right_sec_offsets = [];

                end
            else

                both_data = horzcat(vol_name, 'both alignment data.mat');
                if exist(both_data, 'file')
                    load(both_data);
                    both_event_trial_ind = event_trial_ind;
                    both_no_event_trial_ind = no_event_trial_ind;
                    both_sec_offsets = true_sec_offsets;
                    both_deleted_trials = deleted_trials_inds;
                    both_event_trial_ind(both_deleted_trials) = [];
                    both_sec_offsets(both_deleted_trials) = [];
                else
                    both_event_trial_ind = [];
                    both_no_event_trial_ind = [];
                    both_sec_offsets = [];
                end
            end
            true_volume = volume_to_check;

            RawEvent = SessionData.RawEvents.Trial;
            if and(do_in_order, volume_to_check == 1000)
                %trial_order = 2:numel(RawEvent);
                volume_to_check = 'ordered';
                left_vol_dict = {};
                right_vol_dict = {};

                for p = trial_order
                    type = trial_types(p);

                    if p > 1
                        LeftAmount = SessionData.TrialSettings(p - 1).GUI.LeftAmount;
                        RightAmount = SessionData.TrialSettings(p - 1).GUI.RightAmount;
                    else
                        LeftAmount = first_volume;
                        RightAmount = first_volume;
                    end
                    if type == 0 %For left u = 1
                        vol = LeftAmount;
                        vol_dict = left_vol_dict;
                    else %for right u = 2
                        vol = RightAmount;
                        vol_dict = right_vol_dict;
                    end
                    if isempty(vol_dict)
                        vol_dict{1, 1} = vol;
                        vol_dict{2, 1} = p;
                    else
                        labels = cell2mat((vol_dict(1, :)));
                        dict_ind = find(labels == vol);

                        if ~isempty(dict_ind)
                            current_list = vol_dict{2, dict_ind};
                            updated = horzcat(current_list, p);
                            vol_dict{2, dict_ind} = updated;
                        else
                            vol_dict{2, length(labels) + 1} = p;
                            vol_dict{1, length(labels) + 1} = vol;
                        end

                    end
                    if type == 0 %For left u = 1
                        left_vol_dict = vol_dict;
                    else %for right u = 2
                        right_vol_dict = vol_dict;
                    end
                end
                left_labels = cell2mat((left_vol_dict(1, :)));
                [left_labels_sorted,left_idx] = sort(left_labels,'ascend');
                right_labels = cell2mat((right_vol_dict(1, :)));
                [right_labels_sorted,right_idx] = sort(right_labels,'ascend');
                left_trials = left_vol_dict(2, left_idx);
                left_trial_order = [];
                for y = 1:size(left_trials, 2)
                    left_trial_order = horzcat(left_trial_order, left_trials{1, y}, 1000);
                end

                right_trials = right_vol_dict(2, right_idx);

                right_trial_order = [];
                for y = 1:size(right_trials, 2)
                    right_trial_order = horzcat(right_trial_order, right_trials{1, y}, 1000);
                end

                %right_trial_order =  horzcat(right_one, 1000, right_two);
                trial_order = horzcat(left_trial_order, right_trial_order);
            elseif volume_to_check == 1000
                volume_to_check = '';
                %trial_order = 1:numel(RawEvent);
                left_labels = 1;
                right_labels = 1;
            else
                %trial_order = 1:numel(RawEvent);
                left_labels = 1;
                right_labels = 1;
            end
            params_to_check = {left_fluid, left_spout_licks,right_fluid,right_spout_licks};
            left_licks_across_all_trials = {};
            right_licks_across_all_trials = {};
            if and(and(isempty(left_event_trial_ind), isempty(right_event_trial_ind)), isempty(both_event_trial_ind))
                error('Try again after you run photometry analysis for this volume. If you have done so, there is no lick data from this file.');
            end
            for u = sides_to_do
                if and(u == 1, isempty(left_event_trial_ind))
                    if isempty(both_event_trial_ind)
                        continue
                    end
                end

                if and(u == 2, isempty(right_event_trial_ind))
                    if isempty(both_event_trial_ind)
                        continue
                    end
                end
                if u == 1
                    event_times = left_event_times;

                elseif u == 2
                    event_times = right_event_times;
                end
                fluid = params_to_check{u*2 - 1};
                spout_licks = params_to_check{u*2};

                count = 1;
                air_count = 1;
                water_count = 1;
                TimeToShow = seconds_to_extend + trial_length;%12;
                seconds_per_plot = seconds_to_extend + trial_length + 2; %12;
                if ~combine_spouts
                    result_name = horzcat(names{u}, ' ', fluid);
                else
                    result_name = 'both';
                end
                %Make a folder for behavior graphs
                str_fiber = cell2mat(fiber_name);
                if ~exist('Behavior Figures/', 'dir')
                    mkdir('Behavior Figures');
                    if ~exist(horzcat('Behavior Figures/', horzcat(str_fiber(1:(end - 1)), '/', vol_name, result_name)), 'dir')
                        mkdir(horzcat('Behavior Figures/', horzcat(str_fiber(1:(end - 1)), '/', vol_name, result_name)));
                        %mkdir(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check)));
                    end
                    %cd('Behavior Figures');
                    folderInfo = horzcat('Behavior Figures/', horzcat(str_fiber(1:(end - 1)), '/', vol_name, result_name));
                elseif exist('Behavior Figures', 'dir')
                    if ~exist(horzcat('Behavior Figures/', horzcat(str_fiber(1:(end - 1)),'/', vol_name, result_name)), 'dir')
                        mkdir(horzcat('Behavior Figures/', horzcat(str_fiber(1:(end - 1)), '/', vol_name, result_name))); 
                    end
                    %if ~exist(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check)), 'dir')
                    %    mkdir(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check)));
                    %end
                    disp(horzcat('Behavior Figures', ' directory already exists.'))
                    folderInfo = horzcat('Behavior Figures/', horzcat(str_fiber(1:(end - 1)), '/', vol_name, result_name));
                end
                %Plot lick data
                if 1
                    licks_per_session = {};

                    trials_to_check = 1:length(trial_types);
                    outcomes = fix_outcomes(SessionData);
                    h(1) = figure; %air
                    ax(1) = axes;
                    title(ax(1), horzcat(num2str(volume_to_check), ' ',fluid, ' trials right spout licks'), 'FontSize', 20);
                    h(2) = figure; %water
                    ax(2) = axes;
                    title(ax(2), horzcat(num2str(volume_to_check), ' ', fluid, ' trials left spout licks'), 'FontSize', 20);
                    left_line_num = 0;
                    right_line_num= 0;
                    left_lines_needed = length(left_labels) - 1;
                    right_lines_needed = length(right_labels) - 1;
                    for k=trial_order
                        %For every trial

                        if k == 1000
                            %line(ax(2), [other_lick_on;other_lick_on], [zeros(size(other_lick_on))+0.1+air_count-1;ones(size(other_lick_on))-0.1+air_count-1],'Color','k')
                            %line(ax(1), [other_lick_on;other_lick_on], [zeros(size(other_lick_on))+0.1+water_count-1;ones(size(other_lick_on))-0.1+water_count-1],'Color','k')
                            if u == 1
                                if left_lines_needed <= left_line_num
                                    continue
                                end
                                bigness = water_count - 1;
                            else
                                if right_lines_needed <= right_line_num
                                    continue
                                end
                                bigness = air_count - 1;
                            end

                            if bigness == 0
                                continue
                            end
                            line(ax(2), [0 - pre_time, post_time], [bigness, bigness], 'Color','black','LineStyle','--', 'LineWidth', 3);
                            line(ax(1), [0 - pre_time, post_time], [bigness, bigness], 'Color','black','LineStyle','--', 'LineWidth', 3);
                            %water_count = water_count + 1;
                            %air_count = air_count + 1;
                            left_line_num = left_line_num + 1;
                            right_line_num = right_line_num + 1;
                            continue
                        end

                        type = trial_types(k);
                        if and(type == 0, u == 2)
                            continue
                        elseif and(type == 1, u == 1)
                            continue
                        end
                        if k > 1
                            if isfield(SessionData,'TrialSettings') %(k - 1).GUI, 'LeftAmount')
                                LeftAmount = SessionData.TrialSettings(k - 1).GUI.LeftAmount;
                                RightAmount = SessionData.TrialSettings(k - 1).GUI.RightAmount;
                            elseif isfield(SessionData, 'LeftAmount')
                                LeftAmount = SessionData.LeftAmount;
                                RightAmount = SessionData.RightAmount;
                            else
                                error('Using old protocol');
                            end
                        else
                            LeftAmount = first_volume;
                            RightAmount = first_volume;
                        end
                        if and(~type, 1) %If you're checking a left dispensing trial
                        %count up each time the mouse licked left
                            if combine_spouts
                                res = find(both_event_trial_ind == k); %Find index of trial in left list
                            else
                                res = find(left_event_trial_ind == k); %Find index of trial in left list
                            end
                            if combine_spouts
                                offsets = both_sec_offsets;
                            else
                                offsets = left_sec_offsets;
                            end
                            lick_event = RawEvent{k}.Events;
                            %right_licks = RawEvent{k}.Events;
                            %% Plot what's going on on the right spout
                            if or(LeftAmount == true_volume, true_volume == 1000)
                                if plot_non_dispensing_spout    
                                    other_spout_licks = params_to_check{4};
                                    if isfield(lick_event,other_spout_licks)
                                        evaltext = ['other_lick_on = lick_event.',other_spout_licks,';'];
                                        eval(evaltext);
                                        if isempty(res)
                                            %count = count + 1; %% CHANGE TO ALLOW
                                            %EMPTY LICK TRIALS
                                            %water_count = water_count + 1;
                                            continue;
                                        end
                                        other_lick_on = other_lick_on + offsets(res) + (trial_start_times(k) - event_times(res));
                                        other_lick_on = other_lick_on - setback;
                                        if u == 1
                                            left_licks_across_all_trials{2, water_count} = other_lick_on;
                                        else
                                            right_licks_across_all_trials{2, water_count} = other_lick_on;
                                        end
                                        line(ax(1), [other_lick_on;other_lick_on], [zeros(size(other_lick_on))+0.1+water_count-1;ones(size(other_lick_on))-0.1+water_count-1],'Color','k')
                                        %water_count = water_count + 1;
                                    else
                                        if u == 1
                                            left_licks_across_all_trials{2, water_count} = [];
                                        else
                                            right_licks_across_all_trials{2, water_count} = [];
                                        end
                                    end
                                end

                                %% Plot what's going on on the left spout
                                if isfield(lick_event,spout_licks)
                                    evaltext = ['lick_on = lick_event.',spout_licks,';'];
                                    eval(evaltext);

                                    %if or(LeftAmount == true_volume, true_volume == 1000)
                                        if isempty(res)
                                            %count = count + 1;  %% CHANGED TO
                                            %ALLOW EMPTY LICK TRIALS
                                            %water_count = water_count + 1;
                                            continue;
                                        end
                                        %alt_dif = trial_start_times(k + 1) - trial_start_times(k);
                                        %bpod_dif = SessionData.TrialStartTimestamp(k + 1) - SessionData.TrialStartTimestamp(k);
                                        %licking_offset = bpod_dif - alt_dif;
                                      %disp(horzcat('Dif between alt start times: ', num2str(alt_dif)));
                                      %disp(horzcat('Dif between bpod start times: ', num2str(bpod_dif)));
                                      %disp(horzcat('Dif between bpod and analog between trials: ', num2str(licking_offset)));
                                      %disp(horzcat('Dif between bpod and analog within trial: ', num2str(trial_start_times(k) - (SessionData.TrialStartTimestamp(k) - SessionData.TrialStartTimestamp(1))))); %+ 
                                      %disp(event_times(res));
                                        lick_on = lick_on + offsets(res) + (trial_start_times(k) - event_times(res));
                                        lick_on = lick_on - setback;
                                        if u == 1
                                            left_licks_across_all_trials{1, water_count} = lick_on;
                                        else
                                            right_licks_across_all_trials{1, water_count} = lick_on;
                                        end
                                        line(ax(2), [lick_on;lick_on], [zeros(size(lick_on))+0.1+water_count-1;ones(size(lick_on))-0.1+water_count-1],'Color','k')
                                        water_count = water_count + 1;
                                else
                                    if u == 1
                                        left_licks_across_all_trials{1, water_count} = [];
                                    else
                                        right_licks_across_all_trials{1, water_count} = [];
                                    end
                                end
                                count = count + 1;

                            else
                                if u == 1
                                    %left_licks_across_all_trials{1, water_count} = [];
                                else
                                    %right_licks_across_all_trials{1, water_count} = [];
                                end

                            end

                        elseif and(type, 1) %right
                            %count up each time the mouse licked right
                            if ~combine_spouts
                                offsets = right_sec_offsets;
                                res = find(right_event_trial_ind == k); %Find index of trial in left list
                            else
                                offsets = both_sec_offsets;
                                res = find(both_event_trial_ind == k); %Find index of trial in left list
                            end
                            lick_event = RawEvent{k}.Events;
                            %right_licks = RawEvent{k}.Events;
                            %% Plot what's going on on the left spout
                            if or(RightAmount == true_volume, true_volume == 1000)
                                if plot_non_dispensing_spout    
                                    other_spout_licks = params_to_check{2};
                                    if isfield(lick_event,other_spout_licks)
                                        evaltext = ['other_lick_on = lick_event.',other_spout_licks,';'];
                                        eval(evaltext);
                                        if isempty(res)
                                            %count = count + 1; %%CHANGED
                                            %air_count = air_count + 1;
                                            continue;
                                        end
                                        other_lick_on = other_lick_on + offsets(res) + (trial_start_times(k) - event_times(res));
                                        other_lick_on = other_lick_on - setback;
                                        if u == 1
                                            left_licks_across_all_trials{1, air_count} = other_lick_on;
                                        else
                                            right_licks_across_all_trials{1, air_count} = other_lick_on;
                                        end
                                        line(ax(2), [other_lick_on;other_lick_on], [zeros(size(other_lick_on))+0.1+air_count-1;ones(size(other_lick_on))-0.1+air_count-1],'Color','k')
                                        %water_count = water_count + 1;
                                    else
                                        if u == 1
                                            left_licks_across_all_trials{1, air_count} = [];
                                        else
                                            right_licks_across_all_trials{1, air_count} = [];
                                        end
                                    end
                                end

                                %% Plot what's going on on the right spout
                                if isfield(lick_event,spout_licks)
                                    evaltext = ['lick_on = lick_event.',spout_licks,';'];
                                    eval(evaltext);

                                    %if or(LeftAmount == true_volume, true_volume == 1000)
                                        if isempty(res)
                                            %count = count + 1; %%CHANGED
                                            %air_count = air_count + 1;
                                            continue;
                                        end

                                        lick_on = lick_on + offsets(res) + (trial_start_times(k) - event_times(res));
                                        lick_on = lick_on - setback;
                                        if u == 1
                                            left_licks_across_all_trials{2, air_count} = lick_on;
                                        else
                                            right_licks_across_all_trials{2, air_count} = lick_on;

                                        end
                                        line(ax(1), [lick_on;lick_on], [zeros(size(lick_on))+0.1+air_count-1;ones(size(lick_on))-0.1+air_count-1],'Color','k')
                                        air_count = air_count + 1;
                                else
                                    if u == 1
                                        left_licks_across_all_trials{2, air_count} = [];
                                    else
                                        right_licks_across_all_trials{2, air_count} = [];
                                    end
                                end
                                count = count + 1;

                            else
                                if u == 1
                                    %left_licks_across_all_trials{2, air_count} = [];
                                else
                                    %right_licks_across_all_trials{2, air_count} = [];
                                end
                            end
                            %{
                            %% Plot what's going on on the right spout
                            lick_event = RawEvent{k}.Events;
                            if isfield(lick_event,spout_licks)
                                evaltext = ['lick_on = lick_event.',spout_licks,';'];
                                eval(evaltext);

                                if or(RightAmount == true_volume, true_volume == 1000)
                                    res = find(right_event_trial_ind == k);
                                    if isempty(res)
                                        count = count + 1;
                                        air_count = air_count + 1;
                                        continue;
                                    end
                                    if and(u == 1, ~plot_non_dispensing_spout)
                                        offsets = left_sec_offsets;
                                        count = count + 1;
                                        air_count = air_count + 1;
                                        continue
                                        %% For now, skip plotting licks on non dispensing spout
                                    elseif and(u == 1, plot_non_dispensing_spout)
                                        %Plotting left plus right
                                        offsets = right_sec_offsets;
                                    elseif u == 1 %plotting just left stuff
                                        offsets = left_sec_offsets;
                                    else
                                        offsets = right_sec_offsets;
                                    end
                                    if u == 2
                                        lick_on = lick_on - offsets(res);
                                        licks_across_all_trials{2,air_count} = lick_on;
                                        line(ax(1), [lick_on;lick_on],[zeros(size(lick_on))+0.1+air_count-1;ones(size(lick_on))-0.1+air_count-1],'Color','k')
                                        air_count = air_count + 1;
                                    elseif plot_non_dispensing_spout
                                        lick_on = lick_on - offsets(res);
                                        licks_across_all_trials{2,air_count} = lick_on;
                                        line(ax(1), [lick_on;lick_on],[zeros(size(lick_on))+0.1+air_count-1;ones(size(lick_on))-0.1+air_count-1],'Color','k')
                                        air_count = air_count + 1;
                                    else
                                        %Just don't plot?
                                    end
                                else
                                    licks_across_all_trials{2,air_count} = [];

                                end
                            else
                                licks_across_all_trials{2,air_count} = [];
                                %line(ax(1), [lick_on;lick_on],[zeros(size(lick_on))+0.1+air_count-1;ones(size(lick_on))-0.1+air_count-1],'Color','k')
                                %air_count = air_count + 1;
                            end
                            count = count + 1;
                            %}
                        end
                    end
                    if u == 1
                        licks_per_session = left_licks_across_all_trials;
                    else
                        licks_per_session = right_licks_across_all_trials;
                    end
                    tit = 'Left';
                    if and(or(~isempty(both_event_trial_ind), ~isempty(left_event_trial_ind)), u == 1)
                        %If we're plotting left trials and all spouts, plot right
                        %too
                        before_ten = [];
                        left_lps_first_5 = [];
                        left_lps_second_5 = [];
                        left_lps_third_5 = [];
                        left_lps_pre_5 = [];
                        for t = 1:size(left_licks_across_all_trials, 2)
                            left_first = left_licks_across_all_trials{1, t};
                            which = find(left_first > (0 - .05));
                            first_5 = find(left_first > 0);
                            second_5 = find(left_first > 5);
                            third_5 = find(left_first > 10);
                            pre_5 = find(left_first > -5);

                            left_first_5 = left_first(first_5);
                            left_second_5 = left_first(second_5);
                            left_third_5 = left_first(third_5);
                            left_pre_5 = left_first(pre_5);

                            first_5 = find(left_first_5 < 5);
                            second_5 = find(left_second_5 < 10);
                            third_5 = find(left_third_5 < 15);
                            pre_5 = find(left_pre_5 < 0);

                            left_first_5 = left_first_5(first_5);
                            left_second_5 = left_second_5(second_5);    
                            left_third_5 = left_third_5(third_5);
                            left_pre_5 = left_pre_5(pre_5);

                            left_first = left_first(which);
                            if ~isempty(left_first)
                                left_first = left_first(1) - 0;
                                before_ten = horzcat(before_ten, left_first);
                                left_lps_first_5 = horzcat(left_lps_first_5, length(left_first_5)/5);
                                left_lps_second_5 = horzcat(left_lps_second_5, length(left_second_5)/5);
                                left_lps_third_5 = horzcat(left_lps_third_5, length(left_third_5)/5);
                                left_lps_pre_5 = horzcat(left_lps_pre_5, length(left_pre_5)/5);
                            end
                        end
                        mean_left_lps_first_5 = mean(left_lps_first_5);
                        mean_left_lps_second_5 = mean(left_lps_second_5);
                        mean_left_lps_third_5 = mean(left_lps_third_5);
                        mean_left_lps_pre_5 = mean(left_lps_pre_5);
                        if ~isempty(vol_name)
                            vol_tit = vol_name;
                        else
                            vol_tit = 'all';
                        end
                        save(horzcat(folderInfo, '/left lick data.mat'), 'left_lps_pre_5', 'left_lps_first_5', 'left_lps_second_5', 'left_lps_third_5', 'mean_left_lps_pre_5', 'mean_left_lps_first_5', 'mean_left_lps_second_5', 'mean_left_lps_third_5');
                        total_left_licks_first_5 = length(left_lps_first_5);
                        total_left_licks_first_5_per_trial = total_left_licks_first_5/(water_count - 1);
                        left_lps_first_5 = total_left_licks_first_5_per_trial/5;
                        left_to_setback = mean(before_ten);
                        left_to_setback = 0;% - SessionData.TrialStartTimestamp(1);
                        %left_to_setback = 0 - left_to_setback;
                        if plot_non_dispensing_spout
                            set(ax(1),'TickDir', 'out','xlim',[0 - pre_time + left_to_setback post_time + left_to_setback],'xtick',(0 - pre_time + left_to_setback):2:(post_time + left_to_setback), 'xticklabel',(0 - pre_time):2:post_time,'ylim',[0,water_count - 1], 'ytick',1:2:(water_count - 1), 'yticklabel', flip(0:2:(water_count - 1)))
                            xlabel(ax(1), 'Time (s)','FontSize', 20)
                            ylabel(ax(1), 'Trial #','FontSize', 18)
                            hold(ax(1), 'on');
                            plot(ax(1), [left_to_setback,left_to_setback],[0,water_count - 1], 'color', 'red')
                            hold(ax(1), 'off');

                            set(ax(2),'TickDir', 'out','xlim',[0 - pre_time + left_to_setback post_time + left_to_setback],'xtick',(0 - pre_time + left_to_setback):2:(post_time + left_to_setback), 'xticklabel',(0 - pre_time):2:post_time,'ylim',[0,water_count - 1], 'ytick',1:2:(water_count - 1), 'yticklabel', flip(0:2:(water_count - 1)));
                            xlabel(ax(2), 'Time (s)','FontSize', 20)
                            ylabel(ax(2), 'Trial #','FontSize', 18)
                            hold(ax(2), 'on');
                            plot(ax(2), [left_to_setback,left_to_setback],[0,water_count - 1], 'color', 'red')
                            hold(ax(2), 'off');
                        else
                            set(ax(1),'TickDir', 'out','xlim',[0 - pre_time + left_to_setback post_time + left_to_setback],'xtick',(0 - pre_time + left_to_setback):2:(post_time + left_to_setback), 'xticklabel',(0 - pre_time):2:post_time,'ylim',[0,water_count - 1], 'ytick',1:2:(water_count - 1), 'yticklabel', flip(0:2:(water_count - 1)));
                            xlabel(ax(1), 'Time (s)','FontSize', 20)
                            ylabel(ax(1), 'Trial #','FontSize', 18)
                            hold(ax(1), 'on');
                            plot(ax(1), [left_to_setback, left_to_setback],[0,water_count - 1], 'color', 'red')

                            ax(2) = [];
                            close(h(2));
                            h(2) = [];
                        end
                    end

                    %% TODO these axes were deleted. Just save water instead
                    tit = 'Reward';
                    if and(or(~isempty(right_event_trial_ind), ~isempty(both_event_trial_ind)), u == 2)
                        %If we're plotting left trials and all spouts
                        if plot_non_dispensing_spout
                            if air_count == 1
                                air_count = 2;
                            end
                            %Label left
                            before_ten = [];
                            right_lps_first_5 = [];
                            right_lps_second_5 = [];
                            right_lps_third_5 = [];
                            right_lps_pre_5 = [];
                            for t = 1:size(right_licks_across_all_trials, 2)

                                right_first = right_licks_across_all_trials{2, t};
                                which = find(right_first > (0 - .05));
                                first_5 = find(right_first > 0);
                                second_5 = find(right_first > 5);
                                third_5 = find(right_first > 10);
                                pre_5 = find(right_first > -5);

                                right_first_5 = right_first(first_5);
                                right_second_5 = right_first(second_5);
                                right_third_5 = right_first(third_5);
                                right_pre_5 = right_first(pre_5);

                                first_5 = find(right_first_5 < 5);
                                second_5 = find(right_second_5 < 10);
                                third_5 = find(right_third_5 < 15);
                                pre_5 = find(right_pre_5 < 0);

                                right_first_5 = right_first_5(first_5);
                                right_second_5 = right_second_5(second_5);    
                                right_third_5 = right_third_5(third_5);
                                right_pre_5 = right_pre_5(pre_5);
                                right_first = right_first(which);

                                if ~isempty(right_first)
                                    right_first = right_first(1) - 0;
                                    before_ten = horzcat(before_ten, right_first);
                                    %right_licks_first_5 = horzcat(right_licks_first_5, right_first_5);
                                    right_lps_first_5 = horzcat(right_lps_first_5, length(right_first_5)/5);
                                    right_lps_second_5 = horzcat(right_lps_second_5, length(right_second_5)/5);
                                    right_lps_third_5 = horzcat(right_lps_third_5, length(right_third_5)/5);
                                    right_lps_pre_5 = horzcat(right_lps_pre_5, length(right_pre_5)/5);
                                end
                            end
                            mean_right_lps_first_5 = mean(right_lps_first_5);
                            mean_right_lps_second_5 = mean(right_lps_second_5);
                            mean_right_lps_third_5 = mean(right_lps_third_5);
                            mean_right_lps_pre_5 = mean(right_lps_pre_5);

                            if ~isempty(vol_name)
                                vol_tit = vol_name;
                            else
                                vol_tit = 'all';
                            end

                            save(horzcat(folderInfo, '/right lick data.mat'), 'right_lps_pre_5', 'right_lps_first_5', 'right_lps_second_5', 'right_lps_third_5', 'mean_right_lps_pre_5', 'mean_right_lps_first_5', 'mean_right_lps_second_5', 'mean_right_lps_third_5');
                            right_to_setback = mean(before_ten);
                            right_to_setback = 0;% - SessionData.TrialStartTimestamp(1);
                            %right_to_setback = 0 - right_to_setback;
                            set(ax(1),'TickDir', 'out','xlim',[0 - pre_time + right_to_setback post_time + right_to_setback],'xtick',(0 - pre_time + right_to_setback):2:(post_time + right_to_setback), 'xticklabel',(0 - pre_time):2:post_time,'ylim',[0,air_count - 1], 'ytick',1:2:(air_count - 1), 'yticklabel', flip(0:2:(air_count - 1)))
                            xlabel(ax(1), 'Time (s)','FontSize', 20)
                            ylabel(ax(1), 'Trial #','FontSize', 18)
                            hold(ax(1), 'on')
                            plot(ax(1), [right_to_setback,right_to_setback],[0,air_count - 1], 'color', 'red')
                            %hold off

                            set(ax(2),'TickDir', 'out','xlim',[0 - pre_time + right_to_setback post_time + right_to_setback],'xtick',(0 - pre_time + right_to_setback):2:(post_time + right_to_setback), 'xticklabel',(0 - pre_time):2:post_time,'ylim',[0,air_count - 1], 'ytick',1:2:(air_count - 1), 'yticklabel', flip(0:2:(air_count - 1)));
                            xlabel(ax(2), 'Time (s)','FontSize', 20)
                            ylabel(ax(2), 'Trial #','FontSize', 18)
                            hold(ax(2), 'on');
                            plot(ax(2), [right_to_setback,right_to_setback],[0,air_count - 1], 'color', 'red')
                        else
                            set(ax(2),'TickDir', 'out','xlim',[0 - pre_time + right_to_setback post_time + right_to_setback],'xtick',(0 - pre_time + right_to_setback):2:(post_time + right_to_setback), 'xticklabel',(0 - pre_time):2:post_time,'ylim',[0,air_count - 1], 'ytick',1:2:(air_count - 1), 'yticklabel', flip(0:2:(air_count - 1)));
                            xlabel(ax(2), 'Time (s)','FontSize', 20)
                            ylabel(ax(2), 'Trial #','FontSize', 18)

                            ax(1) = [];
                            close(h(1));
                            h(1) = [];
                        end
                    end

                    if 1
                        savefig(h, horzcat(folderInfo, '/Lick_Data_Plots.fig'));
                        %left then right spout. licks_across_all_trials
                        if u == 1
                            left_spout = left_licks_across_all_trials(1, :);
                            to_setback = left_to_setback;
                        else
                            left_spout = right_licks_across_all_trials(1, :);
                            to_setback = right_to_setback;
                        end

                        left_licks_within_window = [];
                        for w = 1:size(left_spout, 2)
                            trial = left_spout{1, w};
                            bigs = find(trial > (left_window + to_setback));
                            smalls = find(trial(bigs) < (right_window + to_setback));
                            left_licks_within_window = horzcat(left_licks_within_window, length(smalls));
                        end

                        if combine_spouts
                            if u == 1
                                on_spout_licks = horzcat(on_spout_licks, left_licks_within_window);
                            else
                                off_spout_licks = horzcat(off_spout_licks, left_licks_within_window);
                            end
                        end

                        left_no_licks_within_window = length(find(left_licks_within_window == 0));
                        if u == 1
                            right_spout = left_licks_across_all_trials(2, :);
                        else
                            right_spout = right_licks_across_all_trials(2, :);
                        end

                        right_licks_within_window = [];
                        for w = 1:size(right_spout, 2)
                            trial = right_spout{1, w};
                            bigs = find(trial > (left_window + to_setback));
                            smalls = find(trial(bigs) < (right_window + to_setback));
                            right_licks_within_window = horzcat(right_licks_within_window, length(smalls));
                        end

                        if combine_spouts
                            if u == 2
                                on_spout_licks = horzcat(on_spout_licks, right_licks_within_window);
                            else
                                off_spout_licks = horzcat(off_spout_licks, right_licks_within_window);
                            end
                        end

                        right_no_licks_within_window = length(find(right_licks_within_window == 0));

                        if u == 1
                            on_trial_number = on_trial_number + length(left_licks_within_window) - left_no_licks_within_window;
                            trial_num = length(left_licks_within_window) - left_no_licks_within_window;
                        else
                            on_trial_number = on_trial_number + length(right_licks_within_window) - right_no_licks_within_window;
                            trial_num = length(right_licks_within_window) - right_no_licks_within_window;
                        end

                        if ~combine_spouts
                            left_lps_within_window = sum(left_licks_within_window)/trial_num;
                            avg_left_lps_within_window = left_lps_within_window/ (right_window - left_window);
                            left_licks_within_window_per_trial = left_licks_within_window;
                            left_lps_within_window_per_trial = left_licks_within_window./(right_window - left_window);
                            right_licks_within_window_per_trial = right_licks_within_window;
                            right_lps_within_window_per_trial = right_licks_within_window./(right_window - left_window);
                            %Find right lps

                            right_lps_within_window = sum(right_licks_within_window)/trial_num;
                            avg_right_lps_within_window =right_lps_within_window/ (right_window - left_window);
                            total_lps_within_window = (sum(right_licks_within_window) + sum(left_licks_within_window))/(trial_num * (right_window-left_window));
                            total_licks_within_window = sum(right_licks_within_window) + sum(left_licks_within_window);
                            total_licks_within_window_per_trial = left_licks_within_window_per_trial + right_licks_within_window_per_trial;
                            total_lps_within_window_per_trial = total_licks_within_window_per_trial./(right_window - left_window);
                            save(horzcat(folderInfo,'/lick rate data.mat'), 'avg_left_lps_within_window', 'avg_right_lps_within_window', 'left_licks_within_window_per_trial', 'right_licks_within_window_per_trial', 'left_lps_within_window_per_trial', 'right_lps_within_window_per_trial', 'total_licks_within_window_per_trial', 'total_lps_within_window_per_trial', 'total_lps_within_window', 'total_licks_within_window', 'window');
                        end

                    end
                    hold off;

                    baseline = [];
                    base = 7;%SessionData.TrialSettings(1).GUI.BaselineTime;
                    if 1
                        response_window = 2;
                    else
                        response_window = SessionData.TrialSettings(1).GUI.ResponseTimeGo;
                    end
                    air_response = [];
                    right_lick_d = {};
                    a_count = 1;
                    water_response = [];
                    left_lick_d = {};
                    w_count = 1;
                    air_trial_ind = [];
                    %This prints all trial lick graphs
                    %{
                    if strcmp(alignment, 'sound')
                        for q = 1:length(trial_types)
                            trial_licks = licks_across_all_trials{q} - 2;
                            LeftAmount = SessionData.TrialSettings(q).GUI.LeftAmount;
                            RightAmount = SessionData.TrialSettings(q).GUI.RightAmount;
                            if or(and(trial_types(q) == 1, volume_to_check == RightAmount), volume_to_check == 1000)
                                right_lick_d{a_count} = trial_licks;
                                baseline_lick = length(find(trial_licks < base));
                                right_licks = trial_licks(find(trial_licks > base));
                                air_lick_total = length(find(right_licks < (base + response_window)));
                                baseline = horzcat(baseline, baseline_lick/base);
                                air_lps = air_lick_total/response_window;
                                if air_lps > 15
                                    air_lps = 15;
                                end
                                air_response = horzcat(air_response, air_lps);
                                a_count = a_count + 1;
                            elseif or(and(trial_types(q) == 0, volume_to_check == LeftAmount), volume_to_check == 1000)
                                left_lick_d{w_count} = trial_licks;
                                baseline_lick = length(find(trial_licks < base));
                                left_licks = trial_licks(find(trial_licks > base));
                                water_lick_total = length(find(left_licks < (base + response_window)));
                                baseline = horzcat(baseline, baseline_lick/base);
                                water_lps = water_lick_total/response_window;
                                if water_lps > 15
                                    water_lps = 15;
                                end
                                water_response = horzcat(water_response, water_lps);
                                w_count = w_count + 1; 
                            end
                        end
                        %Plot lick rates line graph air
                        right_licks = cell2mat(right_lick_d);
                        right_licks(find(right_licks > seconds_per_plot)) = [];

                        right_H = hist(right_licks, seconds_per_plot*2);
                        right_H = right_H./(a_count - 1);
                        right_H = right_H * 2;
                        left_licks = cell2mat(left_lick_d);
                        left_licks(find(left_licks > seconds_per_plot)) = [];
                        figure
                        hold on
                        left_H = hist(left_licks, seconds_per_plot*2);
                        left_H = left_H./(w_count - 1);
                        left_H = left_H * 2;
                        %Plot CS time
                        x1 = [US_time*2 - .5,US_time*2 - .5+1];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Next color
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%is orange,
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%delete next
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%US, change
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%titles
                        max_licks = 7;
                        if max(left_H) > max_licks
                            max_licks = max(left_H) + 1;
                        end
                        y1 = [max_licks max_licks];
                        h=area(x1,y1,'LineStyle','None');
                        h(1).FaceColor = [0.976 0.894 0.690];
                        %plot([US_time*2 US_time*2],[0 8], 'color', 'black', 'LineWidth', 2);
                        %Plot lick data
                        plot(smooth(left_H, 3), 'LineWidth',5, 'color', 'r');
                        plot(smooth(right_H, 3), 'LineWidth',5, 'color', [0.612 0.961 0.922]);
                        set(gca,'TickDir', 'out','xlim',[1 seconds_per_plot*2],'xtick',0:2:(seconds_per_plot*2), 'xticklabel',0:1:(seconds_per_plot),'ylim',[0,max_licks],'FontSize', 10);
                        title(horzcat('Lick Rate in ', names{u}, ' spout'), 'FontSize', 20);
                        xlabel('Time (s)','FontSize', 20)
                        ylabel('Licks/second','FontSize', 18)
                        legend('Stimulus', horzcat(left_fluid, ' trials'), horzcat(right_fluid, ' trials'))

                        if 1
                            savefig(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check), '/Lick_Data_line'));
                            save(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check), '/lick_rates'), 'right_H', 'left_H', 'right_licks', 'left_licks');
                        end
                        hold off
                        %Plot lick rates
                        avg_right_response = mean(air_response);
                        if isnan(avg_right_response)
                            avg_right_response = 0;
                        end
                        avg_left_response = mean(water_response);
                        if isnan(avg_left_response)
                            avg_left_response = 0;
                        end

                        avg_baseline = mean(baseline);
                        c = categorical({'Baseline','Left Response', 'Right Response'});
                        dat = horzcat(avg_baseline, avg_left_response, avg_right_response);
                        std_base = std(baseline)/sqrt(length(baseline));
                        std_water_resp = std(water_response)/sqrt(length(water_response));
                        std_air_resp = std(air_response)/sqrt(length(air_response));
                        if isnan(std_water_resp)
                            std_water_resp = 0;
                        end
                        if isnan(std_air_resp)
                            std_air_resp = 0;
                        end
                        errhigh = [std_base, std_water_resp, std_air_resp];
                        errlow  = [-std_base, -std_water_resp, -std_air_resp]; 
                        %%TODO SOMETHING IS WRONG WITH BAR
                        bar(c,dat)                
                        hold on
                        er = errorbar(c,dat,errlow,errhigh); 

                        er.Color = [0 0 0];                            
                        er.LineStyle = 'none';  
                        title(horzcat('All', ' Trial Lick Average'), 'FontSize', 20);
                        ylabel('Licks/s','FontSize', 18);
                        ax = gca;
                        ax.FontSize = 16;
                        savefig(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check), '/Lick_rate.fig'));
                        save(horzcat('Behavior Figures/', horzcat(names{u}, ' ', fluid), '/', num2str(volume_to_check), '/lick_rate_avg.mat'), 'dat', 'errhigh');
                        hold off

                    else
                        air_trial_ind = [];
                        for ii = 1:numel(tem)
                            load(horzcat(tem(ii).folder, '/', tem(ii).name));
                        end
                    end
                    %}
                    cd(mouse_folder);
                end
            end
            if combine_spouts
                on_spout_lps = sum(on_spout_licks)/on_trial_number;
                avg_on_spout_lps = on_spout_lps./ (right_window - left_window);
                correct_ons = find(on_spout_licks ~= 0);
                on_licks_per_trial = on_spout_licks(correct_ons);
                on_lps_per_trial = on_spout_licks./(right_window - left_window);
                on_lps_per_trial = on_lps_per_trial(correct_ons);

                off_spout_lps = sum(off_spout_licks)/on_trial_number;
                avg_off_spout_lps = off_spout_lps./ (right_window - left_window);
                off_licks_per_trial = off_spout_licks(correct_ons);
                off_lps_per_trial = off_spout_licks./(right_window - left_window);
                off_lps_per_trial = off_lps_per_trial(correct_ons);
                %Find right lps

                %right_lps = sum(right_licks)/trial_num;
                %avg_right_lps =right_lps/ (right_window - left_window);
                total_lps_within_window = (sum(on_spout_licks) + sum(off_spout_licks))/(on_trial_number * (right_window-left_window));
                total_licks_within_window = sum(off_spout_licks) + sum(on_spout_licks);
                total_licks_within_window_per_trial = on_licks_per_trial + off_licks_per_trial;
                total_lps_within_window_per_trial = total_licks_within_window_per_trial./(right_window - left_window);
                save(horzcat('Behavior Figures/',horzcat(vol_name, result_name), '/', num2str(volume_to_check), ' uL', ' lick rate data.mat'), 'avg_on_spout_lps', 'on_licks_per_trial', 'on_lps_per_trial', 'avg_off_spout_lps', 'off_licks_per_trial', 'off_lps_per_trial', 'total_lps', 'total_licks', 'total_licks_per_trial', 'total_lps_per_trial');
            end
            random_message(1);
            close all
        end
    end
 end
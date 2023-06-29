%% Instructions:
%% Update: Fixed right spout lps. Added lps per trial, licks per trial. 
%saved in 'Behavior Figures' folder. Z score peaks are in 'peak graphs' 
%folder, 'first_peaks' and 'second_peaks' variables. I will make this more 
%intuitive later. 

%{
Set variables to determine how this runs.

trials_to_check: Set which trials to check. It needs to be a list from
    smallest to largest, and it can't be bigger than the amount of trials
    you recorded

fiber: This is set up to handle recording 2 brain areas. Set 1 for the
    first fiber and 2 for the second.

correct_for_bleaching: 1 to correct for bleaching, 0 to not. What this does
    is it fits an exponential curve to your entire dataset. Then it
    subtracts that curve from your data to get values hovering around 0.
    Then it adds back the average fluorescence from all trials. I recommend
    using 1.

moving_avg: Corrects for bleaching in the way Kai suggested. Takes a window
    and slides it along data, changing F value gradually


Written by Sara Boyle

if ~exist(horzcat(file, fiber_name), 'dir')
    mkdir(horzcat(file, fiber_name));
    mkdir(horzcat(file, fiber_name, '/', alignment));
elseif ~exist(horzcat(file, fiber_name, '/', alignment))
    mkdir(horzcat(file, fiber_name, '/', alignment));
end
filesub = horzcat(file, fiber_name, '/', alignment, '/');
cd(filesub);

%}

function  Photometry_Analysis_V11(to_analyze, right_spout, left_spout, left_lick_port, right_lick_port, pre_time, post_time, seconds_per_trial, seconds_to_extend, trials_to_check, fibers, fiber_names, volume_to_check, pre_ITI_length, trigger1, frame_rate, leave_out_no_licks, lick_window, signal_key_word, analog_key_word, bpod_key_word, moving_avg, correct_for_bleaching)
    %% Clear Workspace
    %clear, clc
    if trigger1
        frame_rate = frame_rate/2;
    end
    
    file = pwd;
    file = horzcat(file, '/');
    
    for fiber = fibers
        for signal = 1:2
        close all

        %% Set up Variables

        %Variables to select what you want to do
        stim_duration = .2;
        subtract_autofluorescence = 1; % if 0, don't subtract. If 1, choose the values to subtract carefully
        plot_heatmap_and_average = 1; %Haven't worked this out when water_align variable is on 

        %What type of data do you have?
        stim_time = 7; %doesn't matter

        %Never change below here
        if signal == 1
            if volume_to_check == 1000
                fiber_name = horzcat(fiber_names{fiber}, ' signal');
            else
                fiber_name = horzcat(fiber_names{fiber}, ' ', num2str(volume_to_check), ' uL signal');
            end
        else
            if volume_to_check == 1000
                fiber_name = horzcat(fiber_names{fiber}, ' control');
            else
                fiber_name = horzcat(fiber_names{fiber}, ' ', num2str(volume_to_check), ' uL control');
            end
        end

        if volume_to_check == 1000
            vol_name = '';
        else
            vol_name = horzcat(num2str(volume_to_check), ' ');
        end

        %% Before the for loop   
        
        if signal == 1
            %Load data
            [SessionData, data, start_times, alt_start_times] = load_data(file, bpod_key_word, signal_key_word, analog_key_word, 1, 1, 1);

            [data, auto_flo] = trim_data(data, fiber, subtract_autofluorescence);

            %De-interleave and subtract autofluorescence
            [gcamp, iso] = de_interleave(trigger1, 0, file, fiber_name, data, fiber, auto_flo, subtract_autofluorescence, signal);

            
            
            
            
            
            raw_trials = SessionData.RawEvents.Trial;
            pooled_delivery_times = [];
            for trial_ind = 1:size(raw_trials, 2)
                if ~isnan(raw_trials{1, trial_ind}.States.DeliverLeft(1))
                    delivery = raw_trials{1, trial_ind}.States.DeliverLeft(1) + start_times(trial_ind);
                    pooled_delivery_times = horzcat(pooled_delivery_times, delivery);
                elseif ~isnan(raw_trials{1, trial_ind}.States.DeliverRight(1))
                    delivery = raw_trials{1, trial_ind}.States.DeliverRight(1) + start_times(trial_ind);
                    pooled_delivery_times = horzcat(pooled_delivery_times, delivery);
                end
            end
            
            
            
            
            [gcamp, iso, baseline_gcamp, baseline_iso] = correct_bleaching(fiber_name, correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, trigger1, 0, pooled_delivery_times);
        else
            temp_gcamp = gcamp;
            gcamp = iso;
            iso = temp_gcamp;
        end
        %Is sound on?
        sound_on = 0;
        %{
        %sound_on = input('Is the CS on? Enter for yes, n for no, e for exit. ', 's');

        if strcmp(sound_on, 'n')
            sound_on = 0;
        elseif isempty(sound_on)
            sound_on = 0; %sound CHANGED BECAUSE ALE
            %align_to_US = 0;
        elseif sound_on == 'e'
            disp('User exited the program');
            return
        end
        %}
        sound_only = 0;
        man_target_frame = seconds_per_trial + seconds_to_extend;
        man_target_frame = frame_rate * man_target_frame/2;
        man_target_frame = 100;
        [first_volume, left_recorded_options, right_recorded_options] = find_first_volume(SessionData);
        total_volumes = horzcat(left_recorded_options, right_recorded_options);
        for iteration = 1:size(to_analyze, 2)
            %%%%%%%%%%%%%%% TODO ADD A FOR LOOP FOR ALL THE VOLUMES HERE
            for single_volume = unique(total_volumes)
            %% Don't change anything below here
            
            
            
            
            
            
            if signal == 1
                if single_volume == 1000
                    fiber_name = horzcat(fiber_names{fiber}, ' signal');
                else
                    fiber_name = horzcat(fiber_names{fiber}, ' ', num2str(single_volume), ' uL signal');
                end
            else
                if single_volume == 1000
                    fiber_name = horzcat(fiber_names{fiber}, ' control');
                else
                    fiber_name = horzcat(fiber_names{fiber}, ' ', num2str(single_volume), ' uL control');
                end
            end
            
            
            
            
            
            
            
            
            
            alignment = to_analyze{iteration};
            align_to_US = 1;

            if or(strcmp(alignment,'left'), strcmp(alignment,'right'))
                %To look at water trials
                water_on = 1; 
                air_on = 0;
            elseif strcmp(alignment,'both')
                water_on = 1; 
                air_on = 0;
            end  
            %Variables for saving data
            save_graphs = 1; 
            save_data = 1;
            create_save_folder(fiber_name, alignment, file, right_spout, left_spout, single_volume);
            mode = horzcat(num2str(sound_on), num2str(water_on), num2str(air_on));
            %Record time of each stimulus
            %if ~isnan(SessionData.Outcomes(1))
                %outcomes = SessionData.Outcomes;
            %else
                %outcomes = ones(1, size(trial_types, 2));
            %end
            if isfield(SessionData, 'TrialTypes')
            if ~isempty(SessionData.TrialTypes)
                trial_types = SessionData.TrialTypes;
            else
                trial_types = 1:SessionData.nTrials;
                trial_types = mod(trial_types, 2);
                outcomes = trial_types;
            end
            elseif isfield(SessionData, 'spout')
                spout = SessionData.spout;
                if spout == 3
                    trial_types = 1:SessionData.nTrials;
                    trial_types = mod(trial_types, 2);
                    outcomes = trial_types;
                elseif spout == 2
                    trial_types = ones(1, SessionData.nTrials);
                elseif spout == 1
                    trial_types = zeros(1, SessionData.nTrials);
                end
            end

            outcomes = ones(1, size(trial_types, 2)); %added

            [event_times, event_end_times, event_trial_ind, no_event_trial_ind, event_labels] = find_event_times(SessionData, alt_start_times, alignment, trials_to_check, single_volume, single_volume, leave_out_no_licks, stim_time, sound_on, outcomes, left_lick_port, right_lick_port, lick_window);
             if isempty(event_trial_ind)
                 disp(horzcat(alignment, ' ', num2str(volume_to_check),' uL trials do not exist'))
                 cd(file);
                 if strcmp(alignment, 'left')
                     rmdir(horzcat(file, fiber_name, '/left ', left_spout));
                 else
                     rmdir(horzcat(file, fiber_name, '/right ', right_spout));
                 end
                 continue
             end

            %% De-interleave and subtract autofluorescence      
            gcamp_clone = gcamp;

            %% Align trial data and prepare for display
            end_times = [];
            for x = trials_to_check
                single_trial = SessionData.RawEvents.Trial{1, x};
                %end_times = horzcat(end_times, single_trial.States.EndTime(2) + alt_start_times(x)); %- SessionData.TrialStartTimestamp(1));
            end

            %% Dividing trials based on start times
            if iscell(event_times)
                size_to_check = size(event_times, 2);
            else
                size_to_check = size(event_times, 1);
            end

            if and(size_to_check > 1, ~strcmp(alignment, 'right'))
                sound_delivered_times = event_times{2};
                event_times = event_times{1};
            else
                sound_delivered_times = [];
                event_times = event_times{1};
            end
            %for all trials, not just US receiving
            trials = {};
            iso_trials = {};
            OG_trials = {};
            pre_ITIs = {};
            post_ITIs = {};
            US_frame = {}; 
            no_US_delivered = [];
            US_delivered = [];
            ind = 1;
            US_ind = 1;
            US_frames = [];
            CS_frames = [];
            checked_trial_types = [];
            other_side = [];
            %Just for US receiving
            event_trials = {};
            event_i = 1;
            event_pre_ITIs = {};
            event_post_ITIs = {};
            event_times_temp = event_times;
            event_frame = {};          %Record frame air puff happens                              
            event_delivered = [];

            %Just for no US trials
            no_event_trials = {};
            no_event_i = 1;
            no_event_pre_ITIs = {};
            no_event_post_ITIs = {};
            sound_delivered_times_temp = sound_delivered_times;
            %start_times = alt_start_times; 
            start_times = event_times - pre_time;
            end_times = event_times + post_time;
            evs = 1;
            %TODO INVESTIGATE
             for y = 1:length(start_times)
                %Segment the trials based on start and end times
                %pre_ITI_length = 5;
                post_ITI_length = 150;   
                if evs > length(event_times)
                    continue
                end
                starty = start_times(y);
                endy = end_times(y);
                indices1 = find(gcamp(:,1) >= starty);
                if trigger1
                    indices2 = find(iso(:,1) >= starty);
                    trial2= iso(indices2,:);
                end
                indices3 = find(gcamp_clone(:,1) >= starty);
                trial1=gcamp(indices1,:);
                trial3 = gcamp_clone(indices3, :);

                indices1 = find(trial1(:,1) <= endy);
                if trigger1
                    indices2 = find(trial2(:,1) <= endy);
                    trial2=trial2(indices2,:);
                end
                indices3 = find(trial3(:,1) <= endy);
                trial1=trial1(indices1,:);
                trial3=trial3(indices3,:);
                pre_ITI=gcamp(find(gcamp(:,1) >= (start_times(y) - pre_ITI_length)),:); %Bigger than start - 15
                pre_ITI=pre_ITI(find(pre_ITI(:,1) <= (start_times(y))), :); %smaller than start of trial?
                post_ITI= gcamp(find(gcamp(:,1) >= end_times(y)), :); %Greater than end
                post_ITI_clone= post_ITI(find(post_ITI(:,1) <= (end_times(y) + post_ITI_length)),:);  %Less than end + ITI
                post_ITI= post_ITI(find(post_ITI(:,1) <= (end_times(y) + post_ITI_length)),:);  %Less than end + ITI

                US_ind = event_i;
                if ~isempty(event_times_temp)            %If not too many USs
                    if y ~= 1
                        if  and(event_times_temp(1) < starty, event_times_temp(1) > end_times(y - 1) + seconds_to_extend)
                            disp(horzcat('The mouse licked too late into trial ', num2str(y - 1), '. Skipping that event'));
                            event_times_temp = event_times_temp(2:end);
                            to_take_out = find(event_trial_ind == y - 1);
                            event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                            event_trial_ind(to_take_out) = [];
                            ind = ind - 1;
                        end
                    elseif event_times_temp(1) > end_times(1) + seconds_to_extend
                        disp(horzcat('The mouse licked too late into trial ', num2str(y - 1), '. Skipping that event'));
                        %event_times_temp = event_times_temp(2:end);
                        %to_take_out = find(event_trial_ind == y);
                        %event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                        %event_trial_ind(to_take_out) = [];
                        %ind = ind - 1;
                    end
                    US_frames = find(trial1(:, 1) >= event_times_temp(1)); %Find US frames
                     if isempty(US_frames)
                        US_frames = find(post_ITI_clone(:, 1) >= event_times_temp(1)); %Find US frames
                        US_frames = US_frames + size(trial1, 1);
                        if ~isempty(US_frames)
                            if US_frames(1) > size(trial1,1) + seconds_to_extend*10 -20
                                disp(horzcat('The mouse licked too late into trial ', num2str(event_trial_ind(y)), '. Skipping that event'));
                                event_times_temp = event_times_temp(2:end);
                                to_take_out = find(event_trial_ind == event_trial_ind(y));
                                event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                                event_trial_ind(to_take_out) = [];
                                ind = ind - 1;
                                US_frames = find(trial1(:, 1) >= event_times_temp(1)); %Find US frames
                                if isempty(US_frames)
                                    disp('This mouse licked too late into previous trial then too late into current')
                                end
                            end
                        end
                    end
                end

                if and(sound_on == 1, ~isempty(sound_delivered_times_temp))
                    %Find CS_frames
                    sound_time = sound_delivered_times_temp(1);
                    if sound_time ~= 0
                        CS_frames = find(trial1(:, 1) >= sound_time);
                        sound_delivered_times_temp = sound_delivered_times_temp(2:end); %Delete CS
                    else
                        sound_delivered_times_temp = sound_delivered_times_temp(2:end); %Delete CS
                    end
                end

                if ~isempty(CS_frames)                               %If there's a CS
                    sound_frame{ind} = CS_frames(1,1);
                else
                    sound_frame{ind} = 0;
                end
                CS_frames = [];

                if or(~isempty(US_frames), sound_only)                               %If there's a US
                    if sound_only
                        US_frame{ind} = sound_frame{ind};
                    else
                        if US_frames(1, 1) ~= 1
                            US_frame{ind} = US_frames(1, 1);                 %Find the US onset frame
                        else
                            US_frame{ind} = 0;
                        end
                    end                
                    US_delivered = horzcat(US_delivered, event_trial_ind(y));         %record trial # with US
                    trials{ind} = trial1(:, [1 2]);                  %Record trial data in total list
                    if trigger1
                        iso_trials{ind} = trial2(:, [1 2]);
                    end
                    OG_trials{ind} = trial3(:, [1 2]);
                    pre_ITIs{ind} = pre_ITI;                         %Record pre ITI data
                    post_ITIs{ind} = post_ITI;                       %Record post ITI data

                    if ismember(event_trial_ind(y), US_delivered)
                        event_trials{event_i} = trial1(:, [1 2]);        %Record trial as air trial
                        event_pre_ITIs{event_i} = pre_ITI;               %Record air pre ITI
                        event_post_ITIs{event_i} = post_ITI;             %Record are post ITI
                        event_frame{event_i} = US_frame{ind};          %Record frame air puff happens                              
                        event_delivered = horzcat(event_delivered, event_i);
                        event_i = event_i + 1;
                        event_times_temp = event_times_temp(2:end);

                    end
                    US_frames = [];                                  %Clear US frames
                    ind = ind + 1;
                else                                                 %If there's no US
                    US_frame{ind} = 0;
                    other_side = horzcat(other_side, event_trial_ind(y));
                    no_US_delivered = horzcat(no_US_delivered, ind);

                    if ismember(event_trial_ind(y), no_event_trial_ind)
                        no_event_trials{no_event_i} = trial1(:, [1 2]);        %Record trial as air trial
                        no_event_pre_ITIs{event_i} = pre_ITI;               %Record air pre ITI
                        no_event_post_ITIs{event_i} = post_ITI;             %Record are post ITI
                        no_event_i = no_event_i + 1;
                    end
                    trials{ind} = trial1(:, [1 2]);
                    if trigger1
                        iso_trials{ind} = trial2(:, [1 2]);
                    end
                    OG_trials{ind} = trial3(:, [1 2]);
                    checked_trial_types = horzcat(checked_trial_types, trial_types(event_trial_ind(y)));
                    pre_ITIs{ind} = pre_ITI;
                    post_ITIs{ind} = post_ITI;
                    ind = ind + 1;
                end
            end

            stimulus_times = event_times;
            check_trials = trials;
            if sound_on
                stimulus_times_temp = {sound_delivered_times(event_trial_ind), event_times};

            else
                stimulus_times_temp = {zeros(1, size(event_times, 2)), event_times};
            end
            if sound_on
                if ~isempty(no_event_trial_ind)
                    no_stimulus_times_temp = {sound_delivered_times(no_event_trial_ind), []};
                end
            else
                no_stimulus_times_temp = {zeros(1, size(event_times, 2)), event_times};
            end

            %% Sort and align trials.
            if plot_heatmap_and_average
                align_to_this = alignment;
                if ~isempty(event_trial_ind)
                    good_pre_ITIs = pre_ITIs;
                    good_post_ITIs = post_ITIs;
                    [US_trials, US_target_frame, CS_target_frame, end_CS_frame, true_seconds_per_trial, sec_offsets, true_sec_offsets] = sort_and_align_trials_V7(event_trials, good_pre_ITIs, good_post_ITIs, ...
                        seconds_per_trial, US_frame, sound_frame, align_to_this, align_to_US, frame_rate, seconds_to_extend, stim_duration, 1, pre_time*10);
                    frame_offsets = sec_offsets * frame_rate;
                    if single_volume == 1000
                        vol_name = '';
                    else
                        vol_name = horzcat(num2str(single_volume), ' ');
                    end
                    trial_start_times = alt_start_times;
                    trial_order = trials_to_check;
                    if ~contains(fiber_name, 'control')
                        save(horzcat(file, fiber_name, ' ', alignment, ' spout alignment data.mat'), 'sec_offsets', 'frame_offsets', 'event_trial_ind', 'no_event_trial_ind', 'true_sec_offsets', 'event_times', 'trial_start_times', 'trial_order');
                    end
                end
                no_US_trials = [];
                %if ~isempty(no_event_trial_ind)
                %    no_event_pre_ITIs = pre_ITIs(no_event_trial_ind);
                %    no_event_post_ITIs = post_ITIs(no_event_trial_ind); 
                %    empty_pre_ITIs = pre_ITIs(no_event_trial_ind);
                %    empty_post_ITIs = post_ITIs(no_event_trial_ind);
                %    [no_US_trials, no_US_target_frame, CS_target_frame, end_CS_frame, true_seconds_per_trial, sec_offsets] = sort_and_align_trials_V6(no_event_trials, empty_pre_ITIs, empty_post_ITIs, ...
                %        seconds_per_trial, US_frame(no_event_trial_ind), US_frame(no_event_trial_ind), align_to_this, 0, frame_rate, seconds_to_extend, stim_duration, 1, man_target_frame);
                %    %frame_offsets = sec_offsets * frame_rate;
                %    %save(horzcat(file, alignment, ' no event alignment data.mat'), 'sec_offsets', 'frame_offsets', 'no_event_trial_ind');
                %end
                %man_target_frame = US_target_frame;
                %% Plot average response
                %if no_event_trial_ind
                %    no_US_trials = no_US_trials *100;   
                %end
                if ~isempty(event_trial_ind)
                    US_trials = US_trials * 100;
                    baselines = [];
                    for r = 1:size(US_trials, 1)
                        cur = US_trials(r, :);
                        baseline = cur(10:((pre_time - 1)*frame_rate - 5));
                        baselines(r, :) = baseline;
                    end
                    plot_average(mode, fiber_name, US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                    ylabel('dF/F (%)', 'FontSize', 20);
                    if strcmp(alignment, 'left')
                        title(horzcat(fiber_name, ': % dF/F in ', left_spout, ' trials'));
                    elseif strcmp(alignment, 'right')
                        title(horzcat(fiber_name, ': ', right_spout, ' trials'));
                    else
                        title(horzcat(fiber_name, ': % dF/F in US delivered ', alignment, ' trials'));
                    end

                    if save_graphs
                        if strcmp(alignment , 'left')
                            savefig('left_average_dFF.fig')
                        elseif strcmp(alignment , 'right')
                            savefig('right_average_dFF.fig')
                        else
                            savefig('US_average.fig');
                        end
                    end
                else
                    disp(horzcat('Found no ', alignment, ' trials. Exited.'))
                    cd(file);
                    continue
                end

                %% Plot trial by trial responses
                %if and(~isempty(no_event_trial_ind), ~strcmp(alignment, 'assist'))
                    %baseline_no_US_trials = no_US_trials(:,20:(US_target_frame - 20));
                %else
                    baseline_no_US_trials = [];
                    las_baseline_no_US_trials = [];
                    no_las_baseline_no_US_trials = [];
                    z_no_US_trials = [];
                %end
                if ~isempty(event_trial_ind)            
                    plot_heatmap(US_trials, fiber_name, true_seconds_per_trial, US_target_frame, CS_target_frame, alignment);
                    title(horzcat(fiber_name, ': ', alignment, ' trials'));

                    if strcmp(alignment, 'left')
                        savefig('left_delivery_heatmap.fig');
                    elseif strcmp(alignment, 'right')
                        savefig('right_delivery_heatmap.fig');
                    else
                        savefig('US_heatmap.fig');
                    end              
                end

                increments = length(US_trials);
                if increments == 0
                    increments = length(no_US_trials);
                end
                x_axis = ((1:increments)/increments)*true_seconds_per_trial;
                if single_volume == 1000
                    if align_to_US
                        laser_time = x_axis(US_target_frame);
                    else
                        laser_time = x_axis(CS_target_frame);
                    end
                    x_axis = x_axis - laser_time;
                    if ~isempty(event_trial_ind)
                        if strcmp(alignment, 'left')
                            left_trials = US_trials;
                            save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                        elseif strcmp(alignment, 'right')
                            right_trials = US_trials;
                            save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                        else
                            save('dFF_data.mat', 'US_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                        end
                    end

                    %US
                    if isempty(event_trial_ind)
                        baseline_US_trials = [];
                    else
                        for x = 1:size(US_trials, 2)
                            if isnan(US_trials(1,x))
                                US_trials(1,x) = 0;
                            end
                        end
                        baseline_US_trials = US_trials(:,20:(US_target_frame - 20));
                    end
                    if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                        %baseline_trials = vertcat(baseline_US_trials, baseline_no_US_trials);
                        baseline_trials = baselines;
                        std_avg = mean(std(baseline_trials));
                        avg = mean(mean(baseline_trials));
                        diffs = (no_US_trials - avg);
                        z_no_US_trials = (diffs ./ std_avg);
                        diffs = (US_trials - avg);
                        z_US_trials = (diffs ./ std_avg); 
                        plot_average(mode, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                        title(horzcat(fiber_name, ': ', alignment, ' trials'));
                        ylabel('z-score(dF/F)', 'FontSize', 20);
                        if strcmp(alignment, 'left')
                            savefig('left_z_score_average.fig');
                        elseif strcmp(alignment, 'right')
                            savefig('right_z_score_average.fig');
                        else
                            savefig('US_z_score_average.fig');
                        end
                        if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                            plot_average(mode, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                            title(horzcat(fiber_name, ': ', alignment, ' trials'));
                            ylabel('z-score(dF/F)', 'FontSize', 20);
                            if strcmp(alignment, 'left')
                                z_right_trials = z_no_US_trials;
                                z_left_trials = z_US_trials;
                                save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');                                 
                            elseif strcmp(alignment, 'right')
                                z_right_trials = z_US_trials;
                                z_left_trials = z_no_US_trials;
                                save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                            else
                                save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                            end
                        end
                    else
                        %baseline_trials = vertcat(baseline_US_trials);
                        baseline_trials = baselines;
                        std_avg = mean(std(baseline_trials));
                        avg = mean(mean(baseline_trials));
                        diffs = (US_trials - avg);
                        z_US_trials = (diffs ./ std_avg); 
                        %post_laser = round(mean(end_CS_frame)) + 1;
                        %laser_effect = mean(z_US_trials(:,post_laser:post_laser + 20));
    %                    save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis');

                        plot_average(mode, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                        title(horzcat(fiber_name, ': ', alignment, ' trials'));
                        ylabel('z-score(dF/F)', 'FontSize', 20);
                        if strcmp(alignment, 'left')
                            savefig('left_z_score_average.fig');
                            save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                            z_left_trials = z_US_trials;
                            save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                        elseif strcmp(alignment, 'right')
                            savefig('right_z_score_average.fig');
                            save('dFF_data.mat', 'right_trials', 'US_target_frame', 'x_axis', 'baselines');
                            z_right_trials = z_US_trials;
                            save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                        else
                            savefig('US_z_score_average.fig');
                        end
                        if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                            plot_average(mode, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                            title(horzcat(fiber_name, ': ', alignment, ' trials'));
                            ylabel('z-score(dF/F)', 'FontSize', 20);
                            if strcmp(alignment, 'left') 
                                save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                                z_right_trials = z_no_US_trials;
                                z_left_trials = z_US_trials;
                                save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                            elseif strcmp(alignment, 'right')
                                save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                                z_right_trials = z_US_trials;
                                z_left_trials = z_no_US_trials;
                                save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                            else
                                save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                            end
                        end
                    end    
                else
                    if align_to_US
                        if ~isnan(US_target_frame)
                            laser_time = x_axis(US_target_frame);
                        else
                            laser_time = x_axis(72);
                            US_target_frame = 72;
                        end
                    else
                        laser_time = x_axis(CS_target_frame);
                    end
                    if strcmp(alignment, 'right')
                        right_trials = US_trials;
                    elseif strcmp(alignment, 'left')
                        left_trials = US_trials;
                    end
                    %x_axis = x_axis - laser_time;
                    laser_time = x_axis(US_target_frame);
                    x_axis = x_axis - laser_time;
                    if strcmp(alignment, 'right')
                        save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                    elseif strcmp(alignment, 'left')
                        save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                    else
                        save('dFF_data.mat', 'US_trials', 'no_US_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                    end

                    %US
                    if or(strcmp(alignment, 'right'), strcmp(alignment, 'left'))
                        %baseline_no_US_trials = no_US_trials(:,1:(US_target_frame - 20));
                        baseline_US_trials = US_trials(:,1:(US_target_frame - 20));
                    end
                    %baseline_trials = vertcat(baseline_US_trials, baseline_no_US_trials);
                    baseline_trials = baselines;
                    std_avg = mean(std(baseline_trials));
                    avg = mean(mean(baseline_trials));
                    %diffs = (no_US_trials - avg);
                    %z_no_US_trials = (diffs ./ std_avg);
                    diffs = (US_trials - avg);
                    z_US_trials = (diffs ./ std_avg); 
                    if ~isempty(event_trial_ind)
                        plot_average(mode, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                        if strcmp(alignment, 'left')
                            correct_title = 'left';
                        elseif strcmp(alignment, 'right')
                            correct_title = 'right';
                        else
                            correct_title = alignment;
                        end

                        title(horzcat(fiber_name, ': ', correct_title, ' trials'));
                        ylabel('z-score(dF/F)', 'FontSize', 20);
                        if save_graphs
                            if strcmp(alignment, 'left')
                                    savefig('left_z_score_average.fig');
                                    z_left_trials = z_US_trials;
                            elseif strcmp(alignment, 'right')
                                savefig('right_z_score_average.fig');
                                z_right_trials = z_US_trials;
                            else
                                savefig('US_z_score_average.fig');
                            end
                        end
                    end

                    if strcmp(alignment, 'left')
                        save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                    elseif strcmp(alignment, 'right')
                        save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                    else
                        save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                    end
                end
            end
            % It takes 0.0001 s for arduino to read analog input. .00025 s for bpod
            % state changes. .0001 to do output action
            cd(file);
            end
        random_message(2);
        end
        end
    end
end

function [] = create_save_folder(fiber_name, alignment, file, right_spout, left_spout, vol_name)
    %Sets up folders to save analyzed data, cds to new folder. Reads the bpod file in.

    %% Set up save folders
    %Create the fiber specific folders
    if ~exist(horzcat(file,fiber_name), 'dir')
        mkdir(fiber_name);
        cd(fiber_name);
    elseif exist(horzcat(file,fiber_name), 'dir')
        disp(horzcat(fiber_name, ' trials directory already exists.'))
        current_dir = pwd;
        if ~strcmp(horzcat(current_dir, '/'), horzcat(file, fiber_name, '/'))
            cd(fiber_name);
        else
            disp('You are in the correct directory');
        end
    end
    

    %Create the folder that specifies trial type 
    if strcmp(alignment, 'left')
        alignment_folder = horzcat(alignment, ' ', left_spout);
    elseif strcmp(alignment, 'right')
        alignment_folder = horzcat(alignment, ' ', right_spout);
    elseif strcmp(alignment, 'both')
        alignment_folder = horzcat(vol_name, 'both');
    end
    if ~exist(alignment_folder, 'dir')
            mkdir(alignment_folder);
            cd(alignment_folder);
        elseif exist(alignment_folder, 'dir')
            disp(horzcat(alignment_folder, ' trials directory already exists.'))
            current_dir = pwd;
            if ~strcmp(horzcat(current_dir, '/'), horzcat(file, fiber_name, '/', alignment_folder, '/'))
                cd(alignment_folder);
            else
                disp('You are in the correct directory');
            end
    end
    close all
end

function [data, auto_flo] = trim_data(data, fiber, subtract_autofluorescence)
    %% Cut data from beginning and/or end and subtract baseline fluorescence
    %Plot pre-bleach corrected data.
    if subtract_autofluorescence
        plot(data(:,1), data(:,fiber + 1));
        figure;
        plot(data(1:800,1), data(1:800,fiber + 1));
        title('First 800 seconds. Choose autofluorescence window. Click at beginning, then at end')
        [x_in,y_in]=ginput(2);
        disp(horzcat('Window is from ', num2str(x_in(1)), ' to ', num2str(x_in(2))))
        %disp(y_in)
        close all
        is_good = input('Did you choose? Press enter to continue, e to exit: ', 's');
        if is_good == 'e'
            disp('User exited. Closed program.')
            close all
            return
        end
    end
    %There's weird frame dropping
    if subtract_autofluorescence == 1
        auto = data(find(data(:,1) > x_in(1)), :);
        auto = auto(find(auto(:,1) < x_in(2)), fiber + 1);
        auto_flo = mean(auto);
        if isnan(auto)
            error('You chose autofluorescence times wrong.')
        end
        %auto_flo = 5.28;
    else 
        auto_flo = 0;
    end

    %Cut out beginning of file
    data = data(find(data(:,1) > -12), :);

    %Check if this looks good
    plot(data(:,1), data(:,fiber + 1));
    is_good = input('Does this look right? Press enter to continue without cutting any data, c to cut, e to exit: ', 's');
    if is_good == 'e'
        disp('User exited. Closed program.')
        close all
        return
    elseif is_good == 'c'
        disp('Select interval to analyze. ');
        [x_in,y_in] = ginput(2);

        %Cut anything funky from the end
        data = data(find(data(:,1) > x_in(1)), :);
        data = data(find(data(:,1) < x_in(2)), :);
    end
    close all;
end

function [event_times, event_end_times, event_trial_ind, no_event_trial_ind, event_labels] = find_event_times(SessionData, alt_start_times, alignment, trials_to_check, reward_amount, punish_amount, leave_out_no_licks, stim_time, sound_on, outcomes, left_lick_port, right_lick_port, lick_window)
    no_left_ind = [];
    left_ind = [];
    no_right_ind = [];
    right_ind = [];
    left_lick_event_times = [];
    right_lick_event_times = [];
    left_deliveries = [];
    right_deliveries = [];
    punish_i = 0;
    reward = 7;
    water_delivered_times = [];
    water_earned_times = [];
    water_earned_ind = [];
    no_water_delivered_ind = [];
    water_delivered_ind = [];
    water_assisted_times = [];
    water_assisted_ind = [];
    sound_delivered_times = [];
    sound_delivered_ind = [];
    air_sounds = [];
    air_sound_ind = [];
    water_sounds = [];
    water_sound_ind = [];
    avoid_times = [];
    avoid_ind = [];
    not_avoided = [];
    total_not_avoided = [];
    run_times = [];
    escape_times = [];
    escape_ind = [];
    laser_times = [];
    laser_ind = [];
    air_puff_times = [];
    air_puff_ind = [];
    no_air_delivered = [];
    air_puff_end_times = [];
    air_lengths = [];
    air_fail_ind = [];
    %air_trial_ind = 0;
    event_end_times = [];
    RawEvent = SessionData.RawEvents.Trial;
    [first_volume, left_recorded_options, right_recorded_options] = find_first_volume(SessionData);
    
    %% TODO FIX SO IT CAN READ DECIMALS
    %first_volume = 6;
    start_times = alt_start_times;
    no_event_trial_ind = [];
    event_trial_ind = [];
    licky = [];
    outcomes = fix_outcomes(SessionData);
    trial_types = [];
    if isfield(SessionData, 'TrialTypes')
        if ~isempty(SessionData.TrialTypes)
            trial_types = SessionData.TrialTypes;
        else
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = ~trial_types;
        end
    else
        spout = SessionData.spout;
        if spout == 3
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = ~trial_types;
        elseif spout == 2
            trial_types = ones(1, SessionData.nTrials);
            outcomes = ~trial_types;
        elseif spout == 1
            trial_types = zeros(1, SessionData.nTrials);
            outcomes = ~trial_types;
        end
    end
    %right_stim_time = SessionData.RawEvents.Trial{1, 2}.States.DeliverRight(2);
%    left_stim_time = SessionData.RawEvents.Trial{1, 2}.States.DeliverLeft(2);
    %It's not choosing the right kind of trial
    %3 receives reward with assist. 1 receive reward without assist. 0.
    %punishment. -1, avoided air puff. 2 for ignored reward. 4 escape 
    if ~strcmp(alignment, 'events')
        %trial_types = SessionData.TrialTypes;
        for q = trials_to_check
            trial_start = start_times(q);  %%%% CHANGE THIS to new start times
            if or(strcmp(alignment, 'both'), or(strcmp(alignment, 'left'), strcmp(alignment, 'right')))
                if isfield(RawEvent{q}.Events,left_lick_port)
                    left_lick_event_times = horzcat(left_lick_event_times, eval(horzcat('RawEvent{q}.Events.', left_lick_port)) + trial_start);
                    %licky = left_lick_event_times;
                end
                if isfield(RawEvent{q}.Events,right_lick_port)
                    right_lick_event_times = horzcat(right_lick_event_times, eval(horzcat('RawEvent{q}.Events.', right_lick_port)) + trial_start);
                    %licky = left_lick_event_times;
                end
            else
                if isfield(RawEvent{q}.Events,left_lick_port)
                    lick_event_times = horzcat(left_lick_event_times, eval(horzcat('RawEvent{q}.Events.', left_lick_port)) + trial_start);
                end
            end
            if strcmp(alignment, 'left')
                reward_state = 'DeliverLeft';
            elseif strcmp(alignment, 'right')
                reward_state = 'DeliverRight';
            elseif strcmp(alignment, 'both')
            else
                reward_state = 'DeliverReward';
            end
            %Record left delivery times
            if trial_types(q) == 0
                if or(strcmp(alignment, 'left'), strcmp(alignment, 'right'))
                    reward_state = 'DeliverLeft';
                elseif strcmp(alignment, 'both')
                    reward_state = 'DeliverLeft';
                end
                if isfield(SessionData.RawEvents.Trial{1, q}.States, reward_state)
                    if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)')))
                        reward_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)'));
                        if q ~= 1
                            if isfield(SessionData, 'TrialSettings')
                                if isfield(SessionData.TrialSettings(q - 1).GUI, 'RewardAmount')
                                    reward = SessionData.TrialSettings(q - 1).GUI.RewardAmount;
                                elseif isfield(SessionData.TrialSettings(q - 1).GUI, 'LeftAmount')
                                    reward = SessionData.TrialSettings(q - 1).GUI.LeftAmount;
                                else
                                    reward = first_volume;
                                end
                            else
                                if strcmp(reward_state, 'DeliverLeft')
                                    reward = SessionData.LeftAmount;
                                else
                                    reward = SessionData.RightAmount;
                                end
                            end
                        else
                            reward = first_volume;
                        end
                        if and(~isnan(reward_time), or(reward == reward_amount, reward_amount == 1000))
                            if leave_out_no_licks
                                %Find trials with no licking after the US
                                licky = left_lick_event_times - trial_start;
                                left_stim_time = SessionData.RawEvents.Trial{1, q}.States.DeliverLeft(1);
                                if isfield(SessionData.RawEvents.Trial{1, q}.States, 'LickCheck')
                                    if ~isnan(SessionData.RawEvents.Trial{1, q}.States.LickCheck(2))
                                        left_stim_time = SessionData.RawEvents.Trial{1, q}.States.LickCheck(2);
                                        lick_window = 5; 
                                    end
                                end
                                if ~isempty(licky)
                                    licky = licky(find(licky >= (left_stim_time))); %9
                                end
                                if ~isempty(licky)
                                   licky = licky(find(licky < left_stim_time + lick_window)); %10 changed from
                                end                
                                if ~isempty(licky)
                                    %Exclude trials with no licking
                                    disp(q);
                                    first_lick = licky(1);
                                    reward_time = first_lick;
                                    if outcomes(q) == 3
                                        water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                        water_assisted_ind = horzcat(water_assisted_ind, q);
                                        water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                        water_delivered_ind = horzcat(water_delivered_ind, q);
                                        left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                        left_ind = horzcat(left_ind, q);
                                    else
                                        water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                        left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                        left_ind = horzcat(left_ind, q);
                                        water_delivered_ind = horzcat(water_delivered_ind, q);
                                        water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                        water_earned_ind = horzcat(water_earned_ind, q);
                                    end
                                else
                                    disp(horzcat('DELETED FOR NO LICKING: ', num2str(q)));
                                end
                            else
                                if outcomes(q) == 3
                                    water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                    water_assisted_ind = horzcat(water_assisted_ind, q);
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                    left_ind = horzcat(left_ind, q);
                                else
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                    water_earned_ind = horzcat(water_earned_ind, q);
                                    left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                    left_ind = horzcat(left_ind, q);
                                end
                            end
                        end
                    else
                        no_water_delivered_ind = horzcat(no_water_delivered_ind, q);
                        %left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                        no_left_ind = horzcat(no_left_ind, q);
                    end
                else
                    no_water_delivered_ind = horzcat(no_water_delivered_ind, q);
                end

            else
                %If right or air
                reward_state = 'DeliverRight';
                if or(strcmp(alignment, 'left'), strcmp(alignment, 'right'))
                    reward_state = 'DeliverRight';
                elseif strcmp(alignment, 'both')
                    reward_state = 'DeliverRight';
                end

                if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)')))
                    reward_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)'));
                    if q ~= 1
                        if isfield(SessionData, 'TrialSettings')
                            if isfield(SessionData.TrialSettings(q - 1).GUI, 'RewardAmount')
                                reward = SessionData.TrialSettings(q - 1).GUI.RewardAmount;
                            elseif isfield(SessionData.TrialSettings(q - 1).GUI, 'RightAmount')
                                reward = SessionData.TrialSettings(q - 1).GUI.RightAmount;
                            else
                                reward = first_volume;
                            end
                        else
                            reward = 1;
                        end
                    else
                        reward = first_volume;
                    end
                    if and(~isnan(reward_time), or(reward == reward_amount, reward_amount == 1000))
                        if leave_out_no_licks
                            %Find trials with no licking after the US
                            licky = right_lick_event_times - trial_start;
                            right_stim_time = SessionData.RawEvents.Trial{1, q}.States.DeliverRight(1);
                            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'LickCheck')
                                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.LickCheck(2))
                                    right_stim_time = SessionData.RawEvents.Trial{1, q}.States.LickCheck(2);
                                    lick_window = 5; 
                                end
                            end
                            if ~isempty(licky)
                                licky = licky(find(licky >= (right_stim_time))); %9
                            end
                            if ~isempty(licky)
                               licky = licky(find(licky < right_stim_time + lick_window)); %10 changed from
                            end

                            if ~isempty(licky)
                                %Exclude trials with no licking
                                first_lick = licky(1);
                                reward_time = first_lick;
                                disp(q);
                                if outcomes(q) == 3
                                    water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                    water_assisted_ind = horzcat(water_assisted_ind, q);
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                    right_ind = horzcat(right_ind, q);
                                else
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                    water_earned_ind = horzcat(water_earned_ind, q);
                                    right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                    right_ind = horzcat(right_ind, q);
                                end
                            else
                                disp(horzcat('DELETED FOR NO LICKING: ', num2str(q)));
                            end
                        else
                            if outcomes(q) == 3
                                water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                water_assisted_ind = horzcat(water_assisted_ind, q);
                                water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                water_delivered_ind = horzcat(water_delivered_ind, q);
                                right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                right_ind = horzcat(right_ind, q);
                            else
                                water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                water_delivered_ind = horzcat(water_delivered_ind, q);
                                water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                water_earned_ind = horzcat(water_earned_ind, q);
                                right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                right_ind = horzcat(right_ind, q);
                            end
                        end
                    end
                else
                    no_water_delivered_ind = horzcat(no_water_delivered_ind, q);
                    no_right_ind = horzcat(no_right_ind, q);
                end
            
            end
            %Record sound delivered times
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'StimulusDeliver')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.StimulusDeliver(1))
                    stimulus_time = SessionData.RawEvents.Trial{1, q}.States.StimulusDeliver(1);
                    sound_delivered_times = horzcat(sound_delivered_times, stimulus_time + trial_start);% - session_start);
                    sound_delivered_ind = horzcat(sound_delivered_ind, q);% - session_start);
                    if trial_types(q) == 0 %water
                        water_sounds = [];
                        water_sound_ind = [];
                    elseif trial_types(q) == 1
                        air_sounds = [];
                        air_sound_ind = [];
                    end
                else
                    sound_delivered_times = horzcat(sound_delivered_times, 0); %Record trials with no sound as the start of the trial
                end
            end

            %Avoid time
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'DelayPre')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.DelayPre(1))
                    avoid_time = SessionData.RawEvents.Trial{1, q}.States.DelayPre(1);
                    avoid_times = horzcat(avoid_times, avoid_time + trial_start);% - session_start);
                    avoid_ind = horzcat(avoid_ind, q);
                else
                   not_avoided = horzcat(not_avoided, q);
                   total_not_avoided = horzcat(total_not_avoided, q);
                end
            end

            %Running initiation times
            if isfield(SessionData.RawEvents.Trial{1, q}.Events, 'Wire1High')
                if ~isempty(SessionData.RawEvents.Trial{1, q}.Events.Wire1High) %Might need to change to isfield
                    run_time = SessionData.RawEvents.Trial{1, q}.Events.Wire1High;

                    run_times = horzcat(run_times, run_time + trial_start);% - session_start);

                end
            end

            %Escape time
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'Escape')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.Escape(1))
                    escape_time = SessionData.RawEvents.Trial{1, q}.States.Escape(1);
                    escape_times = horzcat(escape_times, escape_time + trial_start);% - session_start);
                    escape_ind = horzcat(escape_ind, q);
                end
            end

            %Record laser times
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'Laser')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.Laser(1))
                    stimulus_time = SessionData.RawEvents.Trial{1, q}.States.Laser(1);
                    laser_times = horzcat(laser_times, stimulus_time + trial_start);% - session_start);
                    laser_ind = horzcat(laser_ind, q);
                end
            end

            %Record Air puff times
            if trial_types(q) == 1
                if isfield(SessionData.RawEvents.Trial{1, q}.States, 'DeliverPunishment')
                    if ~isnan(SessionData.RawEvents.Trial{1, q}.States.DeliverPunishment(1))
                        %Keeps track of the air puff length
                        if q ~= 1
                            if isfield(SessionData, 'TrialSettings')
                                punish = SessionData.TrialSettings(q - 1).GUI.RewardAmount;
                            else
                                punish = 1;
                            end
                        else
                            punish = first_volume;
                        end
                        punishment_time = SessionData.RawEvents.Trial{1, q}.States.DeliverPunishment(1);
                        punishment_end = SessionData.RawEvents.Trial{1, q}.States.DeliverPunishment(2);
                        if and(~isnan(punishment_time), or(punish == punish_amount, punish_amount == 1000))
                            %record the air puff if length is right
                            air_puff_times = horzcat(air_puff_times, punishment_time + trial_start);
                            air_puff_ind = horzcat(air_puff_ind, q);
                            air_puff_end_times = horzcat(air_puff_end_times, punishment_end + trial_start);
                            air_lengths = horzcat(air_lengths, punishment_end - punishment_time);
                            air_fail_ind = horzcat(air_fail_ind, q);
                        end
                    
                    else
                        no_air_delivered = horzcat(no_air_delivered, q);
                    end
                else
                    no_air_delivered = horzcat(no_air_delivered, q);
                end
            end
            if trial_types(q)
                punish_i = punish_i + 1;
                %air_trial_ind = horzcat(air_trial_ind, q);
            end
            left_lick_event_times = [];
            right_lick_event_times = [];
        end
    else
        trial_types = ones(1, length(trials_to_check)) .* 4;
    end

    %Find running initiation times. This is any running bought at least
    %1 second after the last running bought and continues for 3 or .
    run = -1;
    initiations = [];
    bought_length = 3;
    bought_lengths = [];
    for p = 1:size(run_times,2)
        next_run = run_times(1,p);
        if next_run - run > 1
            if bought_length < 3
                %If the bought length was too short, delete that run
                if ~isempty(initiations)
                    initiations = initiations(1:end-1);
                end
            else
                if or(~isempty(initiations), p == size(run_times,2))
                    bought_lengths = horzcat(bought_lengths, bought_length);
                end
            end
            %Add the next running initiation
            initiations = horzcat(initiations, next_run);
            %reset bought length to 1
            bought_length = 1;
        elseif p == size(run_times,2)
            bought_lengths = horzcat(bought_lengths, bought_length);
        elseif next_run - run < .1
            bought_length = bought_length + 1;
            %Next make it so at the start of the next iteration, if
            %bought length is less than 2, dlete previous bought.
        elseif next_run - run < 1
            disp('I dont know what to do here');

        end
        run = next_run;
    end

    %Cycle through start times. Segment from start_times(q) to q + 1.
    %include everything < 4 seconds as baseline.
    first_start = start_times(1);
    baseline_initiations = [];
    for q = 2:size(start_times, 2)
        next_start = start_times(q);
        %Find all running initiations in the current trial
        trial_segment = initiations(find(and(initiations > first_start, initiations < next_start)));
        %Put into time scale of 10 second trial
        trial_segment = trial_segment - first_start;
        %Find all places where running initiates during baseline
        trial_segment = trial_segment(find(trial_segment < 4)); %trial_segment(find(or(trial_segment < 4, trial_segment > 10)));
        %Convert back to timescale of imaging data
        trial_segment = trial_segment + first_start;
        %Keep track of the initiations
        baseline_initiations = horzcat(baseline_initiations, trial_segment);
        first_start = next_start;
    end
    if sound_on
        event_labels = {'Signal', 'CS', alignment};
    else
        event_labels = {'Signal', alignment};
    end
    
    if strcmp(alignment, 'events')
        [LeftPokes, MidPokes, RightPokes] = ExtractPokeTimes(SessionData);
        which_event = input('Left (l), middle (m), or right (r)? ', 's');
        if strcmp(which_event, 'l')
            event_times = {LeftPokes};
            event_labels = {'Signal', foods{1}};
        elseif strcmp(which_event, 'm')
            event_times = {MidPokes};
            event_labels = {'Signal', foods{2}};
        elseif strcmp(which_event, 'r')
            event_times = {RightPokes};
            event_labels = {'Signal', foods{3}};
        end
        event_trial_ind = 1:size(event_times{1}, 2);
    elseif strcmp(alignment, 'water all')
        event_times = {water_delivered_times};
        no_event_trial_ind = no_water_delivered_ind;
        event_trial_ind = water_delivered_ind;
    elseif strcmp(alignment, 'assist')
        event_times = {water_assisted_times};
        no_event_trial_ind = no_water_delivered_ind;
        event_trial_ind = water_assisted_ind;
    elseif strcmp(alignment, 'water')
        event_times = {water_earned_times};
        no_event_trial_ind = no_water_delivered_ind;
        event_trial_ind = water_earned_ind;
    elseif strcmp(alignment, 'air')
        event_times = {air_puff_times};
        event_end_times = air_puff_end_times;
        no_event_trial_ind = no_air_delivered;
        event_trial_ind = air_puff_ind;
    elseif strcmp(alignment, 'baseline run')
        event_times = {baseline_initiations};
    elseif strcmp(alignment, 'escape')
        event_times = {escape_times};
        event_trial_ind = escape_ind;
    elseif strcmp(alignment, 'avoid')
        event_times = {avoid_times};
        event_trial_ind = avoid_ind;
    elseif strcmp(alignment, 'laser')
        event_times = {laser_times};
        event_trial_ind = laser_ind;
    elseif strcmp(alignment, 'air sound')
        event_times = air_sounds;
        event_trial_ind = air_sound_ind;
    elseif strcmp(alignment, 'left')
        event_times = {};
        event_times{1} = left_deliveries;
        event_times{2} = right_deliveries;
        event_trial_ind = left_ind;
        %no_event_trial_ind = right_ind;
        no_event_trial_ind = setdiff(trials_to_check, left_ind);
    elseif strcmp(alignment, 'right')
        event_times = {};
        event_times{1} = right_deliveries;
        event_times{2} = left_deliveries;
        event_trial_ind = right_ind;
        %no_event_trial_ind = left_ind;
        no_event_trial_ind = setdiff(trials_to_check, right_ind);
    elseif strcmp(alignment, 'both')
        event_times = {};
        total_deliveries = horzcat(left_deliveries, right_deliveries);
        event_times{1} = sort(total_deliveries);
        event_times{2} = [];
        total_ind = horzcat(left_ind, right_ind);
        event_trial_ind = sort(total_ind);
        %no_event_trial_ind = left_ind;
        no_event_trial_ind = setdiff(trials_to_check, total_ind);
    end
    
    if sound_on
        event_times{2} = sound_delivered_times;
    end
    if trials_to_check(1) > 1
        %event_trial_ind = event_trial_ind - trials_to_check(1) + 1;
        %no_event_trial_ind = no_event_trial_ind - trials_to_check(1) + 1;
    end
end

function [gcamp, iso] = de_interleave(trigger1, trigger2, file, fiber_name, data, fiber, auto_flo, subtract_autofluorescence, signal_v)
%% de-interleave and subtract auto fluorescence
    %Use this if you alternate channels. Can handle 2 channels
    if ~and(trigger2, and(exist(horzcat(file, fiber_name, '/', 'signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'control.mat')) == 2))
        if trigger1 == 1
            gcamp = data(:, [1 fiber + 1]);
            avg_F = mean(gcamp(:,2));
            temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit
            close all
            %fitted is the line that forms the threshold to split channels
            adjust_threshold = 0;
            fitted = temp_fit(gcamp(:,1)) + adjust_threshold;
            for tries = 1:5
                figure
                plot(gcamp(:,1),gcamp(:,2)) %plot data and fit on same graph to make sure fit looks ok
                hold on
                fitted = fitted + adjust_threshold;

                plot(gcamp(:,1), fitted);
                title('Threshold to divide channels')
                is_good = input('Check whether the channel separation is good. Press enter if yes, n to adjust threshold, e to exit. ', 's');

                if strcmp(is_good, 'e')
                    close all
                    disp('User exited.')
                    return
                elseif strcmp(is_good, 'n')
                    new_thresh = input('Enter value to add to threshold line to adjust: (+ or -). e to exit ', 's');
                    if strcmp(new_thresh, 'e')
                        break
                    else
                        adjust_threshold = str2num(new_thresh);
                        hold off 
                        close all
                    end
                    continue
                else
                    break
                end
            end
            if tries == 5
                error('Are you sure your threshold is ok?');
            end
            close all;
            %plot(temp_fit);
            original = gcamp(:,2);
            %fitted = temp_fit(gcamp(:,1));
            hold off
            close all
            data_1 = [];
            data_2 = [];
            for x = 1:size(gcamp, 1)
                time = gcamp(x, 1);
                sig = gcamp(x, 2);
                %set fit as boundary threshold to separate channels
                thresh = fitted(x);
                if gcamp(x, 2) > thresh
                    data_1 = vertcat(data_1, [time, sig]);         
                else
                    data_2 = vertcat(data_2, [time, sig]);           
                end
            end
            figure;
            subplot(2,1,1)
            plot(data_1(:,1), data_1(:,2));
            title('Channel 1')
            subplot(2,1,2)
            plot(data_2(:,1), data_2(:,2));
            title('Channel 2')
            if signal_v == 1
                which_channel = input('Is the signal from channel 1 or 2? e to exit, answer either 1 or 2: ', 's');
            else
                which_channel = input('Is the control from channel 1 or 2? e to exit, answer either 1 or 2: ', 's');
            end
            

            if strcmp(which_channel, 'e')
                close all
                return
            end
            if trigger2 == 1
                channel_type = input('Control channel (c) or signal (s)? ','s');
            end
            close all;
            if str2num(which_channel) == 1
                if trigger2 == 0
                    gcamp = data_1;
                    iso = data_2;
                else
                    if strcmp(channel_type, 's')
                        gcamp = data_1;
                        iso = data_2;
                    elseif strcmp(channel_type, 'c')
                        gcamp = data_2;
                        iso = data_1;
                    end

                end
            elseif str2num(which_channel) == 2
                if trigger2 == 0
                    gcamp = data_2;
                    iso = data_1;
                else
                    if strcmp(channel_type, 's')
                        gcamp = data_2;
                        iso = data_1;
                    elseif strcmp(channel_type, 'c')
                        gcamp = data_1;
                        iso = data_2;
                    end
                end

            end
        end
    end

    if trigger2 == 1
        %Specify which channel you're looking at and save
        if and(exist(horzcat(file, fiber_name, '/', 'signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'control.mat')) == 2)
            load(horzcat(file, fiber_name, '/', 'control.mat'))
            load(horzcat(file, fiber_name, '/', 'signal.mat'))

        elseif strcmp(channel_type, 's')
            auto_flo_s = auto_flo;
            if ~exist(horzcat(file, fiber_name, '/', 'signal.mat'), 'var')   
                save(horzcat(file, fiber_name, '/', 'signal.mat'), 'gcamp', 'auto_flo_s')
            end
            disp('Add a control channel. Exited.');
            cd(file)
            return;

        elseif strcmp(channel_type, 'c')
            auto_flo_c = auto_flo;
            if ~exist(horzcat(file, fiber_name, '/', 'control.mat'), 'var')
                save(horzcat(file, fiber_name, '/', 'control.mat'), 'iso', 'auto_flo_c')
            end
            disp('Add a signal channel if you have not yet. Exited.');
            cd(file)
            return;
        end

        %If you have saved control and signal, load them as iso_trials and
        %trials
    elseif trigger1 == 1
        auto_flo_s = auto_flo;
        auto_flo_c = auto_flo;
        %gcamp = data(:, [1 fiber + 1]);
    else
        auto_flo_s = auto_flo;
        auto_flo_c = auto_flo;
        gcamp = data(:, [1 fiber + 1]);
    end
    
    %% Subtract autofluorescence
    if subtract_autofluorescence == 1
        gcamp(:,2) = gcamp(:,2) - auto_flo_s;
        if or(trigger1, trigger2)
            iso(:,2) = iso(:,2) - auto_flo_c;
        end
    else
        auto_flo = 0;
    end
end

function [SessionData, data, start_times, alt_start_times] = load_data(file, bpod_key_word, signal_key_word, analog_key_word, behavior_data, analog_input, record_all_starts) %#ok<STOUT>
    %Load bpod, photometry, and analog data
    
    %Retrieve the bpod behavior file
    tem = dir(horzcat(file, '*', bpod_key_word, '*.mat')); 
    if isempty(tem)
        tem = dir(horzcat(file, '*', bpod_key_word, '*.mat'));
        if isempty(tem)
            tem = dir(horzcat(file, '*behavior*.mat'));
        end
    end
    if isempty(tem)
        error('set signal_key_word to a phrase that appears in the bpod file')
    end
    %Load photometry data
    signal_file = dir(horzcat(file, '*', signal_key_word, '*.csv'));
    data = dlmread(horzcat(file, signal_file.name));
    start_file = dir(horzcat(file, '*', analog_key_word, '*.csv'));
    if ~isempty(start_file)
        photometry_start = dlmread(horzcat(file, start_file.name));
    else
        photometry_start = [0];
    end

    %Convert analog in file to seconds
    photometry_start_clone = photometry_start;
    photometry_start_clone = (photometry_start_clone - photometry_start(1,1))/1000;
    if photometry_start_clone(2) - photometry_start_clone(1) > 1
        photometry_start_clone(1) = [];
        photometry_start(1,:) = [];
        photometry_start_clone = photometry_start_clone - photometry_start_clone(1);
    end
    photometry_start_clone = photometry_start_clone';

    prev = 1;
    cur = 1;
    alt_start_times = [photometry_start_clone(1)];

    if behavior_data
        name = tem(1).name;
        folder = tem(1).folder;
        full_path = horzcat(folder, '/', name);
        disp(name);
        load(full_path);
        
    else
        data(:,1) = data(:,1);
    end

    for s = 2:size(photometry_start_clone,2)
        dif = photometry_start_clone(s) - photometry_start_clone(prev);
        if dif > 5
            alt_start_times = horzcat(alt_start_times, photometry_start_clone(s));
        end
        prev = prev + 1;
    end

    %% Use this to keep track of matlab events
    %set 0 point for imaging beginning. This coincides with trial 1.
    startRef = photometry_start(1);

    %convert time to relative time (s)
    start_times = alt_start_times;
    data(:,1) = ((data(:,1) - startRef)/1000);
    %data(:,1) = data(:,1) - data(1,1);
end
%{
function [gcamp, iso] = correct_bleaching(file, correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, trigger1, smoot)
    %% Correct bleaching- file: the folder/to save in, correct_for_bleaching: if 
    %% 1, fit an exponential, then subtract that, gcamp(n x fiber#) matrix with 
    %% time as first column, signal in others. iso: control channel data, 
    %% moving_avg: if 1, use a sliding window for correction. trigger1: 1 if 
    %% imaging using trigger 1. smoot: 1 if you want to smooth the data 
    %% (I don't recommend smoothing)
    
    first_tau = 10; %1
    second_tau = 12; %3
    
    gcamp_clone = gcamp;
    if correct_for_bleaching
        avg_F = mean(gcamp(:,2));
        plot(gcamp(:,1),gcamp(:,2))
        temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit for gcamp
        close all
        plot(gcamp(:,1),gcamp(:,2)) %plot data and fit on same graph to make sure fit looks ok
        hold on 
        plot(temp_fit);
        original = gcamp(:,2);
        fitted = temp_fit(gcamp(:,1));
        check = (original-fitted)./fitted;
        gcamp(:,2) = check;
        close all
        plot(gcamp(:,1),gcamp(:,2));
        close all
    elseif moving_avg
        roiDATA=gcamp(:, 2);
         %this should be equal to the sampling rate...in your case, if we aquired at 1khz, it sould be 1000/binsize (i think yours was 50ms)
        tau1 = round(frame_rate * first_tau);%0.75 ); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
        tau2 = round(frame_rate * second_tau);
        roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
        avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit

        for t = 1:size(roiDATA,1)
            roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
        end 
        dF_s= (roiDATA-roi_fzero)./roi_fzero;
        gcamp(:,2) = dF_s;
        figure;plot(gcamp(:,1), dF_s);
        title('Windowed average');
        %{
        if and(trigger1, trigger2)
            %Correct bleaching in the control channel
            roiDATA=iso(:, 2);
             %this should be equal to the sampling rate...in your case, if we aquired at 1khz, it sould be 1000/binsize (i think yours was 50ms)
            tau1 = round(frame_rate * 1);%0.75 ); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
            tau2 = round(frame_rate * 3);
            roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
            avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit

            for t = 1:size(roiDATA,1)
                roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
            end 
            dF_c= (roiDATA-roi_fzero)./roi_fzero;
            iso(:,2) = dF_c;
            figure;plot(iso(:,1), dF_c);
            title('Windowed average');
        end
        %}
    end

    times = gcamp(:,1);
    vals = gcamp(:,2);
    time_intervals = [];
    val_intervals = [];
    for x = 1:length(times)
        if x == 1
            prev_x = times(x);
            prev_v = vals(x);
        else
            time_intervals = horzcat(time_intervals, times(x) - prev_x);
            prev_x = times(x);
            val_intervals = horzcat(val_intervals, vals(x) - prev_v);
            prev_v = vals(x);
        end
    end
    plot(val_intervals);
    close all;
    keepers = find(val_intervals < 0.2);
    take_out = find(val_intervals > 0.2);
    times = gcamp_clone(:,1);
    vals = gcamp_clone(:,2);
    gcamp = horzcat(times(keepers), vals(keepers));
    plot(gcamp(:,1),gcamp(:,2));
    if correct_for_bleaching
        plot(gcamp(:,1),gcamp(:,2))
        temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit for gcamp
        close all
        plot(gcamp(:,1),gcamp(:,2)) %plot data and fit on same graph to make sure fit looks ok
        hold on 
        plot(temp_fit);
        original = gcamp(:,2);
        fitted = temp_fit(gcamp(:,1));
        check = (original-fitted)./fitted;
        gcamp(:,2) = check;
        close all
        plot(gcamp(:,1),gcamp(:,2));
        close all
        if trigger1
            temp_fit = fit(iso(:,1),iso(:,2),'exp2'); %find bi-exponential best fit for gcamp
            close all
            plot(iso(:,1),iso(:,2)) %plot data and fit on same graph to make sure fit looks ok
            hold on 
            plot(temp_fit);
            original = iso(:,2);
            fitted = temp_fit(iso(:,1));
            check = (original-fitted)./fitted;
            iso(:,2) = check;
            close all
            plot(iso(:,1),iso(:,2));
            close all
        end
    elseif moving_avg
        roiDATA=gcamp(:, 2);
        tau1 = round(frame_rate * first_tau); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
        tau2 = round(frame_rate * second_tau);
        roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
        avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit
        for t = 1:size(roiDATA,1)
            roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
        end 
        dF= (roiDATA-roi_fzero)./roi_fzero;
        gcamp(:,2) = dF;
        figure;plot(gcamp(:,1), gcamp(:,2));
        title('Windowed average gcamp');
        if trigger1
            roiDATA=iso(:, 2);
            tau1 = round(frame_rate * first_tau); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
            tau2 = round(frame_rate * second_tau);
            roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
            avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit
            for t = 1:size(roiDATA,1)
                roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
            end 
            dF= (roiDATA-roi_fzero)./roi_fzero;
            iso(:,2) = dF;
            figure;plot(iso(:,1), iso(:,2));
            title('Windowed average iso');
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %gcamp_clone = gcamp;
    if trigger1
        if moving_avg
            iso(1:tau1, :) = [];
            gcamp(1:tau1, :) = [];
        end
        if size(gcamp, 1) <= size(iso, 1)
            iso(1:((size(iso, 1) - size(gcamp, 1))), :) = [];
        end   
        gcamp_clone = gcamp;
        green = gcamp(:,2);
        violet = iso(:,2);
        close all
        figure
        plot(gcamp(:,1), green,'g')
        hold on
        plot(iso(:,1), violet,'m')
        %plot(water_delivered_times, zeros(length(water_delivered_times), 1), 'color', 'c', 'MarkerSize', 10, 'Marker', 'o')
        title('gCaMP and Iso plotted together with no fitting')
        legend('gCaMP', 'Iso');
        %if ~exist(file, 'dir')
        %    mkdir(file);
        %end
        savefig(horzcat(file, '_Entire_Session_Plot_no_fit'))
        hold off
        green_cor = green;
        violet_cor = violet;
    elseif trigger1
        iso_size = size(iso, 1);
        gcamp_size = size(gcamp, 1);

        for d = 1:iso_size
            if isnan(iso(d, 2))
                continue
            else
                iso(1:(d - 1), :) = [];
                gcamp(1:(d - 1), :) = [];
                break
            end
        end

        if size(gcamp, 1) <= size(iso, 1)
            difff = iso_size - gcamp_size;
            if iso((difff + 1), 1) < gcamp(1,1)
                iso(1:((difff)), :) = [];
            else
                iso((end - ((difff)) + 1):end, :) = [];
            end
        elseif size(gcamp, 1) > size(iso, 1)
            first_g = gcamp(1,1);
            first_i = iso(1,1);
            if first_g > first_i
                s = 1;
                while first_g > first_i
                    s = s + 1;
                    first_i = iso(s,1);
                end
                iso(1:s, :) = [];
            elseif first_g < first_i
                s = 1;
                while first_g < first_i
                    s = s + 1;
                    first_g = gcamp(s,1);
                end
                gcamp(1:s, :) = [];
            end
            iso_size = size(iso, 1);
            gcamp_size = size(gcamp, 1);
            if iso_size < gcamp_size
                gcamp = gcamp(1:iso_size, :);
            else
                iso = iso(1:gcamp_size, :);
            end
        end 
        green = gcamp(:,2);
        violet = iso(:,2);
        green_cor = green;
        violet_cor = violet;            
        [coeff] = regress(green_cor,[ones(length(violet_cor),1) violet_cor]);
        violet_cor = [ones(length(violet_cor),1) violet_cor]*coeff;
        new_iso = horzcat(iso(:,1), violet_cor);
        iso = new_iso;
    end

    %Smooth the separated channel
    if smoot
        if correct_for_bleaching
            gcamp(:,2) = smooth(gcamp(:,2), 3);
        else
            gcamp(:,2) = smooth(gcamp(:,2), 5);
        end
    end
    close all
    
end
%}

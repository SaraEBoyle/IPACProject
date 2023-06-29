function [first_valve_volume, left_recorded_options, right_recorded_options] = find_first_volume(SessionData)
    %% compares valve times to volume labels to confirm that labels were assigned correctly
    % Returns volume of first trial because I designed something stupidly
    disp('Select bpod file');

    
    % get trial volumes. 1 is right, 0 is left
    if ~isfield(SessionData, 'spout')
    elseif SessionData.spout == 1
        trial_types = ones(1, SessionData.nTrials);
        first_valve_volume = SessionData.FirstLeftAmount;
        return
    elseif SessionData.spout == 2
        trial_types = zeros(1, SessionData.nTrials);
        first_valve_volume = SessionData.FirstRightAmount;
        return
    elseif SessionData.spout == 3
        temp = ones(1, SessionData.nTrials);
        evens = mod(1:SessionData.nTrials, 2);
        trial_types = temp .* ~evens;
    end
    if ~isfield(SessionData, 'TrialTypes')
        if or(SessionData.spout == 2, SessionData.spout == 3)
        if isfield(SessionData, 'FirstRightAmount')
            first_valve_volume = SessionData.FirstRightAmount;
            return
        else
            first_valve_volume = 2;
            return;
        end
        else
            if isfield(SessionData, 'FirstLeftAmount')
            first_valve_volume = SessionData.FirstLeftAmount;
            return
        else
            first_valve_volume = 2;
            return;
        end
        end
    else
        trial_types = SessionData.TrialTypes;
    end
    
    valve_times = [];
    recorded_volumes = [];
    left_recorded_volumes = [];
    right_recorded_volumes = [];
    %valve = SessionData.RawEvents.Trial{1, x}.States.DeliverRight;
    %valve = SessionData.RawEvents.Trial{1, x}.States.DeliverLeft;
    %left = SessionData.TrialSettings(x).GUI.LeftAmount;
    %right = SessionData.TrialSettings(x).GUI.RightAmount;
    
    for x = 1:length(trial_types)
        if x == 1
            if trial_types(x) == 1
                % Right trial
                first_valve_time = SessionData.RawEvents.Trial{1, x}.States.DeliverRight;
                first_valve_time = first_valve_time(2) - first_valve_time(1);
                %recorded_volume = SessionData.TrialSettings(x - 1).GUI.RightAmount;
            else
                % Left trial
                first_valve_time = SessionData.RawEvents.Trial{1, x}.States.DeliverLeft;
                first_valve_time = first_valve_time(2) - first_valve_time(1);
                %recorded_volume = SessionData.TrialSettings(x - 1).GUI.LeftAmount;
            end
            continue
        end
        if trial_types(x) == 1
            % Right trial
            valve_time = SessionData.RawEvents.Trial{1, x}.States.DeliverRight;
            valve_time = valve_time(2) - valve_time(1);
            right_recorded_volume = SessionData.TrialSettings(x - 1).GUI.RightAmount;
            recorded_volume = SessionData.TrialSettings(x - 1).GUI.RightAmount;
            right_recorded_volumes = horzcat(right_recorded_volumes, right_recorded_volume);
        else
            % Left trial
            valve_time = SessionData.RawEvents.Trial{1, x}.States.DeliverLeft;
            valve_time = valve_time(2) - valve_time(1);
            left_recorded_volume = SessionData.TrialSettings(x - 1).GUI.LeftAmount;
            recorded_volume = SessionData.TrialSettings(x - 1).GUI.LeftAmount;
            left_recorded_volumes = horzcat(left_recorded_volumes, left_recorded_volume);
        end
        valve_times = horzcat(valve_times, valve_time);
        recorded_volumes = horzcat(recorded_volumes, recorded_volume);
    end
    
    recorded_options = unique(recorded_volumes);
    right_recorded_options = unique(right_recorded_volumes);
    left_recorded_options = unique(left_recorded_volumes);
    converted_valves = zeros(1,length(valve_times));
    first_ind = 0;
    later_trial_types = trial_types(2:end);
    right_valve_times = later_trial_types .* valve_times;
    left_valve_times = ~later_trial_types .* valve_times;
    right_valve_options = unique(right_valve_times);
    left_valve_options = unique(left_valve_times);
    if right_valve_options(1) == 0
        right_valve_options(1) = [];
    end
    if left_valve_options(1) == 0
        left_valve_options(1) = [];
    end
    
    %% Sometimes bpod records replicate valve times. remove those
    to_delete = [];
    if size(left_valve_options, 2) > 1
        first_time = left_valve_options(1);
        
        for q = 2:size(left_valve_options, 2)
            cur = left_valve_options(q);
            if abs(cur - first_time) < .001
                to_delete = horzcat(to_delete, q);
            end
            first_time = cur;
        end
    end
    left_valve_options(to_delete) = [];
    
    to_delete = [];
    if size(right_valve_options, 2) > 1
        first_time = right_valve_options(1);
        
        for q = 2:size(right_valve_options, 2)
            cur = right_valve_options(q);
            if abs(cur - first_time) < .001
                to_delete = horzcat(to_delete, q);
            end
            first_time = cur;
        end
    end
    right_valve_options(to_delete) = [];
    
    
    for y = 1:length(left_valve_options)
        left_inds = find(valve_times >= left_valve_options(y) - .001);
        right_inds = find(valve_times <= left_valve_options(y) + .001);
        inds = intersect(left_inds, right_inds);
        if and(first_valve_time >= left_valve_options(y) - .001, first_valve_time <= left_valve_options(y) + .001)
            first_ind = 1;
        end
        converted_valves(inds) = left_recorded_options(y);
        disp(horzcat(num2str(left_recorded_options(y)), ' uL left trials: ', num2str(length(inds) + first_ind)));
        if first_ind 
            if and(first_valve_time >= left_valve_options(y) - .001, first_valve_time <= left_valve_options(y) + .001)
                first_valve_volume = left_recorded_options(y);
            end
        end
        first_ind = 0;

    end
    for y = 1:length(right_valve_options)
        left_inds = find(valve_times >= right_valve_options(y) - .001);
        right_inds = find(valve_times <= right_valve_options(y) + .001);
        inds = intersect(left_inds, right_inds);
        if and(first_valve_time >= right_valve_options(y) - .001, first_valve_time <= right_valve_options(y) + .001)
            first_ind = 1;
        end
        disp(horzcat(num2str(right_recorded_options(y)), ' uL right trials: ', num2str(length(inds) + first_ind)));
        converted_valves(inds) = right_recorded_options(y);
        if first_ind
            if and(first_valve_time >= right_valve_options(y) - .001, first_valve_time <= right_valve_options(y) + .001)
                first_valve_volume = right_recorded_options(y);
            end
        end
        first_ind = 0;
    end
    diffs = sum(converted_valves~= recorded_volumes);
    diffs = diffs/length(trial_types) * 100;
    disp(horzcat('Volume assignments agree with valve times ', num2str(100 - diffs), '% of the time' ));
    %if diffs == 0
    %    disp('Volumes were assigned perfectly. Do not worry');
    %else
    %    error('There is something wrong with volume assignment.');
    %end
    
    
    
end
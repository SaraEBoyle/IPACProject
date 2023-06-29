function [] = confirm_volumes()
    %% compares valve times to volume labels to confirm that labels were assigned correctly
    % Choose bpod file to look at.
    disp('Select bpod file');
    
    [file, path] = uigetfile('.mat');
    load(horzcat(path, file));
    
    % get trial volumes. 1 is right, 0 is left
    trial_types = SessionData.TrialTypes;
    valve_times = [];
    recorded_volumes = [];
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
            recorded_volume = SessionData.TrialSettings(x - 1).GUI.RightAmount;
        else
            % Left trial
            valve_time = SessionData.RawEvents.Trial{1, x}.States.DeliverLeft;
            valve_time = valve_time(2) - valve_time(1);
            recorded_volume = SessionData.TrialSettings(x - 1).GUI.LeftAmount;
        end
        valve_times = horzcat(valve_times, valve_time);
        recorded_volumes = horzcat(recorded_volumes, recorded_volume);
    end
    valve_options = unique(valve_times);
    recorded_options = unique(recorded_volumes);
    converted_valves = zeros(1,length(valve_times));
    first_ind = 0;
    for y = 1:length(valve_options)
        inds = find(valve_times == valve_options(y));
        if first_valve_time == valve_options(y)
            first_ind = 1;
        end
        disp(horzcat(num2str(recorded_options(y)), ' uL trials: ', num2str(length(inds) + first_ind)));
        converted_valves(inds) = recorded_options(y);
        if first_valve_time == valve_options(y)
            first_valve_volume = recorded_options(y);
        end
        first_ind = 0;
    end
    diffs = sum(converted_valves~= recorded_volumes);
    diffs = diffs/length(trial_types) * 100;
    disp(horzcat('Volume assignments agree with valve times ', num2str(100 - diffs), '% of the time' ));
    if diffs == 0
        disp('Volumes were assigned perfectly. Do not worry');
    else
        error('There is something wrong with volume assignment.');
    end
    volumes = {converted_valves, recorded_volumes};
    
end
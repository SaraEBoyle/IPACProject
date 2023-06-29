%%HOW TO USE THIS SCRIPT
%Plots the signal from the entire session. If there's signal from multiple
%channels, it plots the fitted control channel as well. Also does bleaching
%correction
function  plot_signal(S)
    %% Set up Variables
    %Select which type of trial to check: Don't change
    correct_motion = 0;
    motion_whole_session = 1;
    handles = guidata(S.fh);
    handles.GUI_35.plotted = 1;
    file = pwd;
    file = horzcat(file, '/');
    fit_control = 1;
    %This includes GUI_35.new_graph_who_dis and plotted
    violet = S.c(1).Value;
    green = S.c(2).Value;
    red = S.c(3).Value;
    plot_control = handles.GUI_35.plot_control;
    data = handles.GUI_35.data;
    analog = handles.GUI_35.analog;
    no_analog = 0;
    subtract_autofluorescence = S.plot_background.Value;
    sync = 1; %When 1, synchronize the 2 channels to get a common timescale
    frame_rate = str2double(S.frame_rate.String);
    fiber = str2double(S.fiber.String); %Choose which fiber to look at. Can add more as needed. This just selects the column of the signal csv to look at.
    area_ind = S.area.DropDown.Value;
    fiber_name = S.area.DropDown.String{area_ind};
    mov_avg = handles.GUI_35.mov_avg;
    bi_ex = handles.GUI_35.bi_ex;
    
    if and(green,violet)
            trigger1 = 1; %Change to 1 if your signal alternates through 2 channels
            trigger2 = 0;
        elseif and(green, red)
            trigger1 = 1; %Change to 1 if your signal alternates through 2 channels
            trigger2 = 1;
        else
            trigger1 = 0; %Change to 1 if your signal alternates through 2 channels
            trigger2 = 0;
    end
    
    %If you've already analyzed this data, load the corrected data
    if and(exist(horzcat(file, fiber_name, '/', 'raw_signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'raw_control.mat')) == 2)
            load(horzcat(file, fiber_name, '/', 'raw_control.mat'))
            load(horzcat(file, fiber_name, '/', 'raw_signal.mat'))
    %Otherwise, compute your dF/F
    else
        if and(mov_avg, bi_ex)
            disp('You chose two bleaching correction methods. Do not do that.')
            bi_ex = 0;
        end

        %data = signal;
        photometry_start = analog;
        photometry_start_clone = photometry_start;
        photometry_start_clone = (photometry_start_clone - photometry_start(1,1))/1000;
        photometry_start_clone = photometry_start_clone';

        for x = 1:length(photometry_start_clone)
            if photometry_start_clone(x) > 0
                photometry_start_ind = x;
                break
            end
        end

        prev = 1;
        cur = 1;
        if no_analog
            alt_start_times = [1 10 20 30 40];
            start_times = alt_start_times;
        else
            alt_start_times = [photometry_start_clone(1)];
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
        if no_analog
            startRef = 0;
        else
            startRef = photometry_start(1); %- start_times;
        end
        %convert time to relative time (s)
        if no_analog
            data(:,1) = ((data(:,1) - data(1,1))/1000);
        else
            data(:,1) = ((data(:,1) - startRef)/1000); 
        end
        %% Cut data from beginning and/or end
        %Plot pre-bleach corrected data.
        if subtract_autofluorescence
            one = figure;
            plot(data(:,1), data(:,fiber + 1));
            two = figure;
            plot(data(1:800,1), data(1:800,fiber + 1));
            title('First 800 seconds. Choose autofluorescence window. Click at beginning, then at end')
            [x_in,y_in]=ginput(2);
            disp(horzcat('Window is from ', num2str(x_in(1)), ' to ', num2str(x_in(2))))
            %disp(y_in)

            is_good = input('Did you choose? Press enter to continue, e to exit: ', 's');
            if is_good == 'e'
                disp('User exited. Closed program.')
                close all
                return
            end
            close(one);
            close(two);
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
        three = figure;
        plot(data(:,1), data(:,fiber + 1));

        is_good = input('Does this look right? Press enter to continue without cutting any data, c to cut, e to exit: ', 's');
        if is_good == 'e'
            disp('User exited. Closed program.')
            return
        elseif is_good == 'c'
            disp('Select interval to analyze. ');
            [x_in,y_in] = ginput(2);
            %Cut anything funky from the end
            data = data(find(data(:,1) > x_in(1)), :);
            data = data(find(data(:,1) < x_in(2)), :);    
        end
        close(three)

        %% de-interleave
        %For trigger1 or trigger2 mode photometry
        if or(and(violet, green), and(red, green))
            gcamp = data(:, [1 fiber + 1]);
            avg_F = mean(gcamp(:,2));
            temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit
            %close all
            %fitted is the line that forms the threshold to split channels
            adjust_threshold = 0;
            fitted = temp_fit(gcamp(:,1)) + adjust_threshold;
            for tries = 1:5
                four = figure;
                plot(gcamp(:,1),gcamp(:,2)) %plot data and fit on same graph to make sure fit looks ok
                hold on
                fitted = fitted + adjust_threshold;

                plot(gcamp(:,1), fitted);
                title('Threshold to divide channels')
                is_good = input('Check whether the channel separation is good. Press enter if yes, n to adjust threshold, e to exit. ', 's');

                if strcmp(is_good, 'e')
                    %close all
                    disp('User exited.')
                    hold off
                    close(four);
                    return
                elseif strcmp(is_good, 'n')
                    new_thresh = input('Enter value to add to threshold line to adjust: (+ or -). e to exit ', 's');
                    if strcmp(new_thresh, 'e')
                        break
                    else
                        adjust_threshold = str2num(new_thresh);
                        hold off 
                        close(four);
                    end
                    continue
                else
                    break
                end
            end
            if tries == 5
                error('Are you sure your threshold is ok?');
            end

            %plot(temp_fit);
            original = gcamp(:,2);
            %fitted = temp_fit(gcamp(:,1));
            hold off
            close(four);
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
            sixteen = figure;
            channel_1 = subplot(2,1,1);
            plot(data_1(:,1), data_1(:,2));
            title('Channel 1')
            channel_2 = subplot(2,1,2);
            plot(data_2(:,1), data_2(:,2));
            title('Channel 2')
            which_channel = input('Is the signal from channel 1 or 2? e to exit, answer either 1 or 2: ', 's');

            if strcmp(which_channel, 'e')
                close(sixteen)
                return
            end
            if red == 1
                channel_type = input('Control channel (c) or signal (s)? ','s');
            end
            close(sixteen);
            if str2num(which_channel) == 1
                if red == 0
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
                if red == 0
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
        elseif red == 1
            channel_type = 's';

        end

            if and(trigger2, plot_control)
                %Specify which channel you're looking at and save
                if and(exist(horzcat(file, fiber_name, '/', 'raw_signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'raw_control.mat')) == 2)
                    load(horzcat(file, fiber_name, '/', 'raw_control.mat'))
                    load(horzcat(file, fiber_name, '/', 'raw_signal.mat'))
                elseif strcmp(channel_type, 's')
                    auto_flo_s = auto_flo;
                    if ~exist(horzcat(file, fiber_name, '/', 'raw_signal.mat'), 'var')   
                        save(horzcat(file, fiber_name, '/', 'raw_signal.mat'), 'gcamp', 'auto_flo_s')
                    end
                    disp('Add a control channel. Exited.');
                    to_check = input('Which fiber has your control channel? e for exit','s');
                    if strcmp(to_check, 'e')
                        disp('user exited');
                        return
                    end
                    S.fiber.String = to_check; 
                    plot_signal(S);
                    cd(file)        
                    return;
                elseif strcmp(channel_type, 'c')
                    auto_flo_c = auto_flo;
                    if ~exist(horzcat(file, fiber_name, '/', 'raw_control.mat'), 'var')
                        save(horzcat(file, fiber_name, '/', 'raw_control.mat'), 'iso', 'auto_flo_c')
                    end
                    if and(exist(horzcat(file, fiber_name, '/', 'raw_signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'raw_control.mat')) == 2)
                        cd(file)
                    else
                        disp('Add a signal channel if you have not yet. Exited.');
                        to_check = input('Which fiber has your signal channel?','s');
                        if strcmp(to_check, 'e')
                            disp('user exited');
                            return
                        end
                        S.fiber.String = to_check;
                        plot_signal(S);
                        cd(file)        
                        return;
                    end

                end
                %If you have saved control and signal, load them as iso_trials and
                %trials
            elseif and(violet, green)
                auto_flo_s = auto_flo;
                auto_flo_c = auto_flo;
            else
                auto_flo_s = auto_flo;
                auto_flo_c = auto_flo;
                gcamp = data(:, [1 fiber + 1]);
            end

            if and(exist(horzcat(file, fiber_name, '/', 'raw_signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'raw_control.mat')) == 2)
                load(horzcat(file, fiber_name, '/', 'raw_control.mat'))
                load(horzcat(file, fiber_name, '/', 'raw_signal.mat'))
            end
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

    gcamp_clone = gcamp;
    %% correct for bleaching
    if bi_ex
        %five is pre bleach correction
        separate = figure;
        avg_F = mean(gcamp(:,2));
        plot(gcamp(:,1),gcamp(:,2))
        temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit for gcamp
        set(separate, 'visible', 'off');
        %close(pre_bleach_correct);
        %six is with fit pre bleach correction
        bleach_correct_w_fit = figure;
        plot(gcamp(:,1),gcamp(:,2)) %plot data and fit on same graph to make sure fit looks ok
        hold on 
        plot(temp_fit);
        original = gcamp(:,2);
        fitted = temp_fit(gcamp(:,1));
        check = (original-fitted)./fitted;
        gcamp(:,2) = check;
        close(bleach_correct_w_fit);
        %seven is post correction
        seven = figure;
        plot(gcamp(:,1),gcamp(:,2));
        close(seven)
    elseif mov_avg
        roiDATA=gcamp(:, 2);
         %this should be equal to the sampling rate...in your case, if we aquired at 1khz, it sould be 1000/binsize (i think yours was 50ms)
        tau1 = round(frame_rate * 1);%0.75 ); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
        tau2 = round(frame_rate * 3);
        roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
        avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit

        for t = 1:size(roiDATA,1)
            roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
        end 
        dF_s= (roiDATA-roi_fzero)./roi_fzero;
        gcamp(:,2) = dF_s;
        eight = figure;
        plot(gcamp(:,1), dF_s);
        title('Windowed average');
        close(eight)
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
    keepers = find(val_intervals < 0.2);
    take_out = find(val_intervals > 0.2);
    times = gcamp_clone(:,1);
    vals = gcamp_clone(:,2);
    gcamp = horzcat(times(keepers), vals(keepers));
    
    %Plot the data pre correction
    pre_bleach_correct = figure;
    hold on 
    plot(gcamp(:,1),gcamp(:,2), 'color', 'g') %plot data and fit on same graph to make sure fit looks ok
    %plot(temp_fit);
    if(plot_control)
        plot(iso(:,1),iso(:,2), 'color', 'm');
        legend('Signal', 'Control');
    else
        legend('Signal');
    end
    set(pre_bleach_correct, 'visible', 'off');
    
    if bi_ex
        temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit for gcamp
        
        
        original = gcamp(:,2);
        fitted = temp_fit(gcamp(:,1));
        check = (original-fitted)./fitted;
        gcamp(:,2) = check;
        
        
        if and(and(violet, green), plot_control)
            temp_fit = fit(iso(:,1),iso(:,2),'exp2'); %find bi-exponential best fit for gcamp
            %close all
            twelve = figure;
            plot(iso(:,1),iso(:,2)) %plot data and fit on same graph to make sure fit looks ok
            hold on 
            plot(temp_fit);
            original = iso(:,2);
            fitted = temp_fit(iso(:,1));
            check = (original-fitted)./fitted;
            iso(:,2) = check;
            %close(post_bleach_correct_w_control)
            set(twelve, 'visible','off');
            post_bleach_correct = figure;
            plot(gcamp(:,1),gcamp(:,2), 'color', 'g');
            hold on
            plot(iso(:,1),iso(:,2), 'color', 'm');
            legend('Signal', 'Control');
            set(post_bleach_correct, 'visible','off');
        else %For only green
            post_bleach_correct = figure;
            plot(gcamp(:,1),gcamp(:,2), 'color', 'g');
            legend('Signal')
            set(post_bleach_correct, 'visible', 'off');
        end
    elseif mov_avg
        roiDATA=gcamp(:, 2);
        tau1 = round(frame_rate *8); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
        tau2 = round(frame_rate * 1);
        roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
        avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit
        for t = 1:size(roiDATA,1)
            roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
        end 
        dF= (roiDATA-roi_fzero)./roi_fzero;
        gcamp(:,2) = dF;
        
        if trigger1
            roiDATA=iso(:, 2);
            tau1 = round(frame_rate *8); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
            tau2 = round(frame_rate * 1);
            roi_fzero = zeros(size(roiDATA)); %this should be the size of your data matrix    
            avg_mat = tsmovavg(roiDATA(:,:),'s',tau1,1); %this filters the data a bit
            for t = 1:size(roiDATA,1)
                roi_fzero(t,:) = min(avg_mat(max(1,t-tau2):t,:));
            end 
            dF= (roiDATA-roi_fzero)./roi_fzero;
            iso(:,2) = dF;
            %figure;plot(iso(:,1), iso(:,2));
            %title('Windowed average iso');
        end
        post_bleach_correct = figure;
        hold on
        plot(gcamp(:,1), gcamp(:,2));
        title('Sliding Average dF/F');
        if plot_control
            plot(iso(:,1), iso(:,2));
        end
        set(post_bleach_correct, 'visible', 'off');
    end

    gcamp_clone = gcamp;
    if and(plot_control, trigger1) %%%%%%%%%%%%IMPORTANT but doesn't fit
        if mov_avg
            iso(1:tau1, :) = [];
            gcamp(1:tau1, :) = [];
        end
        
        if and(size(gcamp, 1) ~= size(iso, 1), fit_control)
            if gcamp(1,1) > iso(1,1)
                %Cut the iso lower than gcamp11
                to_cut = iso(:, 1) < gcamp(1,1);
                iso(to_cut, :) = [];
            else
                %Cut the gcamp lower than iso11
                to_cut = gcamp(:, 1) < iso(1,1);
                gcamp(to_cut, :) = [];
            end
        end
        if and(size(gcamp, 1) <= size(iso, 1), fit_control)
            %if iso(1,1) > gcamp(1,1)
            %    iso(1:((size(iso, 1) - size(gcamp, 1))), :) = [];
            %else
                iso(end + 1 - ((size(iso, 1) - size(gcamp, 1))):end, :) = [];
           % end
        end   
        gcamp_clone = gcamp;
        green = gcamp(:,2);
        violet = iso(:,2);
        %close all
        %ax = axes;
        post_bleach_correct_w_control = figure;
        plot(gcamp(:,1), green,'g')
        hold on
        if plot_control
            plot(iso(:,1), violet,'m')
            legend('Signal', 'Control');
        else
            legend('Signal');
        end
        %plot(ax, water_delivered_times, zeros(length(water_delivered_times), 1), 'color', 'c', 'MarkerSize', 10, 'Marker', 'o')
        title('Raw dF/F')
        
        savefig('Entire_Session_Plot_no_fit')
        
        %close(pre_bleach_correct_w_control);
        set(post_bleach_correct_w_control, 'visible','off');
        hold off
        green_cor = green;
        violet_cor = violet;

        if and(plot_control,fit_control) %%%%%%%%%%%%This is the fitting
            while(length(green_cor) ~= length(violet_cor))
                %If signal has more values than control
                if length(green_cor) > length(violet_cor)
                    %If signal channel came first
                    if iso(1,1) < gcamp(1,1)
                        %Chop off first signal value
                        green_cor = green_cor(2:end);
                    %If control channel came first
                    else
                        %Chop off last signal value
                        green_cor = green_cor(1:end - 1);
                    end
                %If control has more values than signal
                else
                    %If signal channel came first
                    if iso(1,1) < gcamp(1,1)
                        %Chop off last control value
                        violet_cor = violet_cor(1:end - 1);
                    %If control channel came first
                    else
                        %Chop off first control value
                        violet_cor = violet_cor(2:end);
                    end
                end
            end
            [coeff] = regress(green_cor,[ones(length(violet_cor),1) violet_cor]);
            violet_cor = [ones(length(violet_cor),1) violet_cor]*coeff;
            violet_cor = violet_cor + 1;
            green_cor = green_cor + 1;

            %Plot fitted graph
            post_bleach_correct_w_control = figure;
            plot(gcamp(:,1), green,'g')
            hold on
            if plot_control
                plot(iso(:,1), violet_cor - 1,'m')
                legend('Signal', 'Control');
            else
                legend('Signal');
            end
            %plot(water_delivered_times, zeros(length(water_delivered_times), 1), 'color', 'c', 'MarkerSize', 10, 'Marker', 'o')
            title('dF/F with fitting')
            savefig('Entire_Session_Plot_with_fit')
            hold off
            set(post_bleach_correct_w_control, 'visible', 'off');
        end
        if correct_motion
            if sync == 1
                %signal = green_cor - violet_cor;
                ts1 = timeseries(green_cor,gcamp(:,1));
                % Associate another frequency with amplitude
                ts2 = timeseries(violet_cor,iso(:,1));
                % Synchronize the frequency
                [ts3, ts4] = synchronize(ts1,ts2,'Union');
                % Plot the difference
                %figure
                %hold on
                %plot(ts3.Time, ts3.Data-ts4.Data);
                %title("Post correction, with syncing");
                signal = ts3.Data-ts4.Data;
                %new_gcamp = horzcat(ts3.Time, smooth(signal, 3));
                %new_iso = horzcat(ts3.Time, smooth(ts4.Data - 1, 3));
                new_gcamp = horzcat(ts3.Time, signal);
                new_iso = horzcat(ts3.Time, ts4.Data - 1);
                gcamp = new_gcamp;
                iso = new_iso;

            else
                new_iso = horzcat(iso(:,1), violet_cor - 1);
                iso = new_iso;
                signal = movmean(green_cor, 2) - movmean(violet_cor, 2);
                gcamp(:,2) = signal;
                %gcamp(:,2) = movmean(signal,2);
            end
            seventeen = figure;
            plot(gcamp_clone(:,1), green_cor - 1,'g')
            hold on
            plot(gcamp(:,1), gcamp(:,2),'k')
            title('After correction')
            legend('Pre','Post')
            savefig('post_correction_entire.fig')
            close(seventeen);
        end
    elseif and(plot_control, trigger1)
        %I think this is if you want to fit the control channel? Seems
        %useless
        %{
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
        %}
    end

    %stimulus_times_temp = stimulus_times;
    if correct_motion
        %Run motion corrector to plot pre corrected trials
        [synced_trials] = motion_corrector(0, check_trials, OG_trials, ...
        iso_trials, ones(length(trial_types)), types_to_check, motion_whole_session, sync, ...
        save_corrected_trials, trials_to_check, stim_align, stim_time, ...
        stimulus_times);
        %Run it to plot all corrected trials and save output as check
        %trials
        [corrected_trials] = motion_corrector(1, check_trials, OG_trials, ...
        iso_trials, ones(length(trial_types)), types_to_check, motion_whole_session, sync, ...
        save_corrected_trials, trials_to_check, stim_align, stim_time, ...
        stimulus_times);
        if ~isempty(corrected_trials)
            trials = corrected_trials;
        end
    end   
    if plot_control
        handles.GUI_35.pre_correction = pre_bleach_correct;
        handles.GUI_35.post_correction = post_bleach_correct_w_control;
        handles.GUI_35.control_dfF = iso;
        handles.GUI_35.signal_dfF = gcamp;
        save(horzcat(file, fiber_name, '/', 'dfF_signal.mat'), 'gcamp')
        save(horzcat(file, fiber_name, '/', 'dfF_control.mat'), 'iso')
    else
        handles.GUI_35.pre_correction = pre_bleach_correct;
        handles.GUI_35.post_correction = post_bleach_correct;
        handles.GUI_35.signal_dfF = gcamp;
        save(horzcat(file, fiber_name, '/', 'dfF_signal.mat'), 'gcamp')
    end
    %handles.GUI_35.post_correction = post_bleach_correct; % initialized in GUI main code 
    
    guidata(S.fh, handles)

end



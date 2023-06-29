function [gcamp, iso, baseline_gcamp, baseline_iso] = correct_bleaching(file, correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, trigger1, smoot, liquid_delivered_times)
    %% Correct bleaching- file: the folder/to save in, correct_for_bleaching: if 
    %% 1, fit an exponential, then subtract that, gcamp(n x fiber#) matrix with 
    %% time as first column, signal in others. iso: control channel data, 
    %% moving_avg: if 1, use a sliding window for correction. trigger1: 1 if 
    %% imaging using trigger 1. smoot: 1 if you want to smooth the data 
    %% (I don't recommend smoothing)
    
    gcamp_clone = gcamp;
    
    baseline_inds = [];
    baseline_dict = {};
    for e = 1:size(liquid_delivered_times, 2)
        time_scale = gcamp(:, 1);
        event = liquid_delivered_times(e);
        baseline_inds_left = find(time_scale >= (event - 13));
        baseline_inds_right = find(time_scale <= (event - 0));
        baseline_inds = vertcat(baseline_inds, intersect(baseline_inds_left, baseline_inds_right));
        baseline_dict{e} = [intersect(baseline_inds_left, baseline_inds_right)];
    end
    
    gcamp_minus_events = horzcat(gcamp(baseline_inds, 1), gcamp(baseline_inds, 2));
    iso_minus_events = horzcat(iso(baseline_inds, 1), iso(baseline_inds, 2));
    
    
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
        tau1 = round(frame_rate * 1);%0.75 ); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
        tau2 = round(frame_rate * 3);
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
        temp_fit = fit(gcamp_minus_events(:,1),gcamp_minus_events(:,2),'exp2'); %find bi-exponential best fit for gcamp
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
            temp_fit = fit(iso_minus_events(:,1),iso_minus_events(:,2),'exp2'); %find bi-exponential best fit for gcamp
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
        tau1 = round(frame_rate *8); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
        tau2 = round(frame_rate * 1);
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
            tau1 = round(frame_rate *8); %these two parameters you can change (refer to Helmchen Delta F image i showed you)
            tau2 = round(frame_rate * 1);
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
        if ~exist(file, 'dir')
            mkdir(file);
        end
        savefig(horzcat(file, 'Entire_Session_Plot_no_fit'))
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
    
    baseline_gcamp = {};
    baseline_iso = {};
    for e = 1:size(liquid_delivered_times, 2)
        cur_inds = baseline_dict{e};
        baseline_gcamp{e} = gcamp(cur_inds, 2);
        baseline_iso{e} = iso(cur_inds, 2);
    end
    %baseline_gcamp = gcamp(baseline_inds, 2);
    %baseline_iso = iso(baseline_inds, 2);
end
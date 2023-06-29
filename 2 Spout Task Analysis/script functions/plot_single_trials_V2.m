%% TODO: Make this plot both CS and US time

function [synced_trials] = plot_single_trials_V2(US_present, trials, alignment, ...
    event_label, event_times, laser_duration)
    
    
    if strcmp(alignment, 'events')
        file_name = event_label{2};
    else
        file_name = alignment;
    end
    if US_present
        dir_name = horzcat(file_name, ' trials');
    else
        if strcmp(alignment, 'left')
            correct_title = 'right';
        elseif strcmp(alignment, 'right')
            correct_title = 'left';
        else
            correct_title = horzcat('no ', file_name);
        end
        dir_name = horzcat(correct_title, ' trials');
    end
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    else
       % error('single trials directory already exists.')
        delete_dir = (input(horzcat(dir_name, ' directory already exists. Delete contents? y or n: '), 's'));
        if strcmp(delete_dir, 'y')
            rmdir(dir_name, 's')
            mkdir(dir_name);
        else
            return
        end
    end
    
    synced_trials = {};
   
    stim_happened = cell(size(event_times, 2), 1); %2 x 1
    stim_frame = []; %List of the events
    event_frame = cell(size(event_times, 2), 1); %2 x 1
    number_of_events = size(event_times, 2);
    event_times_clone = event_times;
    %if types_to_check == 4
    %    checks = trials_to_check;
    %else
    %    checks = 1:length(trials_to_check);
    %end
    checks = 1:size(trials, 2);
    if or(strcmp(alignment, 'left'), strcmp(alignment, 'right'))
        number_of_events = 1;
    end
    for event_num = 1:number_of_events    
        for p = checks
            %If there are no more events
%            if p > size(event_times_clone{event_num}, 2)
%                stim_happened{event_num, 1} = horzcat(stim_happened{event_num, 1},0);
%                event_frame{event_num} = horzcat(event_frame{event_num}, 0);
%                continue
%            end
            %If there are no more events
            
            if isempty(event_times{event_num})
                stim_happened{event_num, 1} = horzcat(stim_happened{event_num, 1},0);
                event_frame{event_num, 1} = horzcat(event_frame{event_num}, 0);
                continue
            end
            current_event_time = event_times{event_num}(1);  
            current_stim = 0;
            ind = 1;
            tri_num = p;
            %If this isn't the type of trial
            %if trial_types(p) ~= types_to_check
            %    continue
            %end
            g_trial = trials{p};
            times = g_trial(:,1);
            time_scale = g_trial(:,1);
            %If current event is before the start of the trial
            while and(current_event_time < time_scale(1), current_event_time > 0)
                event_times{event_num}(1) = [];
                current_event_time = event_times{event_num}(1);
                stim_happened{event_num, 1} = horzcat(stim_happened{event_num, 1},0);
                event_frame{event_num, 1} = horzcat(event_frame{event_num}, 0);
            end
            
            locate = time_scale > current_event_time;
            stim_ind = 1;
            if and(current_event_time > time_scale(1), current_event_time < time_scale(end))
                for loc = locate'
                    if loc == 0
                        stim_ind = stim_ind + 1;
                        continue
                    else
                        current_stim = 1;
                        stim_happened{event_num, 1} = horzcat(stim_happened{event_num, 1},1);
                        event_frame{event_num} = horzcat(event_frame{event_num}, time_scale(stim_ind) - time_scale(1));
                        event_times{event_num}(1) = [];
                        break
                    end
                end
            elseif current_event_time > time_scale(end)
            %Sometimes a trial is randomly ommitted? Event after trial.
            %Don't delete the trial
                stim_happened{event_num, 1} = horzcat(stim_happened{event_num, 1},0);
                event_frame{event_num} = horzcat(event_frame{event_num}, 0);
                disp('event is too late')
            end
        end
    end
    water_align = 0;
    tit = alignment;
    event_actual_times = event_frame{2};
    for p = checks
        %types_to_check = 0;
        %
        g_trial = trials{p};
        times = g_trial(:,1);
        time_scale = g_trial(:,1);
        figure
        hold on
        if water_align
            plot((g_trial(:,1)- times(1) - current_event_frame), g_trial(:,2), 'g'); %0 at time of stimulus
        else
            plot((g_trial(:,1)- times(1)), g_trial(:,2), 'g');
            
            
            if strcmp(tit, 'laser')
                yLimits = get(gca,'YLim');
                hold off
                close all
                figure
                hold on
                current_event_frame = event_actual_times(p);
                x1 = [current_event_frame, current_event_frame + laser_duration];
                y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
                Laser=area(x1, y1,'LineStyle','None');
                Laser(1).FaceColor = [0.976 0.894 0.690];
                x1 = [current_event_frame, current_event_frame + laser_duration];
                y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
                Laser=area(x1, y1,'LineStyle','None');
                Laser(1).FaceColor = [0.976 0.894 0.690];
                plot((g_trial(:,1)- times(1)), g_trial(:,2), 'g');
                set(gca, 'ylim',yLimits,'FontSize', 10);
                legend(Laser, {'Laser'})

            end
            
            
        end                
        yLimits = get(gca,'YLim');
        actual_event_num = number_of_events;
        for k = 1:number_of_events
            if size(stim_happened{k, 1}, 2) < p
                actual_event_num = actual_event_num - 1;
                continue
            end
            current_stim = stim_happened{k, 1}(p);
            current_event_frame = event_frame{k, 1}(p); %actually it;s time
            if ~strcmp(tit, 'laser')
                if and(current_stim ~= 0, ~water_align)
                    plot([current_event_frame, current_event_frame], yLimits, 'MarkerSize', 20);
                    ylim(yLimits);
                elseif and(water_align, current_stim)
                    plot([0, 0], yLimits, 'MarkerSize', 20);
                    ylim(yLimits);
                end
            end
        end
        %xlim([-10 10])
        label = horzcat(' ', event_label{1, 1:number_of_events});
        if or(strcmp(alignment, 'left'), strcmp(alignment, 'right'))
            label = '';
            tit = '';
        end
        title(horzcat('Trial # ', num2str(p), ', ', tit));
        %if types_to_check == 4
        legend(event_label{1, 1:(actual_event_num + 1)})
        %elseif or(types_to_check == 0, types_to_check == 1)
        %    legend(event_label{1, 1:(number_of_events + 1)})
        %elseif types_to_check == 2
         %   %legend(event_label{1, 1:(number_of_events)})
        %end
        
        current_dir = pwd;    
        test = horzcat(current_dir, horzcat('/', dir_name, '/Trial '), num2str(p), label, '.fig');
        savefig(horzcat(current_dir, horzcat('/', dir_name, '/Trial '), num2str(p), label, '.fig'));
        hold off
        ind = ind + 1;
        %close all
    end
    disp('halp')
    close all
end
function [control] = fit_control_channel(signal, control)
    
    if size(signal, 1) <= size(control, 1)
        control(1:((size(control, 1) - size(signal, 1))), :) = [];
    end   
    signal_clone = signal;
    green = signal(:,2);
    violet = control(:,2);
    close all
    
    if size(green, 1) > size(violet, 1)
        green = green(2:end);
        signal = signal(2:end, :);
    end 
    figure
    plot(signal(:,1), green,'g')
    hold on
    plot(control(:,1), violet,'m')
    
    title('signal and control plotted together with no fitting')
    legend('signal', 'control');
    savefig('Entire_Session_Plot_no_fit')
    hold off
    green_cor = green;
    violet_cor = violet;

    if 1
        

        sync = 1;
        if sync == 1
            %signal = green_cor - violet_cor;
            ts1 = timeseries(green_cor,signal(:,1));
            % Associate another frequency with amplitude
            ts2 = timeseries(violet_cor,control(:,1));
            % Synchronize the frequency
            [ts3, ts4] = synchronize(ts1,ts2,'Union');
            % Plot the difference
            %figure
            %hold on
            %plot(ts3.Time, ts3.Data-ts4.Data);
            %title("Post correction, with syncing");
            new_signal = horzcat(ts3.Time,ts3.Data); %-ts4.Data;
            %new_signal = horzcat(ts3.Time, smooth(signal, 3));
            %new_control = horzcat(ts3.Time, smooth(ts4.Data - 1, 3));
            %new_signal = horzcat(ts3.Time, signal);
            new_control = horzcat(ts4.Time, ts4.Data);
            %signal = new_signal;
            %control = new_control;

        else
            new_control = horzcat(control(:,1), violet_cor - 1);
            control = new_control;
            signal = movmean(green_cor, 2) - movmean(violet_cor, 2);
            signal(:,2) = signal;
            %signal(:,2) = movmean(signal,2);
        end
        figure
        plot(new_signal(:,1), new_signal(:,2),'g')
        hold on
        plot(new_control(:,1), new_control(:,2),'k')
        title('After syncing')
        legend('Signal','Control')
        savefig('post_sync_entire.fig')
        
        green_cor = new_signal(:,2) + 1;
        violet_cor = new_control(:,2) + 1;
        
        
        fun = @(x) sseval(x,green_cor,violet_cor);
        x0=1;     % starting guess; perhaps the rms ratio of the two would be better if far from 1
        bestk=fminsearch(fun,x0);
        violet_cor = violet_cor*bestk;
        green_cor = green_cor - 1;
        violet_cor = violet_cor - 1;
        
        
        %[coeff] = regress(green_cor,[ones(length(violet_cor),1) violet_cor]);
        %[coeff] = regress(green_cor, violet_cor);
        %violet_cor = [ones(length(violet_cor),1) violet_cor]*coeff;
        %violet_cor = violet_cor*coeff;
        %green_cor = green_cor*coeff;
        %violet_cor = violet_cor + 1;
        %green_cor = green_cor + 1;

        %Plot fitted graph
        figure
        plot(ts3.Time, green_cor,'g');
        hold on
        plot(ts4.Time, violet_cor,'m')
        %plot(water_delivered_times, zeros(length(water_delivered_times), 1), 'color', 'c', 'MarkerSize', 10, 'Marker', 'o')
        title('signal and control plotted together with fitting')
        legend('signal', 'control');
        savefig('Entire_Session_Plot_with_fit')
        hold off
    end
    
end
function sse=sseval(k,D1,D2)
    sse=sum((D1-k*D2).^2);
end

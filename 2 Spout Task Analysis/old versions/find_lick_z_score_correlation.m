function [] = find_lick_z_score_correlation()
%% Finds the lick rates per trial

    %% Set these variables
    close all;
    combine_spouts = 0; % IF YOU WANT TO COMBINE SPOUTS, SET THIS
    delete_outliers = 0; %Will still show you outliers but will not remove them
    
    %% DON'T TOUCH PAST HERE
    signal_folder = uigetdir2(pwd, 'select signal folder(s)');
    behavior_folder = dir('*Behavior Figures');
    
    %% Find z scores
    % Get a list of all files and folders in this folder.
    %signal_files = dir(horzcat(signal_folder.folder, '/', signal_folder.name));
    outlier_x = [];
    outlier_y = [];
    names = {};
    P_list = {};
    for r = 1:size(signal_folder, 2)
        signal_name = signal_folder{r};
        parens = strfind(signal_name, '/');
        signal_vol = signal_name(parens(end) + 1);
        if strcmp(signal_vol, 's')
            signal_vol = '';
        else
            signal_vol = horzcat(signal_vol, ' ');
        end
        signal_files = dir(horzcat(signal_folder{r}, '/'));
        behavior_files = dir(horzcat(behavior_folder.folder, '/', behavior_folder.name));

        % Get a logical vector that tells which is a directory.
        signal_dirFlags = [signal_files.isdir];
        behavior_dirFlags = [behavior_files.isdir];
        % Extract only those that are directories.
        signal_subFolders = signal_files(signal_dirFlags);
        behavior_subFolders = behavior_files(behavior_dirFlags);
        % Print folder names to command window.
        signal_z_scores = {};
        lick_data = {};
        lick_ind = 1;
        z_ind = 1;
        for k = 3 : length(signal_subFolders)
            if combine_spouts == 1
                if ~contains(signal_subFolders(k).name, 'both')
                    continue
                end
            else
                if contains(signal_subFolders(k).name, 'both')
                    continue
                end
            end
            cur_subfolder = horzcat(signal_subFolders(k).folder, '/', signal_subFolders(k).name, '/');
            z_score_file = horzcat(signal_subFolders(k).folder, '/', signal_subFolders(k).name, '/', 'artifacts_removed', '/', 'z_scores.mat');
            if exist(z_score_file)
                load(z_score_file);
            else
                error(horzcat('You need to run remove_high_artifact_trials_V3.m on: ', newline, signal_name));
            end
            signal_z_scores{z_ind} = z_peaks; 
            names{(r - 1) * (length(signal_subFolders) - 2) + z_ind} = signal_subFolders(k).name;
            z_ind = z_ind + 1;
            for p = 3:length(behavior_subFolders)
                if ~combine_spouts
                    name_to_check = horzcat(signal_vol, signal_subFolders(k).name);
                else
                    name_to_check = horzcat(signal_subFolders(k).name);
                end
                behavior_name = behavior_subFolders(p).name;
                if strcmp(name_to_check, behavior_name)
                    correct_b_file = behavior_subFolders(p).name;
                    behavior_loc = horzcat(behavior_subFolders(p).folder, '/', behavior_subFolders(p).name, '/');
                    the_b_file = dir(horzcat(behavior_loc, '*mat'));
                    load(horzcat(behavior_loc, the_b_file.name));
                    lick_data{lick_ind} = total_lps_within_window_per_trial;
                    lick_ind = lick_ind + 1;
                    break
                end
            end

            fprintf('Sub folder #%d = %s\n', k, signal_subFolders(k).name);
        end
        colors = {'g', 'b', 'r'};
        
        lick_ind = 1;
        for q = 3:length(signal_subFolders)
            if lick_ind > size(signal_z_scores, 2)
                break
            end
            if combine_spouts == 1
                if ~contains(signal_subFolders(q).name, 'both')
                    continue
                end
            else
                if contains(signal_subFolders(q).name, 'both')
                    continue
                end
            end
            
            t = signal_z_scores{lick_ind};
            if isempty(lick_data)
                error(horzcat('There is no lick data matching ', signal_subFolders(q).name,' ', signal_vol, ' folder'));
            end
            y = lick_data{lick_ind};
            
            plot(t, y, '.', 'color', colors{lick_ind + r - 1}, 'MarkerSize', 10)

            hold on
            TF_z_scores = isoutlier(t);
            TF_licks = isoutlier(y);
            outlier_check = TF_z_scores + TF_licks;
            to_remove = find(outlier_check == 1);

            out_zs = t(to_remove);
            out_licks = y(to_remove);
            if delete_outliers
                t(to_remove) = [];
                y(to_remove) = [];
            end
            outlier_x = horzcat(outlier_x, out_zs);
            outlier_y = horzcat(outlier_y, out_licks);

            [R,P] = corrcoef(horzcat(t',y'));
            P_list{(r - 1) * (length(signal_subFolders) - 2) + q - 2} = num2str(P(1,2));
            [re,m,b] = regression(t,y);
            x_a = min(t);
            x_b = max(t);
            y_a = m*(x_a) + b;
            y_b = m*(x_b) + b;
            plot([x_a x_b], [y_a y_b], 'color', colors{lick_ind + r - 1})
            lick_ind = lick_ind + 1;
        end
        
    end
    if delete_outliers
        plot(outlier_x, outlier_y,'.', 'color', 'black', 'MarkerSize', 20)
    end
    title('Pearson Correlation of lick rate vs peak z-score(dF/F)', 'FontSize', 15);
    xlabel('Z-score(dF/F) peak response');
    ylabel('Lick-rate in Hz')
    if size(names, 2) > 1
        if delete_outliers
            legend(horzcat(names{1}, ': P = ', P_list{1}), '', horzcat(names{2}, ': P = ', P_list{2}), '', 'Outliers')
        else
            legend(horzcat(names{1}, ': P = ', P_list{1}), '', horzcat(names{2}, ': P = ', P_list{2}))
        end
    else
        if delete_outliers
            legend(horzcat(names{1}, ': P = ', P_list{1}),'', 'Outliers')
        else
            legend(horzcat(names{1}, ': P = ', P_list{1}))
        end
    end
    final = gcf;
    if combine_spouts
        savefig(final, horzcat(signal_folder{1}, '/', 'combined_spouts_pearson.fig'));
    else
        savefig(final, horzcat(signal_folder{1}, '/', 'individual_pearson.fig'));
    end
    hold off
end
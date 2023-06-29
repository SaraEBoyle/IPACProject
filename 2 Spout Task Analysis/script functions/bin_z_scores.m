% Now bin the data you've collected.
function bin_z_scores(pre_time, post_time, bin_size)
file = pwd;
filelist = dir(fullfile(file, '*signal'));
control_filelist = dir(fullfile(file, '*control'));

control = 0;
for x = horzcat(filelist, control_filelist)
    signal_folder = horzcat(file, '/', x.name);
    leftfilelist = {};
    rightfilelist = {};
    
    files = dir(signal_folder);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    % Print folder names to command window.
    for k = 3 : length(subFolders)
      fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
      if contains(subFolders(k).name, 'left')
          leftfilelist{1, (size(leftfilelist, 2) + 1)} = subFolders(k).name;
      elseif contains(subFolders(k).name, 'right')
          rightfilelist{1, (size(rightfilelist, 2) + 1)} = subFolders(k).name;
          %rightfilelist = horzcat(rightfilelist, subFolders(k).name);
      end
    end
    
  
    % Load and bin left data
    
    for sub = leftfilelist
        sub = sub{1};
        if ~mod(control, 2) %when not control data
            load(horzcat(signal_folder, '/', sub, '/', 'artifacts_removed/z_scores.mat'))
        else %for control data
            new_signal_folder = strrep(signal_folder,'control', 'signal');
            load(horzcat(new_signal_folder, '/', sub, '/', 'artifacts_removed/z_scores.mat'))
            corrected_z_scores = control_corrected_z_scores;
        end
        
        % Now bin, making sure the stimulus bin is divided correctly
        pre_left_z_score_binned = [];
        post_left_z_score_binned = [];
        binned_x_axis = [];
        bin_origin = US_target_frame - 1;
        bin_origin_time = x_axis(bin_origin);
        post_bin_num = post_time/bin_size;
        pre_bin_num = pre_time/bin_size;
        prev_bin = 0;
        plot(x_axis, mean(corrected_z_scores, 1))
        hold on
        
        for p = 1:pre_bin_num
            single_bin = find(x_axis > (bin_origin_time - p*bin_size));
            updated_x = x_axis(single_bin);
            updated_z_scores = corrected_z_scores(:, single_bin);
            single_bin = find(updated_x < (bin_origin_time - prev_bin*bin_size));
            updated_x_axis = updated_x(single_bin);
            x_mean = mean(updated_x_axis);
            binned_x_axis = horzcat(x_mean, binned_x_axis);
            updated_z_scores = updated_z_scores(:, single_bin);
            z_mean = mean(mean(updated_z_scores));
            pre_left_z_score_binned = horzcat(mean(updated_z_scores, 2), pre_left_z_score_binned);
            plot(x_mean, z_mean, '*');
            prev_bin = p;
        end
        prev_bin = 0;
        for p = 1:post_bin_num
            single_bin = find(x_axis < (bin_origin_time + p*bin_size));
            updated_x = x_axis(single_bin);
            updated_z_scores = corrected_z_scores(:, single_bin);
            single_bin = find(updated_x > (bin_origin_time + prev_bin*bin_size));
            updated_x_axis = updated_x(single_bin);
            updated_z_scores = updated_z_scores(:, single_bin);
            z_mean = mean(mean(updated_z_scores));
            x_mean = mean(updated_x_axis);
            binned_x_axis = horzcat(binned_x_axis, x_mean);
            plot(x_mean, z_mean, '*');
            post_left_z_score_binned = horzcat(post_left_z_score_binned, mean(updated_z_scores, 2));
            prev_bin = p;
        end
        left_z_score_binned = horzcat(pre_left_z_score_binned, post_left_z_score_binned);
        %plot([-9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -.5, .5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], mean(left_z_score_binned, 1))
        savefig('left_binned_avg_data.fig');
        close all
        avg_left_z_score_binned = mean(left_z_score_binned, 1);
        if ~mod(control, 2) %when not control data
            save(horzcat(signal_folder, '/', sub, '/', 'artifacts_removed/binned_data.mat'), 'avg_left_z_score_binned', 'left_z_score_binned', 'binned_x_axis')
        else %for control data
            new_signal_folder = strrep(signal_folder,'control', 'signal');
            control_avg_left_z_score_binned = avg_left_z_score_binned;
            control_left_z_score_binned = left_z_score_binned;
            save(horzcat(new_signal_folder, '/', sub, '/', 'artifacts_removed/binned_control_data.mat'), 'control_avg_left_z_score_binned', 'control_left_z_score_binned', 'binned_x_axis')
        end
    end
    for sub = rightfilelist
        sub = sub{1};
        
        if ~mod(control, 2) %when not control data
            load(horzcat(signal_folder, '/', sub, '/', 'artifacts_removed/z_scores.mat'))
        else %for control data
            new_signal_folder = strrep(signal_folder,'control', 'signal');
            load(horzcat(new_signal_folder, '/', sub, '/', 'artifacts_removed/z_scores.mat'))
            corrected_z_scores = control_corrected_z_scores;
        end
        
        %load(horzcat(signal_folder, '/', sub, '/', 'artifacts_removed/z_scores.mat'))
        % Now bin, making sure the stimulus bin is divided correctly
        pre_right_z_score_binned = [];
        post_right_z_score_binned = [];
        binned_x_axis = [];
        bin_origin = US_target_frame - 1;
        bin_origin_time = x_axis(bin_origin);
        post_bin_num = post_time/bin_size;
        pre_bin_num = pre_time/bin_size;
        prev_bin = 0;
        plot(x_axis, mean(corrected_z_scores, 1))
        hold on
        
        for p = 1:pre_bin_num
            single_bin = find(x_axis > (bin_origin_time - p*bin_size));
            updated_x = x_axis(single_bin);
            updated_z_scores = corrected_z_scores(:, single_bin);
            single_bin = find(updated_x < (bin_origin_time - prev_bin*bin_size));
            updated_x_axis = updated_x(single_bin);
            x_mean = mean(updated_x_axis);
            binned_x_axis = horzcat(x_mean, binned_x_axis);
            updated_z_scores = updated_z_scores(:, single_bin);
            z_mean = mean(mean(updated_z_scores));
            pre_right_z_score_binned = horzcat(mean(updated_z_scores, 2), pre_right_z_score_binned);
            plot(x_mean, z_mean, '*');
            prev_bin = p;
        end
        prev_bin = 0;
        for p = 1:post_bin_num
            single_bin = find(x_axis < (bin_origin_time + p*bin_size));
            updated_x = x_axis(single_bin);
            updated_z_scores = corrected_z_scores(:, single_bin);
            single_bin = find(updated_x > (bin_origin_time + prev_bin*bin_size));
            updated_x_axis = updated_x(single_bin);
            updated_z_scores = updated_z_scores(:, single_bin);
            z_mean = mean(mean(updated_z_scores));
            x_mean = mean(updated_x_axis);
            binned_x_axis = horzcat(binned_x_axis, x_mean);
            plot(x_mean, z_mean, '*');
            post_right_z_score_binned = horzcat(post_right_z_score_binned, mean(updated_z_scores, 2));
            prev_bin = p;
        end
        right_z_score_binned = horzcat(pre_right_z_score_binned, post_right_z_score_binned);
        %plot([-9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -.5, .5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], mean(left_z_score_binned, 1))
        savefig('right_binned_avg_data.fig');
        close all
        avg_right_z_score_binned = mean(right_z_score_binned, 1);
        
        
        if ~mod(control, 2) %when not control data
            save(horzcat(signal_folder, '/', sub, '/', 'artifacts_removed/binned_data.mat'), 'avg_right_z_score_binned', 'right_z_score_binned', 'binned_x_axis')
        else %for control data
            new_signal_folder = strrep(signal_folder,'control', 'signal');
            control_avg_right_z_score_binned = avg_right_z_score_binned;
            control_right_z_score_binned = right_z_score_binned;
            save(horzcat(new_signal_folder, '/', sub, '/', 'artifacts_removed/binned_control_data.mat'), 'control_avg_right_z_score_binned', 'control_right_z_score_binned', 'binned_x_axis')
        end
        
        %save(horzcat(signal_folder, '/', sub, '/', 'artifacts_removed/binned_data.mat'), 'avg_right_z_score_binned', 'right_z_score_binned', 'binned_x_axis')
    end
    control = control + 1;
    %save('binned_data.mat', 'left_z_score_binned', 'avg_left_z_score_binned', 'right_z_score_binned', 'avg_right_z_score_binned', 'binned_x_axis')
end
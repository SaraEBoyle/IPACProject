
%% Code for analyzing photometry responses during the Two Spout Task. 
% Corrects photometry data for bleaching and integrates behavior recorded by bpod. 
% Segments data by type of liquid dispensed, aligning to the first recorded lick
% in each trial.
%
% Setup requirements:
%   1. bpod, with StateMachine-Bpod0_5 firmware or newer
%       - Ports 1 and 2 control liquid dispensing valves, 3 and 4 are
%       custom circuits measuring licks
%   2. Photometry system using Pointgrey CCD camera (We used Neurofiberphotometry)
%   3. Arduino Due used to read TTL signals from the bpod signifying the
%   start of new trials
%
% Inputs: 
%   1. Bpod .m file from delivering liquids using the "CSandUSforAle_WAIT.m/" protocol
%   2. 2 Bonsai .csv files from "2 fibers 1 video 2 spout delivery.bonsai":
%       1 will be a matrix containing signal and timing data, and the other
%       will contain analog input from a secondary arduino used to record
%       trial start signals from the bpod.
% Outputs:
%   1. Segments trial by trial responses to each type of liquid stimulus,
%   synced to the time the animals touch the spout. Outputs trial by trial
%   z scores, average z scores, and records lick rate data.
%% Set up your variables
to_analyze = {'left', 'right'}; %choose left spout, right spout, or both
right_spout = 'Right liquid';
left_spout = 'Left liquid';
left_lick_port = 'Port3In';
right_lick_port = 'Port4In';
pre_time = 10; % time before stimulus to plot
post_time = 20; % time after stimulus to plot
do_AUC = 1; % if 1, calculate area under the curve for stimulus response
trials_to_check = 1:50; %List of trials to analyze
fibers = [1]; %Choose which photometry fiber to look at (supports up to 2)
fiber_names = {'Left', 'Right'};
trigger1 = 1; %Change to 1 if your signal alternates through 2 channels
frame_rate = 20.0;
bin_size = .5; %seconds

%% Set these so the program can identify your input files
signal_key_word = 'SIG'; %Use a word that's in your signal files
analog_key_word = 'ANALOG'; %Use a word in your analog files
bpod_key_word = 'CSand'; %Use a word in the bpod files

%% choose what kind of bleaching correction you want
moving_avg = 0; %Does a moving average to subtract F
correct_for_bleaching = 1; %Subtracts biexponential curve to correct. Best for data that has large, long lasting fluctuations (recommended)

%% Set up variables for removing artifacts
left_window = 0; %time after stimulus to start screening for movement artifacts
right_window = 5; %time top stop screening

%% Run the main code, don't change anything below here
seconds_per_trial = pre_time + post_time;
seconds_to_extend = 0;
pre_ITI_length = 10;
leave_out_no_licks = 1; %Exclude trials where the mouse doesn't lick + aligne to first lick
lick_window = 30; %Length of time to search for first lick for alignment

% First analyze photometry data. This splits the data into signal and control
% channels. It will correct for bleaching, either with a rolling average or
% subtracting a fitted exponential. It then calculate dF/F and z-score
% using the first few seconds of each trial as baseline. This also aligns
% photometry data very closely to delivered stimuli and keeps track of
% alignment data, which is needed for further analysis steps
volume_to_check = 5; %1000 if no variation 
Photometry_Analysis_V11(to_analyze, right_spout, left_spout, left_lick_port, right_lick_port, pre_time, post_time, seconds_per_trial, seconds_to_extend, trials_to_check, fibers, fiber_names, volume_to_check, pre_ITI_length, trigger1, frame_rate, leave_out_no_licks, lick_window, signal_key_word, analog_key_word, bpod_key_word, moving_avg, correct_for_bleaching)


% Then find trials that have high amounts of activity in the control
% channel and remove them. Stores updated trial info in "artifacts_removed"
% folder. Also updates alignment data.
fps = frame_rate/2;
remove_high_artifact_trials_V5()
% This plots and reads lick data.

%% Variables for licking behavior 
left_window = 0;
right_window = 3;
graph_behavior_Ale_V5(right_spout, left_spout, left_lick_port, right_lick_port, left_window, right_window, pre_time, post_time, bpod_key_word);


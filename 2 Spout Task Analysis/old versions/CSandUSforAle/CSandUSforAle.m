function CSandUSforAle
%% Can give a CS followed by water or airpuff 
% 1: CS1 + Water, 2: CS2 plus air puff, 3: Both
% SETUP
% > Connect the left valve in the box to Bpod Port#1.
% > Connect the right valve in the box to Bpod Port#2.
% > Lickometer Left: Bpod Port#3
% > Lickometer Right: Bpod Port#4
% > Imaging: Wire 1 (output)
% > BNC3: laser
% > BNC4: camera trigger

clear
global BpodSystem
global S
global currentTrial

%% Define parameters
type = 3; %1 is left only, 2 is right only, 3 is both intermixed
nTrials =50;
secs_to_wait = 300; %Set this to wait for a lick to move on
LickPortLeft = 'Port3In';
LickPortRight = 'Port4In';
LeftValveState = 1; %That's port 1
RightValveState = 2; %Port 2
lick_threshold = 101; %number of trials with no licks before email alert
no_licks = 0; %Don't change
current_left_trial = 0;
current_right_trial = 0;
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S

if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    S.GUI.LeftAmount = 4; % ul of fluid for left spout
    S.GUI.RightAmount = 4; % ul of fluid for right spout
    S.GUI.VaryLeft = 0; %If 1, vary fluid from 2-10 uL, even numbers only
    S.GUI.VaryRight = 0; %If 1, vary fluid from 2-10 uL, even numbers only
    left_volumes = [4, 4]; %Set this to whichever volumes you want to randomly deliver
    right_volumes = [4, 4]; %ditto
    
    S.GUI.time_between_CS_US_left = 3; %Can't be 0, and this starts at onset of sound. 3 for type 3. 1 for cs and water pairing
    S.GUI.time_between_CS_US_right = 3; %Can't be 0, and this starts at onset of sound. 3 for type 3. 1 for cs and water pairing
    
    S.RndFlag = 1; % If 1, random interspersed trials. If 0, block 
    S.NonrandomMaxTrials = nTrials; %Total number of trials to run.
    if S.RndFlag == 0
        NumTrial = (S.NonrandomMaxTrials)/2;
        S.NonrandomSequence = [zeros(1,NumTrial),ones(1,NumTrial)]; % If nonrandom, sequence of the trials to run. 
    end
	S.GUI.AssistProb = 1; %Chance of reward on failure trials. 1 = always, 0 = never
    S.GUI.BaselineTime = 4; %Time befre beginning of trial until CS
    S.GUI.BaselineTimeLeft = 4; %Time befre beginning of trial until CS
    S.GUI.BaselineTimeRight = 4; %Time befre beginning of trial until CS
    S.GUI.ResponseTimeLeft = 3; % How much time to make a choice, %Doesn't matter because global timer 2 is time between cs onset and us
	S.GUI.ResponseTimeRight = 3;    
    S.GUI.PunishDelayMean = 0;
    S.PunishDelayMax = 0;
    S.PunishDelayMin = 0;
    S.GUI.RewardDelayMean = 0;
    S.RewardDelayMax = 0;
    S.RewardDelayMin = 0;
    S.CueDelay = 0; % the time from cue onset to response window
    S.ITI = 15; % vary ITI between 4 and 8, centered around 6 seconds
    S.ITI_min = 13; S.ITI_max =17;
    S.SoundDuration = 1;
end

% Initialize parameter GUI plugin
BpodParameterGUI('init', S);
TotalLeftDisplay('init');
TotalRightDisplay('init');
%% Define trials
% Randomise the training conditions. 
if type == 3 %If left and right spout 
    if  S.RndFlag
        RandNumber = randperm(nTrials);
        TrialTypes =  mod(RandNumber,2);
        if TrialTypes(1)==2
            RandNumber = randperm(nTrials);
            TrialTypes =  mod(RandNumber,2);
        end
        
        RandNumber = randperm(nTrials/length(left_volumes));
        order_left =  mod(RandNumber,2);
        if order_left(1)==2
            RandNumber = randperm(nTrials/length(left_volumes));
            order_left =  mod(RandNumber,2);
        end
        
        RandNumber = randperm(nTrials/length(right_volumes));
        order_right =  mod(RandNumber,2);
        if order_right(1)==2
            RandNumber = randperm(nTrials/length(right_volumes));
            order_right =  mod(RandNumber,2);
        end
        
    else    
        disp('Using nonrandom trial list.')
        TrialTypes = horzcat(zeros(1,NumTrial),ones(1, NumTrial));
    end
elseif type == 2 % if right only
        TrialTypes = ones(1, nTrials);
elseif type == 1 %if left only
    TrialTypes = zeros(1, nTrials);
elseif type == 4 %This is for just laser stim
    TrialTypes = ones(1, nTrials) .* 2;
end

R = repmat(S.GUI.PunishDelayMean,1,nTrials);
if S.PunishDelayMax>S.PunishDelayMin
    for k=1:nTrials
        candidate_delay = exprnd(S.GUI.PunishDelayMean);
        while candidate_delay>S.PunishDelayMax || candidate_delay<S.PunishDelayMin
            candidate_delay = exprnd(S.GUI.PunishDelayMean);
        end
        R(k) = candidate_delay;
    end
end
PunishDelay = R;

R = repmat(S.GUI.RewardDelayMean,1,nTrials);
if S.RewardDelayMax>S.RewardDelayMin
    for k=1:nTrials
        candidate_delay = exprnd(S.GUI.RewardDelayMean);
        while candidate_delay>S.RewardDelayMax || candidate_delay<S.RewardDelayMin
            candidate_delay = exprnd(S.GUI.RewardDelayMean);
        end
        R(k) = candidate_delay;
    end
end
RewardDelay = R;

R = repmat(S.ITI,1,nTrials);
for k=1:nTrials
    candidate_delay = exprnd(S.ITI);
    while candidate_delay>S.ITI_max || candidate_delay<S.ITI_min
        candidate_delay = exprnd(S.ITI);
    end
    R(k) = candidate_delay;
end
ITI = R;

%e = [];
%r = randi([1 9],1,1000);
%for q = r
%    if or(~mod(q, 3), q == 1)
%        e = horzcat(e, q);
%    end
%end

e_left = datasample(left_volumes,1000);
e_right = datasample(right_volumes, 1000);
BpodSystem.Data.TrialTypes = []; 
BpodSystem.Data.TrialRewarded = [];
BpodSystem.Data.RewardDelay = [];
BpodSystem.Data.PunishDelay = [];
BpodSystem.Data.ITI = [];

%% Initialize plots
BpodSystem.ProtocolFigures.SideOutcomePlotFig = figure('Position', [500 200 800 300],'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.SideOutcomePlot = axes('Position', [.075 .3 .89 .6]);
ThreeOutcomePlot(BpodSystem.GUIHandles.SideOutcomePlot,'init',TrialTypes, nTrials);

% Set soft code handler to trigger sounds
BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler_PlaySoundX';

SF = 192000; % Sound card sampling rate
SinWaveFreq1 = 4000;
sounddata1 = GenerateSineWave(SF, SinWaveFreq1, S.SoundDuration); % Sampling freq (hz), Sine frequency (hz), duration (s)
% sounddata2 = rand(1,S.CueDelay*SF); % white noise
% SinWaveFreq2 = 8000;
% sounddata2 = GenerateSineWave(SF, SinWaveFreq2, S.SoundDuration); % Sampling freq (hz), Sine frequency (hz), duration (s)
SinWaveFreq2 = 12000;
sounddata2 = GenerateSineWave(SF, SinWaveFreq2, S.SoundDuration); % Sampling freq (hz), Sine frequency (hz), duration (s)
% sounddata3 = (rand(1,SF*S.SoundDuration+1)*2) - 1;
% WidthOfFrequencies=1.5; NumberOfFrequencies=7; MeanSoundFreq4 = 6000; SoundRamping=0.2;
% sounddata4 = SoundGenerator(SF, MeanSoundFreq4, WidthOfFrequencies, NumberOfFrequencies, S.SoundDuration, SoundRamping);
sounddata5 = GenerateSineWave(SF, SinWaveFreq1, 0); % Sampling freq (hz), Sine frequency (hz), duration (s)
% Program sound server
PsychToolboxSoundServer('init')
PsychToolboxSoundServer('Load', 1, 0.03*sounddata1);
PsychToolboxSoundServer('Load', 2, 0.05*sounddata2);
PsychToolboxSoundServer('Load', 5, sounddata5);
% PsychToolboxSoundServer('Load', 3, 0.8*sounddata3);
% PsychToolboxSoundServer('Load', 3, sounddata3);
% PsychToolboxSoundServer('Load', 4, sounddata4);

%% Main trial loop

T = timer('period',1,'executionmode','fixedrate', 'TimerFcn', @scriptMonitor);
T.Period = 1;
T.ExecutionMode = 'fixedRate';
start(T);
autoTimerCleanup = onCleanup(@()delete(timerfind));

for currentTrial = 1:nTrials
    tic % Start timer
    S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
    L = GetValveTimes(S.GUI.LeftAmount, 1);
    R = GetValveTimes(S.GUI.RightAmount, 1);
    LeftValveTime = L; % Update reward amounts
    RightValveTime = R; % Update reward amounts
    USDelay = 0; % The delay caused by the unconditioned stimulus. 
    FSDelay = 0; % The difference between the two unconditioned stimuli. 
    RUSD = RewardDelay(currentTrial)+LeftValveTime; % Reward US delay
    PUSD = RightValveTime+PunishDelay(currentTrial); % Punish US delay
    
    switch TrialTypes(currentTrial) % Determine trial-specific state matrix fields
        case 0 % Left spout
            current_left_trial = current_left_trial + 1;
            time_between_CS_US = S.GUI.time_between_CS_US_left;
            soundID = 1;
            BaselineTime = S.GUI.BaselineTimeLeft;
            if S.GUI.VaryLeft == 1
                %S.GUI.LeftAmount = e_left(currentTrial); % ul
                volume_type = order_left(current_left_trial);
                if volume_type
                    S.GUI.LeftAmount = left_volumes(1); % ul
                else
                    S.GUI.LeftAmount = left_volumes(2); % ul
                end
            end
            ActionOutcome = 'LeftPre';
            ResponseTime = S.GUI.ResponseTimeLeft;
			ActionPort = LickPortLeft;
            S.CueDelay = 0;
            ResponseWindow = 'ResponseW';
			if (rand < S.GUI.AssistProb) %#ok<ALIGN> %Changed from >= Sara 1/29/19
				Tup_Action = 'DeliverLeft';
                licksOn = LickPortLeft;
            elseif 1
				Tup_Action = 'EndTime';
            end
          
		case 1 % right spout
            soundID = 2;
            current_right_trial = current_right_trial + 1;
            time_between_CS_US = S.GUI.time_between_CS_US_right;
            if S.GUI.VaryRight == 1
                volume_type = order_right(current_right_trial);
                %S.GUI.RightAmount = e_right(currentTrial); % ul
                if volume_type
                    S.GUI.RightAmount = right_volumes(1); % ul
                else
                    S.GUI.RightAmount = right_volumes(2); % ul
                end
            end
            BaselineTime = S.GUI.BaselineTimeRight;
            ActionOutcome = 'RightPre';
            ResponseTime = S.GUI.ResponseTimeRight;
			ActionPort = LickPortRight;
            S.CueDelay = 0;
            ResponseWindow = 'ResponseW';
			if rand < S.GUI.AssistProb %#ok<ALIGN> %Changed from >= Sara 1/29/19
				%Habituation phase, reward with p=AssistProb.
				Tup_Action = 'DeliverRight';
                licksOn = LickPortRight;
			else
				Tup_Action = 'EndTime';
            end
    end
    % Assemble state matrix
    sma = NewStateMatrix(); 
    sma = SetGlobalTimer(sma,  'TimerID', 1, 'Duration', 10, 'OnsetDelay', 0);%, 'Channel', 'Wire1');
    sma = SetGlobalTimer(sma, 'TimerID', 2, 'Duration', time_between_CS_US, 'OnsetDelay', 0);
% 'BNC1'); %set timer for CS and US interval
    sma = AddState(sma, 'Name', 'TrialStart', ...%1
        'Timer', BaselineTime,... %5 seconds
        'StateChangeConditions', {'Tup', 'StimulusDeliver'},...
        'OutputActions', {'GlobalTimerTrig', 1, 'Wire1', 1});

    sma = AddState(sma, 'Name', 'StimulusDeliver', ...%2
        'Timer', 0,... %1 second
        'StateChangeConditions', {'Tup', 'CueDelay'},...
        'OutputActions', {'GlobalTimerTrig', 2});%'SoftCode', soundID});
    sma = AddState(sma, 'Name', 'CueDelay', ...%3
        'Timer', 0,... %Time between sound onset and response window. 0
        'StateChangeConditions', {'Tup', 'ResponseW'},...
        'OutputActions', {});
    sma = AddState(sma, 'Name', 'ResponseW', ...%4
        'Timer', 3,... %This is 2 seconds
        'StateChangeConditions', {'GlobalTimer2_End', Tup_Action, 'Tup', Tup_Action},...
        'OutputActions', {});      
    sma = AddState(sma, 'Name', 'LeftPre', ...%6
        'Timer', ResponseTime,...
        'StateChangeConditions', {'GlobalTimer2_End', 'DeliverLeft', 'GlobalTimer1_End', 'ITI'},...
        'OutputActions', {});
    sma = AddState(sma, 'Name', 'RightPre', ...%8
        'Timer', ResponseTime,...
        'StateChangeConditions', {'GlobalTimer2_End', 'DeliverRight', 'GlobalTimer1_End', 'ITI'},...
        'OutputActions', {});
    sma = AddState(sma, 'Name', 'DeliverLeft', ...%13
        'Timer', LeftValveTime,...
        'StateChangeConditions', {'Tup', 'LickCheck'},...
        'OutputActions', {'ValveState', LeftValveState});
    sma = AddState(sma, 'Name', 'DeliverRight', ...%13
        'Timer', RightValveTime,...
        'StateChangeConditions', {'Tup', 'LickCheck'},...
        'OutputActions', {'ValveState', RightValveState});
    sma = AddState(sma, 'Name', 'Laser', ...%13
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'LickCheck'},...
        'OutputActions', {'BNC2', 1});
    sma = AddState(sma, 'Name', 'LickCheck', ...%15
        'Timer', secs_to_wait,...
        'StateChangeConditions', {'Tup', 'ITI', licksOn, 'EndTime'},...
        'OutputActions', {});
    sma = AddState(sma, 'Name', 'EndTime', ...%15
        'Timer', 10 - BaselineTime,...
        'StateChangeConditions', {'Tup', 'ITI', 'GlobalTimer1_End', 'ITI'},...
        'OutputActions', {});
    sma = AddState(sma, 'Name', 'ITI', ...%16
        'Timer', ITI(currentTrial),...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {});

    BpodSystem.ProtocolFigures.StartTime = 0;
	BpodSystem.ProtocolFigures.CurrentTrialStartTime = datenum(datetime)*86400-BpodSystem.ProtocolFigures.StartTime;
    SendStateMatrix(sma);
    RawEvents = RunStateMatrix;
	BpodSystem.Data.CurrentTrialLickTimes = [];
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned
		RE = RawEvents.Events;
		RES = RawEvents.States;
		REST = RawEvents.StateTimestamps(2:end);
		RET = RawEvents.EventTimestamps;
		%Find if any timestamps contain a reward deliver (event 13) and
		%record that as the current trial reward delivery timestamp.
		RESR = (RES==13);
		BpodSystem.Data.CurrentTrialRewardTime = round(BpodSystem.ProtocolFigures.CurrentTrialStartTime+sum(RESR.*REST));
		%Filter the timestamps for licks (In3Up, 53/30 depending on computer) 
		%and turn them into a list.
		REE = (RE==32);
		BpodSystem.Data.CurrentTrialLickTimes = round(nonzeros(REE.*RET)+BpodSystem.ProtocolFigures.CurrentTrialStartTime);
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        BpodSystem.Data.TrialTypes(currentTrial) = TrialTypes(currentTrial); % Adds the trial type of the current trial to data
        BpodSystem.Data.RewardDelay(currentTrial) = RewardDelay(currentTrial);
        BpodSystem.Data.PunishDelay(currentTrial) = PunishDelay(currentTrial);
        BpodSystem.Data.ITI(currentTrial) = ITI(currentTrial);
        if or(isfield(BpodSystem.Data.RawEvents.Trial{1,currentTrial}.Events, 'Port3In'), isfield(BpodSystem.Data.RawEvents.Trial{1,currentTrial}.Events, 'Port4In'))
            no_licks = 0;
        else
            no_licks = no_licks + 1;
        end
        
        if no_licks == lick_threshold;
            setpref('Internet','SMTP_Server','smtp.gmail.com');
            setpref('Internet','E_mail','beehaviorlab@gmail.com');
            setpref('Internet','SMTP_Username','beehaviorlab');
            setpref('Internet','SMTP_Password','beehaviorlab85');
            props = java.lang.System.getProperties;
            props.setProperty('mail.smtp.auth','true');
            props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
            props.setProperty('mail.smtp.socketFactory.port','465');
            sendmail('beehaviorlab@gmail.com','Mouse is not licking') ;
        end
        
        %% Calculate Outcome
        if type == 1
            if isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.LeftPre(1))
                BpodSystem.Data.Outcomes(currentTrial) = 3; %Receives left without licking to tone
            elseif ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverLeft(1))
                BpodSystem.Data.Outcomes(currentTrial) = 1; %Receives left
            elseif TrialTypes(currentTrial)==1
                BpodSystem.Data.Outcomes(currentTrial) = -1; %Did not get left
            else
                BpodSystem.Data.Outcomes(currentTrial) = 2; %Ignored no-go
            end
        elseif type == 2
            if ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverRight(1))
                BpodSystem.Data.Outcomes(currentTrial) = 0; %Receives right
            elseif TrialTypes(currentTrial)==1
                BpodSystem.Data.Outcomes(currentTrial) = -1; %doesn't get right
            else
                BpodSystem.Data.Outcomes(currentTrial) = 2; %Ignored reward CS
            end
        elseif type == 3
            if TrialTypes(currentTrial) == 0
                if and(isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.LeftPre(1)), ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverLeft(1)))
                    BpodSystem.Data.Outcomes(currentTrial) = 1; %Receives left reward
                elseif ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.LeftPre(1))
                    BpodSystem.Data.Outcomes(currentTrial) = 1; %Receives left reward
                else
                    BpodSystem.Data.Outcomes(currentTrial) = 2; %Ignored no-go
                end
            elseif TrialTypes(currentTrial) == 1
                 if ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverRight(1))
                    BpodSystem.Data.Outcomes(currentTrial) = 0; %Receives right reward
                 else
                    BpodSystem.Data.Outcomes(currentTrial) = -1; %Ignored go
                 end
            end
        elseif type == 4
          
            BpodSystem.Data.Outcomes(currentTrial) = 1; %Ignored no-go
           
        end
        
        UpdateTotalLeftDisplay(S.GUI.LeftAmount, currentTrial);
        UpdateTotalRightDisplay(S.GUI.RightAmount, currentTrial);
        UpdateThreeOutcomePlot(TrialTypes, BpodSystem.Data, type); %Sara
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
    if BpodSystem.Status.BeingUsed == 0
        delete(T); %added by Sara 1/29/19
        close all;
        return
    end

end

%{
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','SaraElizaBoyle@gmail.com');
setpref('Internet','SMTP_Username','SaraElizaBoyle');
setpref('Internet','SMTP_Password','IDoNotCare4');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('SaraElizaBoyle@gmail.com','Training is done!') ;
%}
graph_behavior_from_SessionData(BpodSystem.Path.CurrentDataFile)
%

function UpdateThreeOutcomePlot(~, Data, type)
	Outcomes = zeros(1,Data.nTrials);
  	for x = 1:Data.nTrials    
        if type == 1
            if and(isnan(BpodSystem.Data.RawEvents.Trial{x}.States.RewardPre(1)), ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.DeliverReward(1)))
                if ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.LickSuccess(1))
                    Outcomes(x) = 1; %Receives reward
                else
                    Outcomes(x) = 3; %Receives reward w assist
                end
                
                
            elseif ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.DeliverReward(1))
                Outcomes(x) = 1; %Receives reward
            elseif TrialTypes(x)==1
                Outcomes(x) = -1; %Ignored go
            else
                Outcomes(x) = 2; %Ignored no-go
            end
        elseif type == 2
            if ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.DeliverPunishment(1))
                Outcomes(x) = 0; %Receives punishement
            elseif TrialTypes(x)==1
                Outcomes(x) = -1; %avoided
            else
                Outcomes(x) = 2; %Ignored rewrd CS
            end
        elseif type == 3
            if TrialTypes(x) == 0
                if and(isnan(BpodSystem.Data.RawEvents.Trial{x}.States.LeftPre(1)), ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.DeliverLeft(1)))
                    Outcomes(x) = 3; %Receives reward w assist
                elseif ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.LeftPre(1))
                    Outcomes(x) = 1; %Receives reward
                else
                    Outcomes(x) = 2; %Ignored no-go
                end
            elseif TrialTypes(x) == 1
                 if ~isnan(BpodSystem.Data.RawEvents.Trial{x}.States.DeliverRight(1))
                        Outcomes(x) = 0; %Receives right liquid
                 else
                    Outcomes(x) = -1; %Ignored go
                 end
            end
        elseif type == 4
            Outcomes(x) = 1; %Ignored no-go
        end
    end
    ThreeOutcomePlot(BpodSystem.GUIHandles.SideOutcomePlot,'update',Data.nTrials+1,TrialTypes,Outcomes);
end

end

function UpdateTotalLeftDisplay(RewardAmount, currentTrial)
% If rewarded based on the state data, update the TotalRewardDisplay
global BpodSystem
    HasReward = 0;
    if (~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverLeft(1)))
        HasReward = 1;
    end
    TotalLeftDisplay('add', RewardAmount*HasReward);
end

function UpdateTotalRightDisplay(RewardAmount, currentTrial)
% If rewarded based on the state data, update the TotalRewardDisplay
global BpodSystem
    HasReward = 0;
    if (~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverRight(1)))
        HasReward = 1;
    end
    TotalRightDisplay('add', RewardAmount*HasReward);
end

function scriptMonitor(~,~)
	%Runs every second to monitor background events.
	global BpodSystem
    global S
    TotalRewardDisplay('update');
    clear BpodSystem;
end





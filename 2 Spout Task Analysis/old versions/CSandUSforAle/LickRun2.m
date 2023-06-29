function LickRun2

% SETUP
% > Connect the water valve in the box to Bpod Port#1.
% > Connect the air valve in the box to Bpod Port#2.
% > Lick: Bpod Port#3
% > Running: Wire 1

global BpodSystem

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    S.GUI.TrainingLevel = 3; % Configurable training level: 'Habituation','Active Avoidance','Imaging Active Avoidance'
    S.GUI.RewardAmount = 10; % ul
    S.GUI.PunishAmount = 0.2; % s (air puff)
    S.RndFlag = 1; % Deploy the random trials or repeat a list of trials. 
    S.NonrandomMaxTrials = 500; %Total number of trials to run. This should not exceed the limit of the grapher (BpodSystem.PluginObjects.maxTrials)
    S.NonrandomSequence = [ones(1,200),zeros(1,100)]; % If nonrandom, sequence of the trials to run. 
    S.ReTeaching = 0;
	S.GUI.AssistProb = 0.5;
    S.GUI.TrialGoProb = 0.5;
    S.GUI.PreGoTrialNum = 0;
    S.GUI.BaselineTime = 10;
    S.arduinoPort = 'COM9'; % Port (USB) of the wheel run sensor Arduino
	S.wheelCirc = 46; % Circumference of the wheel (cm)
	S.wheelTicks = 1024; % Number of ticks recorded by the wheel per rotation.
	S.wheelMult = S.wheelCirc/S.wheelTicks; % Distance per tick (cm)
    %Attempt to connect to the Arduino run detector. 
	BpodSystem.PluginObjects.ardSerial = serial(S.arduinoPort);
    BPAS = BpodSystem.PluginObjects.ardSerial;
	BPAS.InputBufferSize = 16000; % read up to 2kb
	set(BPAS,'DataBits',8);
	set(BPAS,'StopBits',1);
	set(BPAS,'BaudRate',115200);
	set(BPAS,'Parity','none');
	fopen(BPAS);
	autoCleanup = onCleanup(@()fclose(BPAS));
    S.GUI.ResponseTimeGo = 2; % How long until the mouse must make a choice, or forefeit the trial
	S.GUI.ResponseTimeNoGo = 2; % How long until the mouse must make a choice, or forefeit the trial
	S.GUI.ResponseTimeNoGoDelay = 0; %Wait before punishing the mouse.
    
    S.GUIMeta.TrainingLevel.Style = 'popupmenu'; % the GUIMeta field is used by the ParameterGUI plugin to customize UI objects.
    S.GUIMeta.TrainingLevel.String = {'Habituation','Active Avoidance','Imaging Active Avoidance'};
        
    S.GUI.PunishDelayMean = 0.1;
    S.PunishDelayMax = 0.1;
    S.PunishDelayMin = 0.1;
    S.GUI.RewardDelayMean = 0.1;
    S.RewardDelayMax = 0.1;
    S.RewardDelayMin = 0.1;

    S.CueDelay = 1.0; % the time from cue to response
    S.ITI = 20;
    S.ITI_min=16; S.ITI_max=24;
    S.SoundDuration = 1.0;
    %Modify the response timings with relation to the trial type. 
    if S.GUI.TrainingLevel==2
       %Non-imaging case. Reduce wait timings. 
       S.ITI = 6;
       S.ITI_min=5; S.ITI_max=7;
       S.GUI.BaselineTime = 0.1;
    else
       S.GUI.BaselineTime = 10;
    end
    
end

LickPort = 'Port3In';
RunPort = 'Wire1High';
RewardValveState = 1;
PunishValveState = 2;

% Initialize parameter GUI plugin
BpodParameterGUI('init', S);
TotalRewardDisplay('init');

%% Define trials
if (~S.RndFlag && BpodSystem.PluginObjects.maxTrials < S.NonrandomMaxTrials)
    MaxTrials = BpodSystem.PluginObjects.maxTrials;
else
    MaxTrials = S.NonrandomMaxTrials;
end
if S.GUI.TrainingLevel==1
    TrialTypes = ones(1,MaxTrials);
elseif S.GUI.TrainingLevel>=2
    % Randomise the training conditions. 
    TrialTypes = ones(1,MaxTrials);
    for ii=(S.GUI.PreGoTrialNum+1):MaxTrials
        if rand<S.GUI.TrialGoProb
            TrialTypes(ii) = 1;
        else
            TrialTypes(ii) = 0;
        end
    end
    
    if ~S.RndFlag
        disp('Using nonrandom trial list.')
%         seqType = [1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0];
        seqType = S.NonrandomSequence; 
        lenType = length(seqType);
        for ii=1:MaxTrials
            seqPos = mod(ii,lenType);
            if seqPos == 0
               seqPos = lenType;
            end
            TrialTypes(ii) = seqType(seqPos);
        end
    end
    
end

R = repmat(S.GUI.PunishDelayMean,1,MaxTrials);
if S.PunishDelayMax>S.PunishDelayMin
    for k=1:MaxTrials
        candidate_delay = exprnd(S.GUI.PunishDelayMean);
        while candidate_delay>S.PunishDelayMax || candidate_delay<S.PunishDelayMin
            candidate_delay = exprnd(S.GUI.PunishDelayMean);
        end
        R(k) = candidate_delay;
    end
end
PunishDelay = R;

R = repmat(S.GUI.RewardDelayMean,1,MaxTrials);
if S.RewardDelayMax>S.RewardDelayMin
    for k=1:MaxTrials
        candidate_delay = exprnd(S.GUI.RewardDelayMean);
        while candidate_delay>S.RewardDelayMax || candidate_delay<S.RewardDelayMin
            candidate_delay = exprnd(S.GUI.RewardDelayMean);
        end
        R(k) = candidate_delay;
    end
end
RewardDelay = R;

R = repmat(S.ITI,1,MaxTrials);
for k=1:MaxTrials
    candidate_delay = exprnd(S.ITI);
    while candidate_delay>S.ITI_max || candidate_delay<S.ITI_min
        candidate_delay = exprnd(S.ITI);
    end
    R(k) = candidate_delay;
end
ITI = R;

BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.
BpodSystem.Data.TrialRewarded = []; % The trial type of each trial completed will be added here.
BpodSystem.Data.RewardDelay = [];
BpodSystem.Data.PunishDelay = [];
BpodSystem.Data.ITI = [];
BpodSystem.Data.RunEvents = [];
BpodSystem.Data.RunSpeed = {};
BpodSystem.Data.CurrentTrialRunSpeed = zeros(10000,2);
BpodSystem.Data.CurrentTrialRunEvents = 0;
%% Initialize plots
BpodSystem.ProtocolFigures.SideOutcomePlotFig = figure('Position', [500 200 800 300],'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.SideOutcomePlot = axes('Position', [.075 .3 .89 .6]);
% BpodSystem.ProtocolFigures.LickPlotFig = figure('Position', [600 200 600 200],'name','Licking','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
GoNoGoOutcomePlot(BpodSystem.GUIHandles.SideOutcomePlot,'init',TrialTypes);

% Set soft code handler to trigger sounds
BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler_PlaySoundX';

SF = 192000; % Sound card sampling rate
SinWaveFreq1 = 3000;
sounddata1 = GenerateSineWave(SF, SinWaveFreq1, S.SoundDuration); % Sampling freq (hz), Sine frequency (hz), duration (s)
SinWaveFreq2 = 10000;
sounddata2 = GenerateSineWave(SF, SinWaveFreq2, S.SoundDuration); % Sampling freq (hz), Sine frequency (hz), duration (s)
% sounddata3 = (rand(1,SF*S.SoundDuration+1)*2) - 1;
% WidthOfFrequencies=1.5; NumberOfFrequencies=7; MeanSoundFreq4 = 6000; SoundRamping=0.2;
% sounddata4 = SoundGenerator(SF, MeanSoundFreq4, WidthOfFrequencies, NumberOfFrequencies, S.SoundDuration, SoundRamping);

% Program sound server
PsychToolboxSoundServer('init')
PsychToolboxSoundServer('Load', 1, 0.25*sounddata1);
PsychToolboxSoundServer('Load', 2, 1.2*sounddata2);
% PsychToolboxSoundServer('Load', 3, sounddata3);
% PsychToolboxSoundServer('Load', 4, sounddata4);

%% Main trial loop
T = timer('period',1,'executionmode','fixedrate','taskstoexecute',9999999);
T.TimerFcn = @(~,~)scriptMonitor();
T.Period = 1;
T.ExecutionMode = 'fixedRate';
start(T);
autoTimerCleanup = onCleanup(@()delete(timerfind));
for currentTrial = 1:MaxTrials
    S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
    R = GetValveTimes(S.GUI.RewardAmount, [1]);
    RewardValveTime = R; % Update reward amounts
    USDelay = 0; % The delay caused by the unconditioned stimulus. 
    FSDelay = 0; % The difference between the two unconditioned stimuli. 
    RUSD = RewardDelay(currentTrial)+RewardValveTime; % Reward US delay
    PUSD = S.GUI.PunishAmount+PunishDelay(currentTrial); % Punish US delay
    
    switch TrialTypes(currentTrial) % Determine trial-specific state matrix fields
        case 1 % lick trial
            soundID = 1;
            ActionOutcome = 'RewardPre';
            ResponseTime = S.GUI.ResponseTimeGo;
			ActionPort = LickPort;
			if (S.GUI.TrainingLevel >= 2 && rand<S.GUI.AssistProb)
				%Habituation phase, always provide reward with p=AssistProb.
				Tup_Action = 'AssistRewardPre';
			else
				Tup_Action = 'USSync';
            end
            USDelay = RUSD;
            if RUSD < PUSD
                FSDelay = PUSD-RUSD;
            end
            
		case 0 % run to avoid air puff trial
            soundID = 2;
			ActionPort = RunPort;
			ActionOutcome = 'DelayPre';
			ResponseTime = S.GUI.ResponseTimeNoGo;
			Tup_Action = 'Punishment';
            USDelay = PUSD;
            if RUSD > PUSD
                FSDelay = RUSD-PUSD;
            end
    end
    
    sma = NewStateMatrix(); % Assemble state matrix
    
    if S.GUI.TrainingLevel==1 % habituation
        
        sma = AddState(sma, 'Name', 'TrialStart', ...%1
            'Timer', 0.1,... % time before trial start
            'StateChangeConditions', {'Tup', 'ResponseW'},...
            'OutputActions', {});
        sma = AddState(sma, 'Name', 'ResponseW', ...%2
            'Timer', 15,... % reponse time window
            'StateChangeConditions', {ActionPort, 'Reward'},...
            'OutputActions', {});
        sma = AddState(sma, 'Name', 'AssistReward', ...%3
            'Timer', RewardDelay(currentTrial),...
            'StateChangeConditions', {'Tup', 'DeliverReward'},...
            'OutputActions', {});
        sma = AddState(sma, 'Name', 'Reward', ...%4
            'Timer', 0,... % reward delay
            'StateChangeConditions', {'Tup', 'DeliverReward'},...
            'OutputActions', {});
        sma = AddState(sma, 'Name', 'DeliverReward', ...%5
            'Timer', RewardValveTime,... % reward amount
            'StateChangeConditions', {'Tup', 'ITI'},...
            'OutputActions', {'ValveState', RewardValveState}); % 'SoftCode', soundID
        sma = AddState(sma, 'Name', 'ITI', ...%6
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {});
        
    else % unlearning and active avoidance
        % BNC1: imaging trigger
        % BNC2: ChR2 stimulation trigger
        timerTimeout = S.GUI.BaselineTime+ResponseTime+S.CueDelay;
        sma = SetGlobalTimer(sma, 1,timerTimeout);
        %Global timer for synchronising events. Timer1: Sound, Timer2: Actions. 
        sma = AddState(sma, 'Name', 'TrialStart', ...%1
            'Timer', S.GUI.BaselineTime,...
            'StateChangeConditions', {'Tup', 'StimulusDeliver'},...
            'OutputActions', {'BNC1',1,'GlobalTimerTrig', 1});
        sma = AddState(sma, 'Name', 'StimulusDeliver', ...%2
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'CueDelay'},...
            'OutputActions', {'BNC1',1,'SoftCode', soundID});
        sma = AddState(sma, 'Name', 'CueDelay', ...%3
            'Timer', S.CueDelay,...
            'StateChangeConditions', {'Tup', 'ResponseW'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'ResponseW', ...%4
            'Timer', ResponseTime,...
            'StateChangeConditions', {ActionPort, ActionOutcome, 'Tup', Tup_Action},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'AssistRewardPre', ...%5
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'AssistReward'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'RewardPre', ...%6
            'Timer', ResponseTime,...
            'StateChangeConditions', {'Tup', 'Reward','GlobalTimer1_End', 'Reward'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'PunishmentPre', ...%7
            'Timer', ResponseTime,...
            'StateChangeConditions', {'Tup', 'Punishment','GlobalTimer1_End', 'Punishment'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'DelayPre', ...%8
            'Timer', ResponseTime,...
            'StateChangeConditions', {'Tup', 'USSync','GlobalTimer1_End', 'USSync'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'USSync', ...%9
            'Timer', USDelay,...
            'StateChangeConditions', {'Tup', 'FinalSync'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'AssistReward', ...%9
            'Timer', RewardDelay(currentTrial),...
            'StateChangeConditions', {'Tup', 'DeliverReward'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'Reward', ...%10
            'Timer', RewardDelay(currentTrial),...
            'StateChangeConditions', {'Tup', 'DeliverReward'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'Punishment', ...%11
            'Timer', PunishDelay(currentTrial),...
            'StateChangeConditions', {'Tup', 'DeliverPunishment'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'DeliverPunishment', ...%12
            'Timer', S.GUI.PunishAmount,...
            'StateChangeConditions', {'Tup', 'FinalSync'},...
            'OutputActions', {'BNC1',1,'ValveState', PunishValveState});
        sma = AddState(sma, 'Name', 'DeliverReward', ...%13
            'Timer', RewardValveTime,...
            'StateChangeConditions', {'Tup', 'FinalSync'},...
            'OutputActions', {'BNC1',1,'ValveState', RewardValveState});
        sma = AddState(sma, 'Name', 'FinalSync', ...%14
            'Timer', FSDelay,...
            'StateChangeConditions', {'Tup', 'EndTime'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'EndTime', ...%15
            'Timer', S.GUI.BaselineTime,...
            'StateChangeConditions', {'Tup', 'ITI'},...
            'OutputActions', {'BNC1',1});
        sma = AddState(sma, 'Name', 'ITI', ...%16
            'Timer', ITI(currentTrial),...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {'BNC1',0});
   	end
    
	BpodSystem.ProtocolFigures.CurrentTrialStartTime = datenum(datetime)*86400-BpodSystem.ProtocolFigures.StartTime;
    BpodSystem.Data.CurrentTrialRunSpeed = zeros(10000,2);
    BpodSystem.Data.CurrentTrialRunEvents = 0;
    SendStateMatrix(sma);
    RawEvents = RunStateMatrix;
	BpodSystem.Data.CurrentTrialLickTimes = [];
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned
		%disp(RawEvents);
		RE = RawEvents.Events;
        %disp(RE);
		RES = RawEvents.States;
        %disp(RES);
		REST = RawEvents.StateTimestamps(2:end);
		%disp(REST);
		RET = RawEvents.EventTimestamps;
        %disp(RET);
		%Find if any timestamps contain a reward deliver (event 13) and
		%record that as the current trial reward delivery timestamp.
		RESR = (RES==13);
		%disp(RESR);
		BpodSystem.Data.CurrentTrialRewardTime = round(BpodSystem.ProtocolFigures.CurrentTrialStartTime+sum(RESR.*REST));
		%Filter the timestamps for licks (In3Up, 53/30 depending on computer) 
		%and turn them into a list.
		REE = (RE==32);
		BpodSystem.Data.CurrentTrialLickTimes = round(nonzeros(REE.*RET)+BpodSystem.ProtocolFigures.CurrentTrialStartTime);
        %disp(BpodSystem.Data.CurrentTrialLickTimes)
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        BpodSystem.Data.TrialTypes(currentTrial) = TrialTypes(currentTrial); % Adds the trial type of the current trial to data
        BpodSystem.Data.RewardDelay(currentTrial) = RewardDelay(currentTrial);
        BpodSystem.Data.PunishDelay(currentTrial) = PunishDelay(currentTrial);
        BpodSystem.Data.ITI(currentTrial) = ITI(currentTrial);
        %Outcome
        if ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.AssistReward(1))
			BpodSystem.Data.Outcomes(currentTrial) = 3; %Receives reward
		elseif ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.Reward(1))
			BpodSystem.Data.Outcomes(currentTrial) = 1; %Receives reward
        elseif ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.Punishment(1))
            BpodSystem.Data.Outcomes(currentTrial) = 0; %Receives punishement
        elseif TrialTypes(currentTrial)==1
            BpodSystem.Data.Outcomes(currentTrial) = -1; %Ignored go
        else
            BpodSystem.Data.Outcomes(currentTrial) = 2; %Ignored no-go
        end
        
        if S.ReTeaching==1 && S.GUI.TrainingLevel<=2 % full task or running task
            if BpodSystem.Data.Outcomes(currentTrial) == 0
                TrialTypes(currentTrial+1)=0;
            elseif BpodSystem.Data.Outcomes(currentTrial) == -1
                TrialTypes(currentTrial+1)=1;
            end
        end
        
        %Save the run timings and speeds
        BpodSystem.Data.RunEvents(currentTrial) = BpodSystem.Data.CurrentTrialRunEvents;
        BpodSystem.Data.RunSpeed{currentTrial} = BpodSystem.Data.CurrentTrialRunSpeed;
        
        UpdateTotalRewardDisplay(S.GUI.RewardAmount, currentTrial);
        UpdateGoNoGoOutcomePlot(TrialTypes, BpodSystem.Data);
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
    if BpodSystem.Status.BeingUsed == 0
        return
    end
end

function UpdateGoNoGoOutcomePlot(~, Data)
	Outcomes = zeros(1,Data.nTrials);
	for x = 1:Data.nTrials
		if ~isnan(Data.RawEvents.Trial{x}.States.AssistReward(1))
			Outcomes(x) = 3;
		elseif ~isnan(Data.RawEvents.Trial{x}.States.Reward(1))
			Outcomes(x) = 1;
		elseif ~isnan(Data.RawEvents.Trial{x}.States.Punishment(1))
			Outcomes(x) = 0;
		elseif BpodSystem.Data.TrialTypes(x)==1
			Outcomes(x) = -1;
		else
			Outcomes(x) = 2;
		end
    end
    GoNoGoOutcomePlot(BpodSystem.GUIHandles.SideOutcomePlot,'update',Data.nTrials+1,TrialTypes,Outcomes);
end

end

function UpdateTotalRewardDisplay(RewardAmount, currentTrial)
% If rewarded based on the state data, update the TotalRewardDisplay
global BpodSystem
    HasReward = 0;
    if (~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.DeliverReward(1)))
        HasReward = 1;
    end
    TotalRewardDisplay('add', RewardAmount*HasReward);
end
	
function scriptMonitor()
	%Runs every second to monitor background events. 
	global BpodSystem
    S = BpodSystem.ProtocolSettings;
	%Reads the data from the Arduino serial port. 
    BPAS = BpodSystem.PluginObjects.ardSerial;
	dataSize = get(BPAS,'BytesAvailable');
	if (dataSize>0)
		inASCII = fread(BPAS,dataSize);
		inData = transpose(char(inASCII));
		inData = strsplit(inData);
		inData = inData(1:end-1); % Remove the last newline from the array
		for k=1:length(inData) % Convert the array to run speed / run distance and write it to the run speed array. 
			BpodSystem.Data.CurrentTrialRunSpeed(k+BpodSystem.Data.CurrentTrialRunEvents,1:2) = str2double(strsplit(inData{k},'A'));
        end
        BpodSystem.Data.CurrentTrialRunEvents = BpodSystem.Data.CurrentTrialRunEvents + length(inData);
		assignin('base', 'inString', BpodSystem.Data.CurrentTrialRunEvents);
	end
	TotalRewardDisplay('update');
end

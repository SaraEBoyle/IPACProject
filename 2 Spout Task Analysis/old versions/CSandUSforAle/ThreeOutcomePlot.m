%% function OutcomePlot(AxesHandle,TrialTypeSides, OutcomeRecord, CurrentTrial)

function ThreeOutcomePlot(AxesHandle, Action, varargin)
%% 
% Plug in to Plot Go/No-Go side and trial outcome.
% For non-sided trial types, use the TrialTypeOutcomePlot plugin.
% AxesHandle = handle of axes to plot on
% Action = specific action for plot, "init" - initialize OR "update" -  update plot

%Example usage:
% SideOutcomePlot(AxesHandle,'init',TrialTypeSides)
% SideOutcomePlot(AxesHandle,'init',TrialTypeSides,'ntrials',90)
% SideOutcomePlot(AxesHandle,'update',CurrentTrial,TrialTypeSides,OutcomeRecord)

% varargins:
% TrialTypeSides: Vector of 0's (right) or 1's (left) to indicate reward side (0,1), or 'None' to plot trial types individually
% OutcomeRecord:  Vector of trial outcomes
%                 Simplest case: 
%                               1: correct trial (green)
%                               0: incorrect trial (red)
%                 Advanced case: 
%                               NaN: future trial (blue)
%                                -1: withdrawal (red circle)
%                                 0: incorrect choice (red dot)
%                                 1: correct choice (green dot)
%                                 2: did not choose (green circle)
% OutcomeRecord can also be empty
% Current trial: the current trial number

% Adapted from BControl (SidesPlotSection.m) 
% Kachi O. 2014.Mar.17
% Josh S. 2015.Jan.24 - optimized for speed

%% Code Starts Here
global nTrialsToShow %this is for convenience
global BpodSystem

switch Action
    case 'init'
        %initialize pokes plot
        SideList = varargin{1};
        
        nTrialsToShow = 200; %default number of trials to display
        
        if nargin > 3 %custom number of trials
            nTrialsToShow =varargin{2};
        end
        axes(AxesHandle);
        %plot in specified axes
        Xdata = 1:nTrialsToShow; Ydata = SideList(Xdata);
        BpodSystem.GUIHandles.FutureTrialLine = line([Xdata,Xdata],[Ydata,Ydata],'LineStyle','none','Marker','o','MarkerEdge','b','MarkerFace','b', 'MarkerSize',6);
        BpodSystem.GUIHandles.CurrentTrialCircle = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','k','MarkerFace',[1 1 1], 'MarkerSize',6);
        BpodSystem.GUIHandles.CurrentTrialCross = line([0,0],[0,0], 'LineStyle','none','Marker','+','MarkerEdge','k','MarkerFace',[1 1 1], 'MarkerSize',6);
        BpodSystem.GUIHandles.UnpunishedErrorLine = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','r','MarkerFace',[1 1 1], 'MarkerSize',6);
        BpodSystem.GUIHandles.PunishedErrorLine = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','r','MarkerFace','r', 'MarkerSize',6);
        BpodSystem.GUIHandles.RewardedCorrectLine = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','g','MarkerFace','g', 'MarkerSize',6);
        BpodSystem.GUIHandles.UnrewardedCorrectLine = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','g','MarkerFace',[1 1 1], 'MarkerSize',6);
        BpodSystem.GUIHandles.NoResponseLine = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','b','MarkerFace',[1 1 1], 'MarkerSize',6);
        BpodSystem.GUIHandles.EscapedLine = line([0,0],[0,0], 'LineStyle','none','Marker','o','MarkerEdge','c','MarkerFace','c', 'MarkerSize',6);
        set(AxesHandle,'TickDir', 'out','YLim', [-0.5, 1.5], 'YTick', [0 1],'YTickLabel', {'Go','NoGo'}, 'FontSize', 16);
        xlabel(AxesHandle, 'Trial#', 'FontSize', 18);
        hold(AxesHandle, 'on');
        
    case 'update'
        CurrentTrial = varargin{1};
        SideList = varargin{2};
        OutcomeRecord = varargin{3};
        
        if CurrentTrial<1
            CurrentTrial = 1;
        end
        
        % recompute xlim
%         [mn, mx] = rescaleX(AxesHandle,CurrentTrial,nTrialsToShow);
        mx = nTrialsToShow;
        mn = 1;
        
        %axes(AxesHandle); %cla;
        %plot future trials
        FutureTrialsIndx = CurrentTrial:mx;
        Xdata = FutureTrialsIndx; Ydata = SideList(Xdata);
        set(BpodSystem.GUIHandles.FutureTrialLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
        %Plot current trial
        set(BpodSystem.GUIHandles.CurrentTrialCircle, 'xdata', [CurrentTrial,CurrentTrial], 'ydata', [SideList(CurrentTrial-1),SideList(CurrentTrial-1)]);
        set(BpodSystem.GUIHandles.CurrentTrialCross, 'xdata', [CurrentTrial,CurrentTrial], 'ydata', [SideList(CurrentTrial-1),SideList(CurrentTrial-1)]);
        
        %Plot past trials
        if ~isempty(OutcomeRecord)
            indxToPlot = mn:CurrentTrial-1;
            %Plot Error, unpunished
            EarlyWithdrawalTrialsIndx =(OutcomeRecord(indxToPlot) == -1);
            Xdata = indxToPlot(EarlyWithdrawalTrialsIndx); Ydata = SideList(Xdata);
            set(BpodSystem.GUIHandles.UnpunishedErrorLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
            %Plot Error, punished
            InCorrectTrialsIndx = (OutcomeRecord(indxToPlot) == 0);
            Xdata = indxToPlot(InCorrectTrialsIndx); Ydata = SideList(Xdata);
            set(BpodSystem.GUIHandles.PunishedErrorLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
            %Plot Correct, rewarded
            CorrectTrialsIndx = (OutcomeRecord(indxToPlot) == 1);
            Xdata = indxToPlot(CorrectTrialsIndx); Ydata = SideList(Xdata);
            set(BpodSystem.GUIHandles.RewardedCorrectLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
            %Plot Correct, unrewarded
            UnrewardedTrialsIndx = (OutcomeRecord(indxToPlot) == 2);
            Xdata = indxToPlot(UnrewardedTrialsIndx); Ydata = SideList(Xdata);
            set(BpodSystem.GUIHandles.UnrewardedCorrectLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
            %Plot DidNotChoose
            DidNotChooseTrialsIndx = (OutcomeRecord(indxToPlot) == 3);
            Xdata = indxToPlot(DidNotChooseTrialsIndx); Ydata = SideList(Xdata);
            set(BpodSystem.GUIHandles.NoResponseLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
            
            EscapedTrialsIndx = (OutcomeRecord(indxToPlot) == 4);
            Xdata = indxToPlot(EscapedTrialsIndx); Ydata = SideList(Xdata);
            set(BpodSystem.GUIHandles.EscapedLine, 'xdata', [Xdata,Xdata], 'ydata', [Ydata,Ydata]);
        end
end

end

% function [mn,mx] = rescaleX(AxesHandle,CurrentTrial,nTrialsToShow)
% FractionWindowStickpoint = 0.75; % After this fraction of visible trials, the trial position in the window "sticks" and the window begins to slide through trials.
% mn = max(round(CurrentTrial - FractionWindowStickpoint*nTrialsToShow),1);
% mx = mn + nTrialsToShow-1;
% set(AxesHandle,'XLim',[mn-1 mx+1]);
% end



%% function OutcomePlot(AxesHandle,TrialTypeSides, OutcomeRecord, CurrentTrial)

function LickPlot(AxesHandle, Action, varargin)
%% 
% Plug in to Plot Go/No-Go lick outcome.
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
global TimeToShow %this is for convenience

switch Action
    case 'init'
        %initialize pokes plot      
        TimeToShow = 20; %default time window to display        
        if nargin > 3 %custom time window to display 
            TimeToShow =varargin{3};
        end
        axes(AxesHandle);
        %plot in specified axes
        xlim([0,TimeToShow]);
        set(AxesHandle,'TickDir', 'out','YLim', [-1, 2], 'YTick', [0 1],'YTickLabel', {'NoLick','Lick'}, 'FontSize', 16);
        xlabel(AxesHandle, 'Time (s)', 'FontSize', 18);
        hold(AxesHandle, 'on');
        
    case 'update'
        CurrentTrial = varargin{1};
        LickEvent = varargin{2};
        
        if CurrentTrial<1
            CurrentTrial = 1;
        end
        
        % axes(AxesHandle); %cla;
        %plot future trials
        cla(AxesHandle)
        if ~isempty(LickEvent)
            stem([LickEvent(1,:),LickEvent(2,:)],ones(1,2*size(LickEvent,2)),'Marker','none','Color','k')
            hold on
            line(LickEvent,ones(size(LickEvent)),'Color','k')
            hold off
            set(AxesHandle,'TickDir', 'out','YLim', [-1, 2], 'YTick', [0 1],'YTickLabel', {'NoLick','Lick'}, 'FontSize', 16);
        xlabel(AxesHandle, 'Time (s)', 'FontSize', 18);
        end
               
end

end

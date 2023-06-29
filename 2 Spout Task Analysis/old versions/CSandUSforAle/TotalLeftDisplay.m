%{
----------------------------------------------------------------------------

This file is part of the Bpod Project
Copyright (C) 2014 Joshua I. Sanders, Cold Spring Harbor Laboratory, NY, USA

----------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed  WITHOUT ANY WARRANTY and without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
function TotalLeftDisplay(varargin)
% TotalRewardDisplay('init') - initializes a window that displays total reward
% TotalRewardDisplay('add', Amount) - updates the total reward display with
% a new reward, adding to the total amount.

global BpodSystem
Op = varargin{1};
if nargin > 1
    AmountToAdd = varargin{2};
end
Op = lower(Op);
switch Op
    case 'init'
        BpodSystem.PluginObjects.TotalLeftDelivered = 0;
        BpodSystem.ProtocolFigures.TotalLeftDisplay = figure('Position', [450 550 150 150],'name','Left','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off', 'Color', [.6 .6 1]);
        BpodSystem.GUIHandles.TotalLeftDisplay = struct;
        BpodSystem.GUIHandles.TotalLeftDisplay.Label = uicontrol('Style', 'text', 'String', 'Left', 'units', 'normalized', 'Position', [.15 .7 .7 .15], 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Arial', 'BackgroundColor', [.7 .7 1]);
        BpodSystem.GUIHandles.TotalLeftDisplay.Amount = uicontrol('Style', 'text', 'String', ['0 ' char(181) 'l'], 'units', 'normalized', 'Position', [.1 .25 .8 .25], 'FontSize', 24, 'FontName', 'Arial', 'BackgroundColor', [.7 .7 1]);
        
    case 'add'
        BpodSystem.PluginObjects.TotalLeftDelivered = BpodSystem.PluginObjects.TotalLeftDelivered + AmountToAdd;
        if BpodSystem.PluginObjects.TotalLeftDelivered > 1000
            DisplayAmount = BpodSystem.PluginObjects.TotalLeftDelivered/1000;
            DisplayUnits = ' ml';
        else
            DisplayAmount = BpodSystem.PluginObjects.TotalLeftDelivered;
            DisplayUnits = [' ' char(181) 'l'];
        end
        set(BpodSystem.GUIHandles.TotalLeftDisplay.Amount, 'String', [num2str(DisplayAmount) DisplayUnits]);
end

%3 receives reward with assist. 1 receive reward without assist. 0.
%punishment. -1, avoided air puff. 2 for ignored reward. 

function [outcomes] = fix_outcomes(SessionData)
    %outcomes = SessionData.Outcomes;
    outcomes = [];
    spout = 1;
    if isfield(SessionData, 'TrialTypes')
        if ~isempty(SessionData.TrialTypes)
            trial_types = SessionData.TrialTypes;
        else
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = trial_types;
        end
    else
        
        if spout == 3
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = trial_types;
        elseif spout == 2
            trial_types = ones(1, SessionData.nTrials);
        elseif spout == 1
            trial_types = zeros(1, SessionData.nTrials);
        end
    end
    
    for x = 1:length(trial_types)
        if trial_types(x) == 1 %If the trial is an air trial
            if isfield(SessionData.RawEvents.Trial{1,x}.States, 'DeliverPunishment')
                punish_start = SessionData.RawEvents.Trial{1,x}.States.DeliverPunishment(1);
                punish_end = SessionData.RawEvents.Trial{1,x}.States.DeliverPunishment(2);

                if ~isnan(punish_start)            %If not too many USs
                    punish_length = punish_end - punish_start;
                    punish_amount = SessionData.TrialSettings(x).GUI.PunishAmount;
                    if punish_length < punish_amount
                        %Escaped
                        outcomes = horzcat(outcomes, 4);
                    else
                        %Got full air puff
                        outcomes = horzcat(outcomes, 0);
                    end
                else
                    %Avoids the air puff correctly
                    outcomes = horzcat(outcomes, -1);
                end
            else
                outcomes = horzcat(outcomes, -1);
            end
            
        elseif trial_types(x) == 0 %If it's a water trial
            if isfield(SessionData.RawEvents.Trial{1,x}.States, 'DeliverReward')
                reward_time = SessionData.RawEvents.Trial{1,x}.States.DeliverReward(1);
                response_s = SessionData.RawEvents.Trial{1,x}.States.ResponseW(1);
                response_e = SessionData.RawEvents.Trial{1,x}.States.ResponseW(2);
                window = response_e - response_s;
                if ~isnan(reward_time)            %If not too many USs
                    if window < 2.9995
                        outcomes = horzcat(outcomes, 1);
                    else
                        outcomes = horzcat(outcomes, 3);
                    end
                else
                    outcomes = horzcat(outcomes, 2);
                end
            else
                outcomes = horzcat(outcomes, 1);
            end
            
            
        elseif trial_types(x) == 2 %If it's a laser trial
            outcomes = horzcat(outcomes, 5);
        end
    end
end

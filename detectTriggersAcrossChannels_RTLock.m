function [trl, event] = detectTriggersAcrossChannels_RTLock(cfg)
% read events in data, locked to RT.

if isfield(cfg, 'dataset')
    data = cfg.dataset;
elseif isfield(cfg, 'datafile')
    data = cfg.datafile;
else
    error('no recognized data file')
end

eventOrig = ft_read_event(data,'trigindx', cfg.trialdef.trigChannels,'threshold',3,'detectflank','down');

event = eventOrig;
for e = 1:length(event)
    event(e).sample = event(e).sample + round(cfg.RT(e));
end

pretrig = -1*cfg.trialdef.prestim;
posttrig = cfg.trialdef.poststim;

sample_h = cell(length(cfg.trialdef.trigChannels),1);
for i=1:length(cfg.trialdef.trigChannels)
    trigChannelStr = num2str(cfg.trialdef.trigChannels(i));
    
    % search for "trigger" events according to 'trigchannel' defined outside the function
    sample_h{i} = [event(strcmp(trigChannelStr, {event.type})).sample]';
end

% sortsample by chronological order
sample = sort(vertcat(sample_h{:}));

trl = nan(length(event),4);
% creating your own trialdefinition based upon the events
for j = 1:length(sample);
    trlbegin = sample(j) + pretrig;
    trlend   = sample(j) + posttrig;
    offset   = pretrig;
    newtrl   = [ trlbegin trlend offset j];
    trl(j,:) = newtrl;
end

end
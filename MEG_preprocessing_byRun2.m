function res = MEG_preprocessing_byRun2(par, flags)
%%   MEG_preprocessing_byRun
%   Preprocess runs of MEG data using mainly fieldtrip routines
%   Alan Gordon, Stanford University, 08/19/2013
%   Adapted from scripts by Andy Heusser, NYU

if (nargin<2)
    flags = 'cf';
end

ft_defaults

%% read header
res.hdrInit = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

%res.deadChannels = find(checkForDeadChannels(par.dataFiles{1}));

% for each run
for f = 1:length(par.dataRuns)
    thisRun = par.dataRuns{f};
    [~, thisRunName] = fileparts(thisRun);
    
    %% Read in continuous data
    if ismember('c', flags)
        % read in continuous data
        cfg = [];
        cfg.dataset = thisRun;
        cfg.trialdef.eventtype  = '?';
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials     = Inf;
        cfg = ft_definetrial(cfg);
        res.continuous_data = ft_preprocessing(cfg);
        res.continuous_data.label = res.hdrInit.label; %for some reason, the run-specific .sqds don't contain labels fields...
    end
    
    %% filter timeseries
    if ismember('f', flags)
        %high-pass filter
        cfg = [];
        cfg.hpfilter = par.doHPFilt;
        cfg.hpfreq = par.HPFreq;
        cfg.channel = par.MEGChannels;
        res.continuous_data_hp = ft_preprocessing(cfg,res.continuous_data);
        
        %low-pass filter
        cfg = [];
        cfg.lpfilter = par.doLPFilt;
        cfg.lpfreq = par.LPFreq;
        cfg.channel = par.MEGChannels;
        res.continuous_data_bpf = ft_preprocessing(cfg,res.continuous_data_hp);
    end
    
    
    %% get events
    cfg = [];
    
    
    cfg.dataset = thisRun;
    cfg.trialdef.trigChannels = par.trigChannels;
    cfg.trialfun = par.trialFun;
    
    cfg.continuous = par.continuous;
    cfg.trialdef.prestim = 200;
    cfg.trialdef.poststim = 2000;
    
    cfg = ft_definetrial(cfg);
    res.event = cfg.event;        
    
    %% save data and events
    res.continuous_data_bpf_all = res.continuous_data;
    res.continuous_data_bpf_all.trial{1}(par.idxMEGChan,:) = res.continuous_data_bpf.trial{1};
    
    thisFile = fullfile(par.preprocRunsDir,[thisRunName, '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '_all_filt.mat']);
    
    ft_write_data(thisFile, res.continuous_data_bpf_all.trial{1}, 'header', res.hdrInit, 'dataformat', 'fcdc_matbin');
    ft_write_event(thisFile, res.event)
    
    %% bin data into trials
    cfg = [];
    cfg.channel = par.MEGChannels;
    cfg.dataset = thisFile;
    cfg.trialdef.trigChannels = par.trigChannels;
    cfg.trialfun = par.trialFun;
    
    cfg.continuous = par.continuous;
    cfg.trialdef.prestim = par.preStimOrigSamples;
    cfg.trialdef.poststim = par.postStimOrigSamples;
    
    cfg = ft_definetrial(cfg);
    res.data_epochs = ft_preprocessing(cfg);
    
    %% resample data
    cfg = [];
    cfg.detrend = 'no';
    cfg.resamplefs = par.resampleRate;
    res.data_epochs = ft_resampledata(cfg, res.data_epochs);
    

    %% discard unwanted data
    res = rmfield(res, 'continuous_data');
    res = rmfield(res, 'continuous_data_hp');
    res = rmfield(res, 'continuous_data_bpf');
    res = rmfield(res, 'continuous_data_bpf_all');
    
    %% fix header
    res.hdr = res.hdrInit;
    res.hdr.Fs = res.data_epochs.fsample;
    res.hdr.label = res.data_epochs.label;
    res.hdr.chantype = res.hdr.chantype(par.idxMEGChan);
    res.hdr.chanunit = res.hdr.chanunit(par.idxMEGChan);
    res.hdr.nSamples = size(res.data_epochs.trial{1},2);
    res.hdr.nChans = size(res.data_epochs.trial{1},1);
        
    %% save data    
    thisFile = fullfile(par.preprocRunsDir,[thisRunName, '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '_filt_epochs.mat']);
    save(thisFile, 'res')
end
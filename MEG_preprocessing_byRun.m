function res = MEG_preprocessing_byRun(par, flags)
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
    
    %% resample data
    cfg = [];
    cfg.detrend = 'no';
    cfg.resamplefs = par.resampleRate;
    res.continuous_data_bpf_rs = ft_resampledata(cfg, res.continuous_data_bpf);
    
    %% get onset info
    cfg = [];
    cfg.dataset = thisRun;
    cfg.trialdef.trigChannels = par.trigChannels;
    cfg.trialfun = par.trialFun;
    
    cfg.continuous = par.continuous;
    cfg.trialdef.prestim = 0;
    cfg.trialdef.poststim = 0;
    
    cfg = ft_definetrial(cfg);
    res.event = cfg.event;        
    
    %% discard unwanted data
    res = rmfield(res, 'continuous_data');
    res = rmfield(res, 'continuous_data_hp');
    res = rmfield(res, 'continuous_data_bpf');
    
    %% fix header
    res.hdr = res.hdrInit;
    res.hdr.Fs = res.continuous_data_bpf_rs.fsample;
    res.hdr.label = res.continuous_data_bpf_rs.label;
    res.hdr.chantype = res.hdr.chantype(par.idxMEGChan);
    res.hdr.chanunit = res.hdr.chanunit(par.idxMEGChan);
    res.hdr.nSamples = size(res.continuous_data_bpf_rs.trial{1},2);
    res.hdr.nChans = size(res.continuous_data_bpf_rs.trial{1},1);
        
    %% save data
    %save (fullfile(par.preprocRunsDir, [thisRunName, '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '_filt.mat']), 'res');
    
    thisFile = fullfile(par.preprocRunsDir,[thisRunName, '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '_filt.mat']);
    
    ft_write_data(thisFile, res.continuous_data_bpf_rs.trial{1}, 'header', res.hdr, 'dataformat', 'fcdc_matbin');
    ft_write_event(thisFile, res.event)
end
function res = MEG_preprocessing_FreqDecomp(par)
%%   MEG_preprocessing_byRun
%   Preprocess runs of MEG data using mainly fieldtrip routines
%   Alan Gordon, Stanford University, 08/19/2013
%   Adapted from scripts by Andy Heusser, NYU


ft_defaults

%% read header
res.hdrInit = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

idx = MEG_behav_analysis(par);
%res.deadChannels = find(checkForDeadChannels(par.dataFiles{1}));

% for each run
for f = 1:length(par.dataRuns)
    thisRun = par.dataRuns{f};
    [~, thisRunName] = fileparts(thisRun);
    
    %% Read in continuous data

    % read in continuous data
    cfg = [];
    cfg.dataset = thisRun;
    cfg.trialdef.eventtype  = '?';
    cfg.trialdef.triallength = Inf;
    cfg.trialdef.ntrials     = Inf;
    cfg = ft_definetrial(cfg);
    res.continuous_data = ft_preprocessing(cfg);
    res.continuous_data.label = res.hdrInit.label; %for some reason, the run-specific .sqds don't contain labels fields...
 
    %% get events
    cfg = [];
       
    cfg.dataset = thisRun;
    cfg.trialdef.trigChannels = par.trigChannels;
    cfg.trialfun = par.trialFun;
    
    cfg.continuous = par.continuous;
    cfg.trialdef.prestim = par.preStimOrigSamples;
    cfg.trialdef.poststim = par.postStimOrigSamples;
    cfg.RT = par.acquisitionRate * idx.test.RT(idx.test.miniblock==f);
    
    cfg = ft_definetrial(cfg);
    res.event = cfg.event;
    
    %% save data and events
    
    thisFile = fullfile(par.preprocRunsDir,[thisRunName, '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '_no_filt_stimLockedEvents.mat']);
    ft_write_data(thisFile, res.continuous_data.trial{1}, 'header', res.hdrInit, 'dataformat', 'fcdc_matbin');
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
    cfg.RT = par.acquisitionRate * idx.test.RT(idx.test.miniblock==f);
    
    cfg = ft_definetrial(cfg);
    res.data_epochs = ft_preprocessing(cfg);    
    
    
    %% spectral decomposition
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'MEG';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = par.foi;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = par.t_ftimwin;   % length of time window = 0.5 sec
    cfg.toi          = (-1*par.preStimFreqEpochsInSec):par.FreqEpochWindow:par.postStimFreqEpochsInSec; % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
    res.continuous_data_freq = ft_freqanalysis(cfg, res.data_epochs);

    
    %% resample data
%     cfg = [];
%     cfg.detrend = 'no';
%     cfg.resamplefs = par.resampleRate;
%     res.continuous_data_freq_rs = ft_resampledata(cfg, res.continuous_data_freq);
    

    %% plot spectral decomposition
%     cfg = [];
%     cfg.layout = par.dataFiles{1};
%     res.layout = ft_prepare_layout(cfg);
%     
%     cfg = [];
%     cfg.baseline     = [-0.5 -0.1];
%     cfg.baselinetype = 'absolute';
%     cfg.zlim         = [-8e-27 8e-27];
%     cfg.showlabels   = 'yes';
%     cfg.layout       = res.layout;
%     ft_multiplotTFR(cfg, res.continuous_data_freq);
    
    %% fix header
    res.hdr = res.hdrInit;
    res.hdr.Fs = res.data_epochs.fsample;
    res.hdr.label = res.data_epochs.label;
    res.hdr.chantype = res.hdr.chantype(par.idxMEGChan);
    res.hdr.chanunit = res.hdr.chanunit(par.idxMEGChan);
    res.hdr.nSamples = size(res.data_epochs.trial{1},2);
    res.hdr.nChans = size(res.data_epochs.trial{1},1);
        
   %% discard unwanted data
    res = rmfield(res, 'continuous_data');
    res = rmfield(res, 'data_epochs');
    
    %% save data    
    thisFile = fullfile(par.preprocRunsDir,[thisRunName, '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) ...
        '_filt_epochs_' num2str(par.preStimOrigSamples) '_to_' num2str(par.postStimOrigSamples) '_' par.lock 'Lock_spectral' '.mat']);
    save(thisFile, 'res')
end
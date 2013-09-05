function [ output_args ] = MEG_preprocess_acrossRun(par, flags)

%% MEG sample data preprocessing
%   Preprocess MEG using mainly fieldtrip routines
%   Alan Gordon, Stanford University, 09/02/2013
%   Adapted from scripts by Andy Heusser, NYU

if (nargin<2)
    flags = 'tai';
end

%concatData = load(par.dataConcat);

ft_defaults

res.hdr = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

cfg = [];
cfg.layout = par.dataFiles{1};
res.layout = ft_prepare_layout(cfg);
%% organize data into trials

if ismember('t', flags)
    
    cfg = [];
    %cfg.dataset
    cfg.datafile     = par.dataConcat;
    cfg.headerfile   = par.dataConcat;
    cfg.trialdef.trigChannels = par.trigChannels;
    cfg.trialfun = par.trialFun;
    
    cfg.continuous = par.continuous;
    cfg.trialdef.prestim = par.trialDefPreStim;
    cfg.trialdef.poststim = par.trialDefPostStim;
    
    % enter trl in cfg
    cfg = ft_definetrial(cfg);
    
    % read data
    % note that 1 second prestimulus is included in each trial of data, due
    % to processing in MEG_splitSqdIntoRuns
    res.data_epochs = ft_preprocessing(cfg);
    
    preprocData = res.data_epochs;
end



cfgOrig = cfg;

%% identify jump artifacts
if ismember('a', flags)
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel    = 'MEG';
    cfg.artfctdef.zvalue.cutoff     = 20;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative    = 'yes';
    cfg.artfctdef.zvalue.medianfilter  = 'yes';
    cfg.artfctdef.zvalue.medianfiltord = 9;
    cfg.artfctdef.zvalue.absdiff       = 'yes';
    
    % make the process interactive
    cfg.artfctdef.zvalue.interactive = 'yes';
    
    [~, art.jump] = ft_artifact_zvalue(cfg);
    
    %% identify muscle artifacts
    
    cfg = cfgOrig;
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel = 'MEG';
    cfg.artfctdef.zvalue.cutoff      = 4;
    cfg.artfctdef.zvalue.trlpadding  = 0;
    cfg.artfctdef.zvalue.fltpadding  = 0;
    cfg.artfctdef.zvalue.artpadding  = 0.1;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter    = 'yes';
    cfg.artfctdef.zvalue.bpfreq      = [60 90];
    cfg.artfctdef.zvalue.bpfiltord   = 9;
    cfg.artfctdef.zvalue.bpfilttype  = 'but';
    cfg.artfctdef.zvalue.hilbert     = 'yes';
    cfg.artfctdef.zvalue.boxcar      = 0.2;
    
    % make the process interactive
    cfg.artfctdef.zvalue.interactive = 'yes';
    
    [~, art.muscle] = ft_artifact_zvalue(cfg);
    
    %% manually identify bad trials
    
    cfg          = [];
    cfg.method   = 'summary';
    cfg.alim     = 1e-12;
    cfg.megscale = 1;
    cfg.eogscale = 5e-8;
    dummy        = ft_rejectvisual(cfg,res.data_epochs);
    
    rejectVisual_h = setdiff(1:length(cfgOrig.trl), dummy.trialinfo);
    art.rejectVisual = cfgOrig.trl(rejectVisual_h,1:2);
    
    clear dummy;
    
    %% reject all artifacts
    cfg = cfgOrig;
    cfg.artfctdef.zvalue.artifact = vertcat(art.jump, art.muscle, art.rejectVisual);
    
    [res.data_epochs_ar] = ft_rejectartifact(cfg, res.data_epochs);
    
    preprocData = res.data_epochs_ar;
    %idx.goodTrials = ismember(cfgOrig.trl(:,1), res.data_epochs_ar.cfg.trl(:,1));
end

%% calculate components
if ismember('o', flags)
    cfg        = [];
    cfg.channel = 'MEG';
    cfg.numcomponent = par.NComponents;
    res.components = ft_componentanalysis(cfg, res.data_epochs_ar);
    
    %% plot scalp maps of components
    cfg             = [];
    cfg.grad        = res.grad;
    lay             = ft_prepare_layout(cfg, res.components);
    cfg.layout      = lay;
    cfg.component   = par.componentsToPlot;
    cfg.comment     = 'no';
    ft_topoplotIC(cfg, res.components)
    
    %% plot timeseries from bad components
    cfg = [];
    cfg.channel   = par.artComponents;
    cfg.grad = res.grad;
    cfg.layout = res.layout;
    cfg.viewmode ='component';
    cfg.zlim =  'maxabs';
    cfg.compscale =  'local';
    ft_databrowser(cfg, res.components)
    
    %% reject components
    cfg = [];
    cfg.component = par.componentsToRemove;
    cfg.channel = {'all'};
    res.data_epochs_ar_postica = ft_rejectcomponent(cfg, res.components);
    preprocData = res.data_epochs_ar_postica;
end

%% save data
thisFile = fullfile(par.concatRunsDir,['MEGDatFiltRSConcat', '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '_artRemoved.mat']);
save(thisFile, 'preprocData')

end


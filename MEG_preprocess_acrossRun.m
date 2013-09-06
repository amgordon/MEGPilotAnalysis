function [ output_args ] = MEG_preprocess_acrossRun(par, flags)

%% MEG sample data preprocessing
%   Preprocess MEG using mainly fieldtrip routines
%   Alan Gordon, Stanford University, 09/02/2013
%   Adapted from scripts by Andy Heusser, NYU

if (nargin<2)
    flags = 'tami';
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


    
%% manually identify bad trials
 
if ismember('m', flags)
    cfg          = [];
    cfg.method   = 'summary';
    cfg.alim     = 1e-12;
    cfg.megscale = 1;
    cfg.eogscale = 5e-8;
    preprocData        = ft_rejectvisual(cfg,res.data_epochs);
    
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


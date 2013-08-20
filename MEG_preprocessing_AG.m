%% MEG sample data preprocessing
%   Preprocess MEG using mainly fieldtrip routines
%   Alan Gordon, Stanford University, 08/19/2013
%   Adapted from scripts by Andy Heusser, NYU


%%  parameters

par.exptdir = '/biac4/wagner/biac3/wagner5/alan/MEG';
par.scriptDir = fullfile(par.exptdir, 'scripts');
par.dataDir = fullfile(par.exptdir, 'sampleData');

par.dataName = fullfile(par.dataDir, 'R0485_Checkerboard_7.9.13.sqd');

par.fileformat = 'besa_pos';
par.neighborCalcMethod = 'triangulation';
par.neighborfeedback = 'no';

par.tscpa_shifts = [-100:100];
par.refChannels = {'RM158' 'RM159' ' RM160'};

par.MEGChannels = {'MEG'};
par.LPFreq = 101;
par.HPFreq = 1;
par.doLPFilt = true;
par.doHPFilt = true;

par.trigChannel = '161';
par.trialFun = 'detectTrigger';
par.continuous = 'yes';
par.trialDefPreStim = 1;
par.trialDefPostStim = 1.5;

par.MEGDataLabels = hdr.label(1:157);

par.NComponents = 60;
par.componentsToPlot = 1:25;

par.artComponents = [1, 4, 20];

par.componentsToRemove = [1,4,20];

par.fontSize = 8;
%% read in data and load scripts into path
addpath(genpath(par.scriptDir))
addpath(genpath(par.dataDir));

ft_defaults

%% Map out spatial layout of sensors

% define header and gradiometer
hdr = ft_read_header(par.dataName);
grad = ft_read_sens(par.dataName,'filename',par.fileformat);

% define layout
cfg = [];
cfg.layout = par.dataName;
res.layout = ft_prepare_layout(cfg);

% calculate neigbours
cfg          = [];
cfg.feedback = par.neighborfeedback;
cfg.layout = grad;
cfg.method = par.neighborCalcMethod;
res.neighbours = ft_prepare_neighbours(cfg,hdr);

%% Read in continuous data

% read in continuous data
cfg = [];
cfg.dataset = dataName;
cfg.trialdef.eventtype  = '?';
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials     = Inf;

% enter trl in cfg
cfg = ft_definetrial(cfg);
res.continuous_data = ft_preprocessing(cfg);

%% replace dead channels with interpolated data from nearby channels.

res.deadChannels = find(checkForDeadChannels(dataName));

cfg = [];
cfg.badchannel = res.continuous_data.label(res.deadChannels);
cfg.neighbours = res.neighbours;
res.continuous_data_rc = ft_channelrepair(cfg,res.continuous_data);

%% denoise sensors with tspca
cfg = [];
cfg.shifts = par.tscpa_shifts;
cfg.refchannel = par.refChannels;
cfg.channel = par.MEGChannels;
res.continuous_data_rc_tspca = denoise_tsr(cfg,continuous_data_rc);

%% filter timeseries

cfg = [];
cfg.hpfilter = par.doHPFilt;
cfg.hpfreq = par.HPFreq;
cfg.lpfilter = par.doLPFilt;
cfg.lpfreq = par.LPFreq;
par.filterChannels = par.MEGChannels;
res.continuous_data_rc_tspca_bpf = ft_preprocessing(cfg,res.continuous_data_rc_tspca);

%% organize data into trials
cfg = [];
cfg.dataset = dataName;
cfg.trialdef.trigchannel = par.trigChannel;
cfg.trialfun = par.trialFun;

cfg.continuous = par.continuous;
cfg.trialdef.prestim = par.trialDefPreStim;
cfg.trialdef.poststim = par.trialDefPostStim;

% enter trl in cfg
cfg = ft_definetrial(cfg);

% read data
res.data_epochs = ft_preprocessing(cfg);

%% identify and reject bad trials
cfg          = [];
cfg.method   = 'trial';
cfg.channel = par.MEGDataLabels;
cfg.alim     = 1e-12;
data_epochs_tspca_hr  = ft_rejectvisual(cfg,res.data_epochs);


%% calculate components
cfg        = [];
cfg.channel = par.MEGChannels;
cfg.numcomponent = par.NComponents;
res.components = ft_componentanalysis(cfg, data_epochs_tspca_hr);

%% plot scalp maps of components
cfg             = [];
cfg.grad        = grad;
lay             = ft_prepare_layout(cfg, res.components);
cfg.layout      = lay;
cfg.component   = par.componentsToPlot;
cfg.comment     = 'no';
ft_topoplotIC(cfg, res.components)

%% plot timeseries from bad components
cfg = [];
cfg.channel   = par.artComponents;
cfg.grad = grad;
cfg.layout = layout;
cfg.viewmode ='component';
cfg.zlim =  'maxabs';
cfg.compscale =  'local';
ft_databrowser(cfg, res.components)

%% reject components
cfg = [];
cfg.component = par.componentsToRemove;
cfg.channel = {'all'};
res.data_epochs_postica_tspca = ft_rejectcomponent(cfg, res.components);


%% save data
save preprocessedData res par
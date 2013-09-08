function par = MEG_par(substr)

%% reading in input
if isnumeric(substr) 
    substr = MEG_num2Sub(substr);
end
par.substr = substr;

%% which sub-task: 'Study' or 'Test'?
par.subTask = 'Study';

%% directory info
par.exptdir = '/biac4/wagner/biac3/wagner5/alan/MEG';
par.scriptDir = fullfile(par.exptdir, 'scripts/MEGPilotAnalysis');
par.dataDir = fullfile(par.exptdir, 'PilotData');
par.subDir = fullfile(par.dataDir, par.substr);
par.rawDir = fullfile(par.subDir, 'raw');
par.runsDir = fullfile(par.subDir, 'runs');
par.preprocRunsDir = fullfile(par.subDir, 'preprocRuns'); 
par.concatRunsDir = fullfile(par.subDir, 'concatRuns'); 
par.behavDir = fullfile(par.subDir, 'behav');
par.headerDir = fullfile(par.subDir, 'header');

%% subject-specific parameters

switch substr
    case 'meg_081913'
        par.dataFiles{1} = fullfile(par.rawDir, 'R0485_MEGClass_8.19.13.sqd');
        par.nCycles = 25;
        par.componentsToPlot = 1:25;
        par.artComponents = [1,4,20];        
        par.componentsToRemove = [];
        par.badTrials = [164 368 437 583 877 878 879 880];
    case 'meg_082013'
        par.dataFiles{1} = fullfile(par.rawDir, 'R0504_MEGclass_8.20.13.sqd');
        par.nCycles = 20;
        par.componentsToPlot = 1:25;
        par.artComponents = [1,4,20];        
        par.componentsToRemove = [];
        par.badTrials = [];
    case 'meg_082113'
        par.dataFiles{1} = fullfile(par.rawDir, 'R0487_MEGclassblock1_8.21.13.sqd');
        par.dataFiles{2} = fullfile(par.rawDir, 'R0487_MEGclassblock2_8.21.13.sqd');
        par.dataFiles{3} = fullfile(par.rawDir, 'R0487_MEGclassblock3_8.21.13.sqd');
        par.componentsToPlot = 1:25;
        par.artComponents = [1,4,20];
        par.componentsToRemove = [];
        par.nCycles = 25;
        par.badTrials = [];
    otherwise
        error('unrecognized subject');
end


%% data runs
dr = dir(fullfile(par.runsDir, 'MEGDat*sqd'));
drp = dir(fullfile(par.preprocRunsDir, 'MEGDat*HP1*LP30*epochs.mat'));
drpc = dir(fullfile(par.concatRunsDir, 'MEGDat*.mat'));
drpca = dir(fullfile(par.concatRunsDir, 'MEGDat*artRemoved.mat'));

b = dir(fullfile(par.behavDir, 'MEG*study.mat'));

for i=1:length(dr)
   par.dataRuns{i} = fullfile(par.runsDir, dr(i).name); 
end

for i=1:length(drp)
   par.dataRunsPreproc{i} = fullfile(par.preprocRunsDir, drp(i).name); 
end

par.dataConcat = fullfile(par.concatRunsDir, drpc(1).name);
par.headerFile = fullfile(par.headerDir, 'hdr.mat');
if length(drpca)>0
    par.dataConcatNoArt = fullfile(par.concatRunsDir, drpca(1).name);
end
par.behavFile = fullfile(par.behavDir, b(1).name);

%%
par.miniBlockLength = 12;

%% splitting sqds into runs
par.nPretriggerSamplesForRun = 1000;

%% channel info
par.triggerRange = 160:166;
par.idxMEGChan = 1:157;
par.idxNonMEGChan = 158:192;
par.trigChannels = 160:166;
par.MEGChannels = {'MEG'};

%% frequency info
par.acquisitionRate = 1000; %in Hz
par.resampleRate = 200; % in Hz

%% filtering 
par.LPFreq = 30;
par.HPFreq = 1;
par.doLPFilt = 'yes';
par.doHPFilt = 'yes';

%% event-related analysis
% in samples where Fs = par.resampleRate

% also note that epochs have been prepadded with
% par.nPretriggerSamplesForRun during the MEG_splitSqdIntoRuns script.

par.preStimOrigSamples = 199;
par.postStimOrigSamples = 2000;
par.trialDefPreStim = 40; 
par.trialDefPostStim = 400;

%%
par.fileformat = 'besa_pos';
par.neighborCalcMethod = 'triangulation';
par.neighborfeedback = 'no';

par.tscpa_shifts = [-100:100];
par.refChannels = {'RM158' 'RM159' ' RM160'};



par.trialFun = 'detectTriggersAcrossChannels';
par.continuous = 'yes';


MEGDataLabels_h = load(fullfile(par.scriptDir, 'MEGLabels.mat'));
par.MEGDataLabels = MEGDataLabels_h.MEGLabels;

par.NComponents = 60;


par.fontSize = 8;
end


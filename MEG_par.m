function par = MEG_par(substr)
% creates subject-specific parameter information.

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
        par.badTrials = [164 368 437 583 648 877 878 879 880];
        par.subNo=1;
    case 'meg_082013'
        par.dataFiles{1} = fullfile(par.rawDir, 'R0504_MEGclass_8.20.13.sqd');
        par.nCycles = 20;
        par.componentsToPlot = 1:25;
        par.artComponents = [1,4,20];        
        par.componentsToRemove = [];
        par.badTrials = [];
        par.subNo=2;
    case 'meg_082113'
        par.dataFiles{1} = fullfile(par.rawDir, 'R0487_MEGclassblock1_8.21.13.sqd');
        par.dataFiles{2} = fullfile(par.rawDir, 'R0487_MEGclassblock2_8.21.13.sqd');
        par.dataFiles{3} = fullfile(par.rawDir, 'R0487_MEGclassblock3_8.21.13.sqd');
        par.componentsToPlot = 1:25;
        par.artComponents = [1,4,20];
        par.componentsToRemove = [];
        par.nCycles = 25;
        par.badTrials = [];
        par.subNo=3;
    otherwise
        error('unrecognized subject');
end


%% data runs
par.lock = 'stim';

dr = dir(fullfile(par.runsDir, 'MEGDat*sqd'));
drp = dir(fullfile(par.preprocRunsDir, 'MEGDat*HP1*LP30*epochs.mat'));
drps = dir(fullfile(par.preprocRunsDir, 'MEGDat*Study*HP1*LP30*epochs*200*5000*.mat'));
%drpt = dir(fullfile(par.preprocRunsDir, 'MEGDat*Test*HP1*LP30*epochs_RTLock.mat'));
drpt = dir(fullfile(par.preprocRunsDir, 'MEGDat*Test*HP1*LP30*epochs*200*5000*.mat'));
drpc = dir(fullfile(par.concatRunsDir, 'MEGDat*.mat'));
drpca = dir(fullfile(par.concatRunsDir, 'MEGDat*artRemoved.mat'));


drps_spec = dir(fullfile(par.preprocRunsDir, 'MEGDat*Study*HP1*LP30*700*spectral.mat'));
drpt_spec = dir(fullfile(par.preprocRunsDir, 'MEGDat_Test*HP1_LP30*epochs*700*stimLock_spectral.mat'));

b = dir(fullfile(par.behavDir, 'MEG*study.mat'));

for i=1:length(dr)
   par.dataRuns_h{i} = fullfile(par.runsDir, dr(i).name); 
end

for i=1:length(drp)
   par.dataRunsPreproc{i} = fullfile(par.preprocRunsDir, drp(i).name); 
end

for i=1:length(drps)
   par.dataRunsPreprocStudy{i} = fullfile(par.preprocRunsDir, drps(i).name); 
end

for i=1:length(drpt)
   par.dataRunsPreprocTest{i} = fullfile(par.preprocRunsDir, drpt(i).name); 
end

for i=1:length(drps_spec)
   par.dataRunsPreprocFreqStudy{i} = fullfile(par.preprocRunsDir, drps_spec(i).name); 
end

for i=1:length(drpt_spec)
   par.dataRunsPreprocFreqTest{i} = fullfile(par.preprocRunsDir, drpt_spec(i).name); 
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
par.occChannels = {'AG015' 'AG045' 'AG062' 'AG014' 'AG026' 'AG043' ...
'AG025' 'AG023' 'AG013' 'AG063' 'AG032' 'AG047' 'AG044' 'AG046' ...
'AG036' 'AG008' 'AG001' 'AG020' 'AG002' 'AG051' 'AG018' 'AG005' ...
'AG004' 'AG003' 'AG006' 'AG019' 'AG050'};

%% frequency info
par.acquisitionRate = 1000; %in Hz
par.resampleRate = 200; % in Hz

%% classification
par.classificationFs = 10;

%% filtering 
par.LPFreq = 30;
par.HPFreq = 1;
par.doLPFilt = 'yes';
par.doHPFilt = 'yes';

%% event-related analysis
% in samples where Fs = par.resampleRate
%
% note that epochs have been prepadded with
% par.nPretriggerSamplesForRun during the MEG_splitSqdIntoRuns script.

par.preStimOrigSamples = 200;
par.postStimOrigSamples = 5000;
par.preStimFreqEpochs = 200;
par.postStimFreqEpochs = 5000;

par.trialDefPreStim = 40; 
par.trialDefPostStim = 400;

par.preStimFreqEpochsInSec = par.preStimFreqEpochs/par.acquisitionRate;
par.postStimFreqEpochsInSec = par.postStimFreqEpochs/par.acquisitionRate;
par.FreqEpochWindow = .05;
par.foi = 2:2:30;
par.t_ftimwin    = ones(length(par.foi),1).*0.5; 

%% miscellaneous params
par.fileformat = 'besa_pos';
par.neighborCalcMethod = 'triangulation';
par.neighborfeedback = 'no';

par.tscpa_shifts = [-100:100];
par.refChannels = {'RM158' 'RM159' ' RM160'};

if strcmp(par.lock, 'stim')
    par.trialFun = 'detectTriggersAcrossChannels';
    par.dataRuns = par.dataRuns_h;
elseif strcmp(par.lock, 'RT')
    par.trialFun = 'detectTriggersAcrossChannels_RTLock';
    par.dataRuns = par.dataRuns_h((par.nCycles*2+1):end);
end
par.continuous = 'yes';

MEGDataLabels_h = load(fullfile(par.scriptDir, 'MEGLabels.mat'));
par.MEGDataLabels = MEGDataLabels_h.MEGLabels;

par.NComponents = 60;
par.fontSize = 8;

end


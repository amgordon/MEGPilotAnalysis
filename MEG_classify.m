function [ results ] = MEG_classify(par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


ft_defaults

res.hdr = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

cfg = [];
cfg.layout = par.dataFiles{1};


%note: put layout stuff in par structure...
for f = 1:length(par.dataRunsPreproc)
    disp(sprintf('loading run %g', f))
    data_h = load(par.dataRunsPreproc{f});

    T{f} = data_h.res.data_epochs.time; 
    R{f} = data_h.res.data_epochs.trial;    
end

allTrials = horzcat(R{:});
allTime = horzcat(T{:});

res = data_h.res;
res.data_epochs.trial = allTrials;
res.data_epochs.time = allTime;
res.data_epochs.trialinfo = (1:length(res.data_epochs.trial))';
res.layout = ft_prepare_layout(cfg);

%MEGDat = load(par.dataConcatNoArt);
behavDat = load(par.behavFile);

nTrialsPerTask = par.miniBlockLength * par.nCycles;
idx.study1 = [true(1, par.miniBlockLength * par.nCycles), false(1, 2*par.miniBlockLength * par.nCycles)];
idx.study2 = [false(1, par.miniBlockLength * par.nCycles), true(1, par.miniBlockLength * par.nCycles), false(1, par.miniBlockLength * par.nCycles)];
idx.test = [false(1, 2*par.miniBlockLength * par.nCycles), true(1, par.miniBlockLength * par.nCycles)];


for i=1:length(behavDat.theData.Test.conID);
    
    switch behavDat.theData.Study1.conID{i}(1)
        case 'F'
            study1Cond(i) = 1;
        case 'S'
            study1Cond(i) = 2;
        case 'O'
            study1Cond(i) = 3;
    end
    
    
    switch behavDat.theData.Study2.conID{i}(1)
        case 'F'
            study2Cond(i) = 1;
        case 'S'
            study2Cond(i) = 2;
        case 'O'
            study2Cond(i) = 3;
    end
    
    switch behavDat.theData.Test.conID{i}(1)
        case 'F'
            testCond(i) = 1;
        case 'S'
            testCond(i) = 2;
        case 'O'
            testCond(i) = 3;
    end
    
end
idx.cond = [study1Cond(1:nTrialsPerTask) study2Cond(1:nTrialsPerTask) testCond(1:nTrialsPerTask)];

idx.goodTrials = true(size(idx.cond));
idx.goodTrials(par.badTrials) = false;
% idx.goodTrial(MEGDat.preprocData.trialinfo) = true;

%%
par.timePointsForClassification = 4:13;
par.classificationFs = 10;

cfg = [];
cfg.detrend = 'no';
cfg.resamplefs = par.classificationFs;
res.continuous_data_bpf_rs = ft_resampledata(cfg, res.data_epochs);

%% Look at the ERPs

cfg = [];

idx.viewedFaces = (idx.cond==1) & (~idx.test) & idx.goodTrials;
idx.viewedScenes = (idx.cond==2) & (~idx.test) & idx.goodTrials;

viewedFaces = res.continuous_data_bpf_rs;
viewedScenes = res.continuous_data_bpf_rs;

viewedFaces.trial = viewedFaces.trial(idx.viewedFaces);
viewedScenes.trial = viewedScenes.trial(idx.viewedScenes);

avgRespFaces = ft_timelockanalysis(cfg, viewedFaces);
avgRespScenes = ft_timelockanalysis(cfg, viewedScenes);

cfg = [];
cfg.showlabels = 'yes'; 
cfg.fontsize = 6; 
cfg.layout = res.layout;
cfg.ylim = [-1e-14 1e-14];
ft_multiplotER(cfg, avgRespFaces, avgRespScenes); 

%% Classify!!

S.nXvals = 10;
S.balanceTrainingSet = false;
S.FS = false;
S.nFeats = 1000;
S.trainOptsLibLinear = '-s 0 -B 1 -q -c ';  
S.ValidationLambda = .01;

%% reshape the data matrix
X_h1 = cat(3,res.continuous_data_bpf_rs.trial{:});

X_h2 = X_h1(:,par.timePointsForClassification,:);

X = reshape(X_h2, [size(X_h2,1)*size(X_h2,2), size(X_h2,3)]);
sparseX = 10^14 * sparse(X)';


%% classify study vs. test
idx.study1AndTest = (~idx.study2) & idx.goodTrials;
cond = idx.study1 - idx.test;

results.classifyStudyVsTest = MEG_x_validation(sparseX(idx.study1AndTest,:),cond(idx.study1AndTest),S,S.ValidationLambda,'discrete', 'liblinear');

%% classify one timepoint vs. another
t1 = squeeze(X_h2(:,10,:));
t2 = squeeze(X_h2(:,1,:));

idx.notNan = [idx.goodTrials&idx.study1 idx.goodTrials&idx.study1];

sparseT = 10^14 * sparse(horzcat(t1,t2))';
idxT1vsT2 = [ones(1,900), -1*ones(1,900)];

results.oneTimepointVsAnother = MEG_x_validation(sparseT(idx.notNan,:),idxT1vsT2(idx.notNan),S,S.ValidationLambda,'discrete', 'liblinear');

%% classify face vs. house viewing
idx.faceVsHouseSt1 = (idx.cond<3) & idx.study1 & idx.goodTrials;
cond = (2*idx.cond)-3;

results.classifyFaceVsHouseStudy1 = MEG_x_validation(sparseX(idx.faceVsHouseSt1,:),cond(idx.faceVsHouseSt1),S,S.ValidationLambda,'discrete', 'liblinear');


%% classify face vs. house remembering
idx.faceVsHouseSt1 = (idx.cond<3) & idx.test & idx.goodTrials;
cond = (2*idx.cond)-3;

results.classifyFaceVsHouseTest = MEG_x_validation(sparseX(idx.faceVsHouseSt1,:),cond(idx.faceVsHouseSt1),S,S.ValidationLambda,'discrete', 'liblinear');


%% train on face vs. house viewing, test on face vs. house remembering
idx.faceVsHouseTrain = (idx.cond<3) & idx.goodTrials & ~idx.test;
idx.faceVsHouseTest = (idx.cond<3) & idx.goodTrials & idx.test;
cond = (2*idx.cond)-3;

model = train(cond(idx.faceVsHouseTrain)', sparseX(idx.faceVsHouseTrain,:), '-s 0 -B 1 -q -c 1');
[YP, X0E, results.trainOnStudyTestOnTest] = predict(cond(idx.faceVsHouseTest)', sparseX(idx.faceVsHouseTest,:), model);
end



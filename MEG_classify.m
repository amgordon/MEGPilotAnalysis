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

idxB = MEG_behav_analysis(par);

idx.cond = [idxB.study1.cond idxB.study2.cond idxB.test.cond];

idx.goodTrials = true(size(idx.cond));
idx.goodTrials(par.badTrials) = false;

idx.study = [true(size(idxB.study1.cond)) true(size(idxB.study2.cond)) false(size(idxB.test.cond))];
idx.test = ~idx.study;

idx.remember = false(size(idx.test));
idx.remember(idx.test) = idxB.test.remember;
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
S.FS = true;
S.nFeats = 1000;
S.trainOptsLibLinear = '-s 0 -B 1 -q -c ';  
S.ValidationLambda = [.1];

%% reshape the data matrix
X_h1 = cat(3,res.continuous_data_bpf_rs.trial{:});

X_h2 = X_h1(:,par.timePointsForClassification,:);

X = 10^14 * reshape(X_h2, [size(X_h2,1)*size(X_h2,2), size(X_h2,3)]);
sparseX = sparse(X)';


% %% classify study vs. test
cond = idx.study - idx.test; 
results.classifyStudyVsTest = MEG_x_validation(sparseX(idx.goodTrials,:),cond(idx.goodTrials),S,S.ValidationLambda,'discrete', 'liblinear');
 

%% classify face vs. house viewing
idx.faceVsHouseSt1 =  idx.study & idx.goodTrials;
cond = idx.cond;

results.classifyFaceVsHouseStudy1 = MEG_x_validation(sparseX(idx.faceVsHouseSt1,:),cond(idx.faceVsHouseSt1),S,S.ValidationLambda,'discrete', 'liblinear');


%% classify face vs. house remembering
idx.faceVsHouseSt1 =  idx.remember & idx.goodTrials;
cond = idx.cond;

results.classifyFaceVsHouseTest = MEG_x_validation(sparseX(idx.faceVsHouseSt1,:),cond(idx.faceVsHouseSt1),S,S.ValidationLambda,'discrete', 'liblinear');


%% train on face vs. house viewing, test on face vs. house remembering
idx.faceVsHouseTrain =  idx.goodTrials & ~idx.test;
idx.faceVsHouseTest =  idx.goodTrials & idx.test;
cond = idx.cond;

model = train(cond(idx.faceVsHouseTrain)', sparseX(idx.faceVsHouseTrain,:), '-s 0 -B 1 -q -c .0001');
[results.trainOnStudyTestOnTest] = predict(cond(idx.faceVsHouseTest)', sparseX(idx.faceVsHouseTest,:), model);

%% class-specific RSA analyses

    idx.faceTrain =  idx.goodTrials & ~idx.test & (idx.cond==1);
    idx.faceTest =  idx.goodTrials & idx.test & (idx.cond==1);
    idx.houseTrain =  idx.goodTrials & ~idx.test & (idx.cond==2);
    idx.houseTest =  idx.goodTrials & idx.test & (idx.cond==2);
    
for i=1:size(X_h2,2)
    [rFF(i,:,:)] = corr(squeeze(X_h2(:,i,idx.faceTrain)), squeeze(X_h2(:,i,idx.faceTest)));
    [rHH(i,:,:)] = corr(squeeze(X_h2(:,i,idx.houseTest)), squeeze(X_h2(:,i,idx.houseTest)));
    [rFH(i,:,:)] = corr(squeeze(X_h2(:,i,idx.faceTrain)), squeeze(X_h2(:,i,idx.houseTest)));
    [rHF(i,:,:)] = corr(squeeze(X_h2(:,i,idx.houseTest)), squeeze(X_h2(:,i,idx.faceTest)));
end

rFFM = squeeze(max(rFF));
rHHM = squeeze(max(rHH));
rFHM = squeeze(max(rFH));
rHFM = squeeze(max(rHF));

%% stim-specific RSA analysis




end


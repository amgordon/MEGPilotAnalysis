function [ output_args ] = MEG_classify( par )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

MEGDat = load(par.dataConcatNoArt);
behavDat = load(par.behavFile);

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
    
    idx.cond = [study1Cond study2Cond testCond];
end

idx.goodTrial = false(size(idx.cond));
idx.goodTrial(MEGDat.preprocData.trialinfo) = true;

%%
par.timePointsForClassification = 6:26;
par.classificationFs = 20;

cfg = [];
cfg.detrend = 'no';
cfg.resamplefs = par.classificationFs;
res.continuous_data_bpf_rs = ft_resampledata(cfg, MEGDat.preprocData);

%% Classify!!

S.nXvals = 10;
S.balanceTrainingSet = false;
S.FS = false;
S.nFeats = 1000;
S.trainOptsLibLinear = '-s 0 -B 1 -q -c ';  
S.ValidationLambda = .01;

%% reshape the data matrix
X_h0 = cat(3,res.continuous_data_bpf_rs.trial{:});
sz = size(X_h0);
X_h1 = nan(sz(1), sz(2), length(idx.cond));
X_h1(:,:,idx.goodTrial) = X_h0;

X_h2 = X_h1(:,par.timePointsForClassification,:);

X = reshape(X_h2, [size(X_h2,1)*size(X_h2,2), size(X_h2,3)]);
sparseX = 10^14 * sparse(X)';


%% classify one timepoint vs. another
t1 = squeeze(X_h2(:,1,:));
t2 = squeeze(X_h2(:,5,:));

idx.notNan = [idx.goodTrial idx.goodTrial];

sparseT = 10^14 * sparse(horzcat(t1,t2))';
idxT1vsT2 = [ones(1,900), -1*ones(1,900)];

results.oneTimepointVsAnother = MEG_x_validation(sparseT(idx.notNan,:),idxT1vsT2(idx.notNan),S,S.ValidationLambda,'discrete', 'liblinear');
%% classify face vs. house viewing
idx.faceVsHouseSt1 = (idx.cond<3) & idx.study1 & idx.goodTrial;
cond = (2*idx.cond)-3;

results.classifyFaceVsHouseStudy1 = MEG_x_validation(sparseX(idx.faceVsHouseSt1,:),cond(idx.faceVsHouseSt1),S,S.ValidationLambda,'discrete', 'liblinear');

end



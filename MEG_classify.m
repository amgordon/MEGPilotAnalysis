function [results, idx] = MEG_classify(par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


ft_defaults

res.hdr = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

cfg = [];
cfg.layout = par.dataFiles{1};
res.layout = ft_prepare_layout(cfg);

%note: put layout stuff in par structure...
for f = 1:length(par.dataRunsPreprocStudy)
    disp(sprintf('loading study run %g', f))
    data_hs = load(par.dataRunsPreprocStudy{f});

    Ts{f} = data_hs.res.data_epochs.time; 
    Rs{f} = data_hs.res.data_epochs.trial;    
end

for g = 1:length(par.dataRunsPreprocTest)
    disp(sprintf('loading test run %g', g))
    data_ht = load(par.dataRunsPreprocTest{g});

    Tt{g} = data_ht.res.data_epochs.time; 
    Rt{g} = data_ht.res.data_epochs.trial;    
end

%%
allTrials_study = horzcat(Rs{:});
allTime_study = horzcat(Ts{:});

resS = data_hs.res;
resS.data_epochs.trial = allTrials_study;
resS.data_epochs.time = allTime_study;
resS.data_epochs.trialinfo = (1:length(resS.data_epochs.trial))';
resS.layout = ft_prepare_layout(cfg);

%%
allTrials_test = horzcat(Rt{:});
allTime_test = horzcat(Tt{:});

resT = data_ht.res;
resT.data_epochs.trial = allTrials_test;
resT.data_epochs.time = allTime_test;
resT.data_epochs.trialinfo = (1:length(resT.data_epochs.trial))';
resT.layout = ft_prepare_layout(cfg);


%%
par.timePointsForClassification = 4:13;
par.classificationFs = 10;

cfg = [];
cfg.detrend = 'no';
cfg.resamplefs = par.classificationFs;

res.continuous_data_bpf_rs_study = ft_resampledata(cfg, resS.data_epochs);
res.continuous_data_bpf_rs_test = ft_resampledata(cfg, resT.data_epochs);

%% Classify!!

S.nXvals = 10;
S.balanceTrainingSet = true;
S.FS = false;
S.nFeats = 1000;
S.trainOptsLibLinear = '-s 0 -B 1 -q -c ';  
S.trainOptsSVM = '-s 3 -q -c ';  
S.ValidationLambda = [1];

%% reshape the data matrix
% fix trials that are missing data points because the RT was so late.
NBinsStudy_h = cellfun(@(x) size(x,2), res.continuous_data_bpf_rs_study.trial, 'UniformOutput', false);
NBinsStudy = max([NBinsStudy_h{:}]);
for i=1:length(res.continuous_data_bpf_rs_study.trial)
    if (size(res.continuous_data_bpf_rs_study.trial{i},2)<NBinsStudy)
        res.continuous_data_bpf_rs_study.trial{i} = nan(size(res.continuous_data_bpf_rs_study.trial{i},1), NBinsStudy);
        idx.insufficientTimePointsStudy(i) = true;
    else
        idx.insufficientTimePointsStudy(i) = false;
    end
end
Xs_h1 = cat(3,res.continuous_data_bpf_rs_study.trial{:});
Xs_h2 = Xs_h1(:,par.timePointsForClassification,:);

Xs = 10^16 * reshape(Xs_h2, [size(Xs_h2,1)*size(Xs_h2,2), size(Xs_h2,3)]);
sparseXs = sparse(Xs)';


% fix trials that are missing data points because the RT was so late.
NBinsTest_h = cellfun(@(x) size(x,2), res.continuous_data_bpf_rs_test.trial, 'UniformOutput', false);
NBinsTest = max([NBinsTest_h{:}]);
for i=1:length(res.continuous_data_bpf_rs_test.trial)
    if (size(res.continuous_data_bpf_rs_test.trial{i},2)<NBinsTest)
        res.continuous_data_bpf_rs_test.trial{i} = nan(size(res.continuous_data_bpf_rs_test.trial{i},1), NBinsTest);
        idx.insufficientTimePointsTest(i) = true;
    else
        idx.insufficientTimePointsTest(i) = false;
    end
end
Xt_h1 = cat(3,res.continuous_data_bpf_rs_test.trial{:});
Xt_h2 = Xt_h1(:,par.timePointsForClassification,:);

Xt = 10^16 * reshape(Xt_h2, [size(Xt_h2,1)*size(Xt_h2,2), size(Xt_h2,3)]);
sparseXt = sparse(Xt)';


% %% classify study vs. test
% cond = idx.study - idx.test; 
% results.classifyStudyVsTest = MEG_x_validation(sparseX(idx.goodTrials,:),cond(idx.goodTrials),S,S.ValidationLambda,'discrete', 'liblinear');
%  
%%
idxB = MEG_behav_analysis(par);
NTrials = length(idxB.trialNum);

idx.cond = [idxB.study1.cond idxB.study2.cond idxB.test.cond];
idx.goodTrials = true(size(idx.cond));
idx.goodTrials(par.badTrials) = false;
idx.study = [true(size(idxB.study1.cond)) true(size(idxB.study2.cond)) false(size(idxB.test.cond))];
idx.study1 = [true(size(idxB.study1.cond)) false(size(idxB.study2.cond)) false(size(idxB.test.cond))];
idx.study2 = [false(size(idxB.study1.cond)) true(size(idxB.study2.cond)) false(size(idxB.test.cond))];
idx.test = ~idx.study;
idx.remember = false(size(idx.test));
idx.remember(idx.test) = idxB.test.remember;
idx.testRemember = idxB.test.remember;
idx.SM = [idxB.study1.SM idxB.study2.SM];

idx.condStudy = [idxB.study1.cond idxB.study2.cond];
idx.condTest = idxB.test.cond;
idx.goodTrialsStudy = idx.goodTrials(idx.study) & ~idx.insufficientTimePointsStudy;
idx.goodTrialsTest = idx.goodTrials(idx.test) &  ~idx.insufficientTimePointsTest;

%% Look at the Study ERPs

cfg = [];
idx.viewedFaces = (idx.condStudy==1) & idx.goodTrialsStudy;
idx.viewedScenes = (idx.condStudy==2) & idx.goodTrialsStudy;

viewedFaces = res.continuous_data_bpf_rs_study;
viewedScenes = res.continuous_data_bpf_rs_study;

viewedFaces.trial = viewedFaces.trial(idx.viewedFaces);
viewedFaces.time = viewedFaces.time(idx.viewedFaces);
viewedScenes.trial = viewedScenes.trial(idx.viewedScenes);
viewedScenes.time = viewedScenes.time(idx.viewedScenes);

avgRespFaces = ft_timelockanalysis(cfg, viewedFaces);
avgRespScenes = ft_timelockanalysis(cfg, viewedScenes);

cfg = [];
cfg.showlabels = 'yes'; 
cfg.fontsize = 6; 
cfg.layout = res.layout;
cfg.ylim = [-1e-14 1e-14];
figure; ft_multiplotER(cfg, avgRespFaces, avgRespScenes); 

%% Look at the Test ERPs
idx.RTBins = makeQuantileLabels(idxB.test.RT,3);

cfg = [];

idx.RT1 = (idx.RTBins==1) &  idx.goodTrialsTest;
idx.RT3 = (idx.RTBins==3) & idx.goodTrialsTest;
%idx.viewedFaces = (idx.cond==1) & (~idx.test) & idx.goodTrials;
%idx.viewedScenes = (idx.cond==2) & (~idx.test) & idx.goodTrials;

testTrialRT1 = res.continuous_data_bpf_rs_test;
testTrialRT3 = res.continuous_data_bpf_rs_test;
%viewedScenes = res.continuous_data_bpf_rs_test;

testTrialRT1.trial = testTrialRT1.trial(idx.RT1);
testTrialRT1.time = testTrialRT1.time(idx.RT1);
testTrialRT3.trial = testTrialRT3.trial(idx.RT3);
testTrialRT3.time = testTrialRT3.time(idx.RT3);
%viewedScenes.trial = viewedScenes.trial(idx.viewedScenes);

avgRespTestRT1 = ft_timelockanalysis(cfg, testTrialRT1);
avgRespTestRT3 = ft_timelockanalysis(cfg, testTrialRT3);
%avgRespScenes = ft_timelockanalysis(cfg, viewedScenes);

cfg = [];
cfg.showlabels = 'yes'; 
cfg.fontsize = 6; 
cfg.layout = res.layout;
cfg.ylim = [-1e-14 1e-14];
figure; ft_multiplotER(cfg, avgRespTestRT1, avgRespTestRT3); 

%% classify face vs. house viewing
% idx.faceVsHouseSt1 =  idx.study & idx.goodTrials;
% cond = idx.cond;
% 
% 
% 
%% classify face vs. house viewing
% conds = idx.condStudy';
%  
% thisX_h = reshape(Xs_h1(:,:,idx.goodTrialsStudy), size(Xs_h1,1)* size(Xs_h1,2), sum(idx.goodTrialsStudy));
% thisX = 10^16 * sparse(thisX_h)';
% results.classifyStudyStims =  MEG_x_validation(...
%     thisX,conds(idx.goodTrialsStudy),S,S.ValidationLambda,'discrete', 'liblinear');
% 
% 
% yProb = results.classifyStudyStims.mod.YProb;
% trueY = results.classifyStudyStims.mod.Y;
% labSet = results.classifyStudyStims.mod.model.Label;
% [~, labSetInv] = sort(labSet);
% for i=1:length(yProb)
%    guessCorClass(i) = yProb(i,labSetInv(trueY(i)));    
% end
% 
% [~,~,~,results.encStrengthAndSM] = ttest2(guessCorClass(~idx.SM(idx.goodTrialsStudy)), guessCorClass(idx.SM(idx.goodTrialsStudy)));
% 

%% classify Subs Correct vs. Incorrect
% conds = double(idx.SM');
% thisX_h = reshape(Xs_h1(:,:,idx.goodTrialsStudy), size(Xs_h1,1)* size(Xs_h1,2), sum(idx.goodTrialsStudy));
% thisX = 10^16 * sparse(thisX_h)';
% results.classifySubsMem =  MEG_x_validation(...
%     thisX,conds(idx.goodTrialsStudy),S,S.ValidationLambda,'discrete', 'liblinear');

%% train on face vs. house viewing, test on face vs. house remembering
% idx.faceVsHouseTrain =  idx.goodTrials & idx.study1;
% idx.faceVsHouseTest =  idx.goodTrials & idx.study2;
% cond = idx.cond;
% 
% model = train(cond(idx.faceVsHouseTrain)', sparseX(idx.faceVsHouseTrain,:), '-s 0 -B 1 -q -c .0001');
% [results.trainOnStudy1TestOnStudy2] = predict(cond(idx.faceVsHouseTest)', sparseX(idx.faceVsHouseTest,:), model);

%% train and test on each study timepoint
% conds = idx.condStudy';
% 
% for j = 1:size(Xs_h1,2)
%     thisX = zscore(sparse(squeeze(Xs_h1(:,j, idx.goodTrialsStudy))))';
%     results.classifyFaceVsHouseStudyTimebins(j) =  MEG_x_validation(...
%         thisX,conds( idx.goodTrialsStudy),S,S.ValidationLambda,'discrete', 'liblinear');
% end

%% train on face vs. house viewing, test on face vs. house remembering
% labels = res.layout.label(1:157);
% idx.ChansToInclude = labels;

% idx.faceVsHouseTrain =  idx.goodTrialsStudy;
% idx.faceVsHouseTest =  idx.goodTrialsTest;
% condsTrain = idx.condStudy;
% condsTest = idx.condTest;
% 
% thisXTrain = 10^16 * sparse(squeeze(mean(Xs_h1(:,5:6,idx.faceVsHouseTrain),2)))';
% model = train(condsTrain(idx.faceVsHouseTrain)', thisXTrain, '-s 0 -B 1 -q -c .01');
% for j = 1:size(Xt_h1,2)
%     thisXTest = 10^16 * sparse(squeeze(Xt_h1(:,j,idx.faceVsHouseTest)))';
%     [a1{j}, results.trainOnStudy1TestOnTest2(j), a3{j}] = predict(condsTest(idx.faceVsHouseTest)', thisXTest, model);
% end
% 
%% train and test on each test timepoint
%   conds = idx.condTest';
%   for j = 1:size(Xt_h1,2)
%       thisX = 10^16 * sparse(squeeze(Xt_h1(:,j,idx.goodTrialsTest)))';
%       results.classifyFaceVsHouseTestTimebins(j) =  MEG_x_validation(...
%           thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
%   end

%% classify timepoints



for i=1:20
    for j=1:20
        thisX_t1 = 10^16 * sparse(squeeze(Xt_h1(:,i,idx.goodTrialsTest)))';
        thisX_t2 = 10^16 * sparse(squeeze(Xt_h1(:,j,idx.goodTrialsTest)))';
        thisX = cat(1,thisX_t1,thisX_t2);
        conds = [ones(1,size(thisX_t1,1)), zeros(1,size(thisX_t1,1))];
        conds = conds([idx.goodTrialsTest idx.goodTrialsTest])';
        
        results.classifyTwoTimepoints(i,j) =  MEG_x_validation(...
            thisX,conds,S,S.ValidationLambda,'discrete', 'liblinear');
    end
end


% %% remembered vs. forgotten by timpoint
% conds = double(2-idxB.test.remember)';
% for j = 1:size(Xt_h1,2)
%     thisX = 10^16 * sparse(squeeze(Xt_h1(:,j,idx.goodTrialsTest)))';
%     results.classifyRememberVsForgottenTimepoints(j) =  MEG_x_validation(...
% thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
% end
% 
% %% slow vs. fast RT discrete by timepoint
% conds = makeQuantileLabels(idxB.test.RT,3)';
% idx.goodTrialsTest = idx.goodTrialsTest & (conds'~=2);
% for j = 1:size(Xt_h1,2)
%     thisX = 10^16 * sparse(squeeze(Xt_h1(:,j,idx.goodTrialsTest)))';
%     results.classifyRTDiscreteByTimepoint(j) =  MEG_x_validation(...
%         thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
% end
% 
% %% slow vs. fast RT discrete by channel
% conds = makeQuantileLabels(idxB.test.RT,3)';
% idx.goodTrialsTest = idx.goodTrialsTest & (conds'~=2);
% for j = 1:size(Xt_h1,1)
%     thisX = 10^16 * sparse(squeeze(Xt_h1(j,:,idx.goodTrialsTest)))';
%     results.classifyRTDiscreteByChannel(j) =  MEG_x_validation(...
%         thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
% end
% 
% %% logistic regression to predict RT
% conds = makeQuantileLabels(idxB.test.RT,3)';
% idx.goodTrialsTest = idx.goodTrialsTest & (conds'~=2);
% % 
%  results.classifyTestRTDiscrete =  MEG_x_validation(...
%      sparseXt(idx.goodTrialsTest,:),conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
% 
%  
% %% linear regression to predict RT
% conds = idxB.test.RT';
% % 
%  results.classifyTestRT =  MEG_x_validation(...
%      zscore(sparseXt(idx.goodTrialsTest,:)),conds(idx.goodTrialsTest),S,10,'continuous', 'svr');


% 
% %% class-specific RSA analyses
% 
% idx.faceTrain =  idx.goodTrials & ~idx.test & (idx.cond==1);
% idx.faceTest =  idx.goodTrials & idx.test & (idx.cond==1);
% idx.houseTrain =  idx.goodTrials & ~idx.test & (idx.cond==2);
% idx.houseTest =  idx.goodTrials & idx.test & (idx.cond==2);
% 
% for i=1:size(X_h2,2)
%     [rFF(i,:,:)] = corr(squeeze(X_h2(:,i,idx.faceTrain)), squeeze(X_h2(:,i,idx.faceTest)));
%     [rHH(i,:,:)] = corr(squeeze(X_h2(:,i,idx.houseTest)), squeeze(X_h2(:,i,idx.houseTest)));
%     [rFH(i,:,:)] = corr(squeeze(X_h2(:,i,idx.faceTrain)), squeeze(X_h2(:,i,idx.houseTest)));
%     [rHF(i,:,:)] = corr(squeeze(X_h2(:,i,idx.houseTest)), squeeze(X_h2(:,i,idx.faceTest)));
% end
% 
% rFFM = squeeze(max(rFF));
% rHHM = squeeze(max(rHH));
% rFHM = squeeze(max(rFH));
% rHFM = squeeze(max(rHF));

%% stim-specific RSA analysis
% for i=1:length(idxB.test.word)
%     thisWord = idxB.test.word{i};
%     idxThisWordTest = i + NTrials*2;
%     
%     idxThisWordStudy1 = find(strcmp(thisWord, idxB.study1.word));
%     idxThisWordStudy2 = find(strcmp(thisWord, idxB.study2.word)) + NTrials;
%     
%     idxOtherWordsStudy1 = find(~strcmp(thisWord, idxB.study1.word));
%     idxOtherWordsStudy2 = find(~strcmp(thisWord, idxB.study2.word)) + NTrials;
%     
%     for j=1:size(X_h2,2)
%         for k=1:size(X_h2,2)
%             r_S1_T_sameWord(i,j,k) = corr(squeeze(X_h2(:,j,idxThisWordTest)), squeeze(X_h2(:,k,idxThisWordStudy1)));
%             r_S2_T_sameWord(i,j,k) = corr(squeeze(X_h2(:,j,idxThisWordTest)), squeeze(X_h2(:,k,idxThisWordStudy2)));
%             r_S1_S2_sameWord(i,j,k) = corr(squeeze(X_h2(:,j,idxThisWordStudy1)), squeeze(X_h2(:,k,idxThisWordStudy2)));
%             
%             r_S1_T_diffWord(i,j,k) = mean(corr(squeeze(X_h2(:,j,idxThisWordTest)), squeeze(X_h2(:,k,idxOtherWordsStudy1))));
%             r_S2_T_diffWord(i,j,k) = mean(corr(squeeze(X_h2(:,j,idxThisWordTest)), squeeze(X_h2(:,k,idxOtherWordsStudy2))));
%             r_S1_S2_diffWord(i,j,k) = mean(corr(squeeze(X_h2(:,j,idxThisWordStudy1)), squeeze(X_h2(:,k,idxOtherWordsStudy2))));
%             
%             word(i,j,k) = i;
%             timeBinT(i,j,k) = j;
%             timeBinS(i,j,k) = k;
%         end
%     end
% end
% 
% r_S1_T_sameWordAvg = squeeze(mean(r_S1_T_sameWord,1));
% r_S2_T_sameWordAvg = squeeze(mean(r_S2_T_sameWord,1));
% r_S1_S2_sameWordAvg = squeeze(mean(r_S1_S2_sameWord,1));
% 
% S_T_sameWord = .5*r_S1_T_sameWord + .5*r_S2_T_sameWord;
% 
% r_S1_T_diffWordAvg = squeeze(mean(r_S1_T_diffWord,1));
% r_S2_T_diffWordAvg = squeeze(mean(r_S2_T_diffWord,1));
% r_S1_S2_diffWordAvg = squeeze(mean(r_S1_S2_diffWord,1));
% 
% S_T_diffWord = .5*r_S1_T_diffWord + .5*r_S2_T_diffWord;
% 
% for i=1:length(idxB.test.word)
%     sameWordDiag(i,:) = diag(squeeze(S_T_sameWord(i,:,:)));
%     diffWordDiag(i,:) = diag(squeeze(S_T_diffWord(i,:,:)));
% end
% 
% for j=1:size(X_h2,2)
%     for k=1:size(X_h2,2)
%         [a1{j,k}, a2{j,k} a3{j,k} a4(j,k)] = ttest(squeeze(r_S1_T_sameWord(:,j,k)), squeeze(r_S1_T_diffWord(:,j,k)));
%         b(j,k) = a4(j,k).tstat;
%     end
% end
% 
% TOI = {2:3, 3:4};
% 
% results.sameWord_studyTest = squeeze(mean(mean(S_T_sameWord(:,TOI{1}, TOI{2}),3),2));
% results.diffWord_studyTest = squeeze(mean(mean(S_T_diffWord(:,TOI{1}, TOI{2}),3),2));
% 

% toPrint(:,1) = ['S1_T_sameWord' num2cell(r_S1_T_sameWord(1:end))];
% toPrint(:,2) = ['S1_T_diffWord' num2cell(r_S1_T_diffWord(1:end))];
% toPrint(:,3) = ['S2_T_sameWord' num2cell(r_S2_T_sameWord(1:end))];
% toPrint(:,4) = ['S2_T_diffWord' num2cell(r_S2_T_diffWord(1:end))];
% toPrint(:,5) = ['word' num2cell(word(1:end))];
% toPrint(:,6) = ['timeBinT' num2cell(timeBinT(1:end))];
% toPrint(:,7) = ['timeBinS' num2cell(timeBinS(1:end))];
% 
% cell2csv('/biac4/wagner/biac3/wagner5/alan/MEG/summaryDataForR/dataS1.csv', toPrint, ',', 2000);
save (['/biac4/wagner/biac3/wagner5/alan/MEG/PilotData/classifications/sub' num2str(par.subNo) 'stimLockAmp'], 'results')
disp('done');

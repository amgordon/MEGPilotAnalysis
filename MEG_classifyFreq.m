function [results, idx] = MEG_classifyFreq(par)
% classify MEG data based on spectral properties.  

ft_defaults

res.hdr = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

cfg = [];
cfg.layout = par.dataFiles{1};
res.layout = ft_prepare_layout(cfg);

%note: put layout stuff in par structure...
for f = 1:length(par.dataRunsPreprocFreqStudy)
    disp(sprintf('loading study run %g', f));
    data_hs = load(par.dataRunsPreprocFreqStudy{f});
    
    Ts{f} = data_hs.res.continuous_data_freq.powspctrm;
end

allDat_study = 10^26*cat(1, Ts{:});
clear Ts;

for g = 1:length(par.dataRunsPreprocFreqTest)
    disp(sprintf('loading test run %g', g));
    data_ht = load(par.dataRunsPreprocFreqTest{g});
    
    Tt{g} = data_ht.res.continuous_data_freq.powspctrm;
end


allDat_test = 10^26*cat(1,Tt{:});
clear Tt;

%%
par.timePointsForClassification = 4:13;

%% Classify!!
S.nXvals = 10;
S.balanceTrainingSet = false;
S.FS = false;
S.nFeats = 1000;
S.trainOptsLibLinear = '-s 0 -B 1 -q -c ';  
S.trainOptsSVM = '-s 3 -q -c ';  
S.ValidationLambda = [1];

%% reshape the data matrix
Xs_h1 = shiftdim(allDat_study,3);
Xt_h1 = shiftdim(allDat_test,3);

szS = size(Xs_h1);
szT = size(Xt_h1);

Xs_h2 = nan(size(Xs_h1));
Xt_h2 = nan(size(Xt_h1));
for i=1:szS(3)
    for j=1:szS(4)
        thisChFreqXs = squeeze(Xs_h1(:,:,i,j));
        zOverTimeXs = zscore(thisChFreqXs(1:end));
        Xs_h2(:,:,i,j) = reshape(zOverTimeXs,szS(1), szS(2));
        
        thisChFreqXt = squeeze(Xt_h1(:,:,i,j));
        zOverTimeXt = zscore(thisChFreqXt(1:end));
        Xt_h2(:,:,i,j) = reshape(zOverTimeXt,szT(1), szT(2));
    end
end

Xs = reshape(Xs_h1, [szS(1), szS(2), szS(3)*szS(4)]);
Xt = reshape(Xt_h1, [szT(1), szT(2), szT(3)*szT(4)]);


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
idx.SM = [idxB.study1.SM idxB.study2.SM];

idx.condStudy = [idxB.study1.cond idxB.study2.cond];
idx.condTest = idxB.test.cond;
idx.goodTrialsStudy = idx.goodTrials(idx.study);
idx.goodTrialsTest = idx.goodTrials(idx.test);
idx.rememberTest = idx.remember(idx.test);

a = find(isnan(allDat_test));
[i1,~,~,~] = ind2sub(size(allDat_test),a);

idx.goodTrialsTest(i1) = false;
% % %% study vs. test by channel and timepoint
% conds = double(idx.test');
% X = cat(1,allDat_study,allDat_test);
% 
% for j = 1:size(X,3)
%     for k = 1:size(X,4)
%         thisX = sparse(squeeze(X(idx.goodTrials,:,j,k)));
%         results.classifyFaceVsHouseStudyTimebins(j,k) =  MEG_x_validation(...
%             thisX,conds(idx.goodTrials),S,S.ValidationLambda,'discrete', 'liblinear');
%     end
% end


% %% viewed class by channel and timepoint
% conds = idx.condStudy';
% X = Xs_h2;
% 
% for j = 1:size(X,1)
%     for k = 1:size(X,4)
%         thisX = sparse(squeeze(X(j,idx.goodTrialsStudy,:,k)));
%         results.classifyViewedClassEachChanAndBin(j,k) =  MEG_x_validation(...
%             thisX,conds(idx.goodTrialsStudy),S,S.ValidationLambda,'discrete', 'liblinear');
%     end
% end
% 
% %% viewed class by timepoint
% conds = idx.condStudy';
% X = Xs_h2;
% 
% idx.goodFreqs = 4:6;
% 
% for k = 1:size(X,1)
%     thisX_h1 = squeeze(X(k,idx.goodTrialsStudy,:,idx.goodFreqs));
%     szXh1 = size(thisX_h1);
%     thisX = sparse(reshape(thisX_h1, szXh1(1), szXh1(2)*szXh1(3)));
%     results.classifyViewedClass(k) =  MEG_x_validation(...
%         thisX,conds(idx.goodTrialsStudy),S,S.ValidationLambda,'discrete', 'liblinear');
% end
% 
% %% viewed class all timepoints
% conds = idx.condStudy';
% Xh = shiftdim(Xs_h2(:,idx.goodTrialsStudy,:,:),1);
% szXh = size(Xh);
% thisX = sparse(reshape(Xh, szXh(1), prod(szXh(2:4))));
% 
% results.classifyViewedClassAllDat =  MEG_x_validation(...
%     thisX,conds(idx.goodTrialsStudy),S,S.ValidationLambda,'discrete', 'liblinear');


% %% train on max study timepoint, test on test timepoints

idx.faceVsHouseTrain =  idx.goodTrialsStudy;
idx.faceVsHouseTest =  idx.goodTrialsTest;
condsTrain = idx.condStudy;
condsTest = idx.condTest;

if strcmp(par.substr, 'meg_081913')
    idx.goodTimes = 10:14;
else
    idx.goodTimes = 5;
end

thisXTrain_h = squeeze(mean(Xs_h1(idx.goodTimes,idx.goodTrialsStudy,:,:),1));
szXTrH = size(thisXTrain_h);
thisXTrain = sparse(reshape(thisXTrain_h, szXTrH(1), szXTrH(2)*szXTrH(3)));

model = train(condsTrain(idx.goodTrialsStudy)', thisXTrain, '-s 0 -B 1 -q -c .01');
for j = 1:size(Xt,1)
    thisXTest_h = squeeze(Xt_h1(j,:,:,:));
    szXTeH = size(thisXTest_h);
    thisXTest = sparse(reshape(thisXTest_h, szXTeH(1), szXTeH(2)*szXTeH(3)));
    [z1{j}, results.trainOnStudyTestOnTest(j), z3{j}] = predict(condsTest', thisXTest, model);
end

% 
% %% retrieved class by channel and timepoint
% conds = idx.condStudy';
% X = Xt_h1;
% 
% for j = 1:size(X,1)
%     for k = 1:size(X,4)
%         thisX = sparse(squeeze(X(j,idx.goodTrialsTest,:,k)));
%         results.classifyViewedClassEachChanAndBin(j,k) =  MEG_x_validation(...
%             thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
%     end
% end
% 
% 
% %% train and test on each test timepoint
% conds = idx.condTest';
% X = Xt_h1;
% 
% for k = 1:size(X,1)'
%     thisX_h1 = squeeze(X(k,idx.goodTrialsTest,:,:));
%     szXh1 = size(thisX_h1);
%     thisX = sparse(reshape(thisX_h1, szXh1(1), szXh1(2)*szXh1(3)));
%     results.classifyRetrievedClass(k) =  MEG_x_validation(...
%         thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
% end


%% remembered vs. forgotten by timpoint
% conds = double(idx.rememberTest)';
% for j = 1:size(Xt_h1,1)
%     thisX_h1 = squeeze(Xt_h1(j,idx.goodTrialsTest,:,:));
%     szXh1 = size(thisX_h1);
%     thisX = sparse(reshape(thisX_h1, szXh1(1), szXh1(2)*szXh1(3)));
%     results.classifyRememberVsForgottenTimepoints(j) =  MEG_x_validation(...
%         thisX,conds(idx.goodTrialsTest),S,S.ValidationLambda,'discrete', 'liblinear');
% end

%% linear regression to predict RT
% conds = idxB.test.RT';
% 
% results.classifyTestRT =  MEG_x_validation(...
%     sparseXt(idx.goodTrialsTest,:),conds(idx.goodTrialsTest),S,S.ValidationLambda,'continuous', 'svr');


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

save (['/biac4/wagner/biac3/wagner5/alan/MEG/PilotData/classifications/sub' num2str(par.subNo) 'RTLock'], 'results')

disp('done');

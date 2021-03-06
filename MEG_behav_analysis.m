function [idx] = MEG_behav_analysis(par)
% analyze behavioral data for a subject

% load raw behavioral data
behavDat = load(par.behavFile);

idx.trialNum = 1:length(behavDat.theData.Test.trialOnset);


study1Cond = conID2Cond(behavDat.theData.Study1.imageID);
study2Cond = conID2Cond(behavDat.theData.Study2.imageID);
testCond = conID2Cond(behavDat.theData.Test.imageID);

% image conditions
idx.study1.cond = study1Cond(idx.trialNum);
idx.study2.cond = study2Cond(idx.trialNum);
idx.test.cond = testCond(idx.trialNum);

% word
idx.study1.word = behavDat.theData.Study1.word(idx.trialNum);
idx.study2.word = behavDat.theData.Study2.word(idx.trialNum);
idx.test.word = behavDat.theData.Test.word(idx.trialNum);

% image label
idx.study1.image = behavDat.theData.Study1.image(idx.trialNum);
idx.study2.image = behavDat.theData.Study2.image(idx.trialNum);
idx.test.image = behavDat.theData.Test.image(idx.trialNum);

% mini block
idx.study1.miniblock = [behavDat.theData.Study1.miniblock{idx.trialNum}];
idx.study2.miniblock = [behavDat.theData.Study2.miniblock{idx.trialNum}];
idx.test.miniblock = [behavDat.theData.Test.miniblock{idx.trialNum}];

% onset
idx.study1.trialOnset = behavDat.theData.Study1.trialOnset(idx.trialNum);
idx.study2.trialOnset = behavDat.theData.Study2.trialOnset(idx.trialNum);
idx.test.trialOnset = behavDat.theData.Test.trialOnset(idx.trialNum);

% retrieval status
idx.test.remember = logical(behavDat.theData.Test.remember);
idx.test.forgot = logical(behavDat.theData.Test.forgot);

% retrieval RT
idx.test.RT = behavDat.theData.Test.RT;

% first or second study presentation?
for i=1:length(idx.study1.word)
    thisWordIdxTest = (strcmp(idx.study1.word(i), idx.test.word));
    idx.study1.SM(i) = idx.test.remember(thisWordIdxTest);    
end

for i=1:length(idx.study2.word)
    thisWordIdxTest = (strcmp(idx.study2.word(i), idx.test.word));
    idx.study2.SM(i) = idx.test.remember(thisWordIdxTest);    
end

end





function cond = conID2Cond(imageIDVec)
% given an image, convert the image to a condition label.  

cond = nan(size(imageIDVec));
for i=1:length(imageIDVec)
    switch (imageIDVec{i}(1))
        case 'f'
            cond(i) = 1;
        case 's'
            cond(i) = 2;
        case 'o'
            cond(i) = 3;
    end
end
end







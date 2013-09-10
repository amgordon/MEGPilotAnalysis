function [idx] = MEG_behav_analysis(par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


behavDat = load(par.behavFile);

idx.trialNum = 1:length(behavDat.theData.Test.trialOnset);

study1Cond = conID2Cond(behavDat.theData.Study1.imageID);
study2Cond = conID2Cond(behavDat.theData.Study2.imageID);
testCond = conID2Cond(behavDat.theData.Test.imageID);

idx.study1.cond = study1Cond(idx.trialNum);
idx.study2.cond = study2Cond(idx.trialNum);
idx.test.cond = testCond(idx.trialNum);

idx.study1.word = behavDat.theData.Study1.word(idx.trialNum);
idx.study2.word = behavDat.theData.Study2.word(idx.trialNum);
idx.test.word = behavDat.theData.Test.word(idx.trialNum);

idx.study1.image = behavDat.theData.Study1.image(idx.trialNum);
idx.study2.image = behavDat.theData.Study2.image(idx.trialNum);
idx.test.image = behavDat.theData.Test.image(idx.trialNum);

idx.study1.miniblock = behavDat.theData.Study1.miniblock(idx.trialNum);
idx.study2.miniblock = behavDat.theData.Study2.miniblock(idx.trialNum);
idx.test.miniblock = behavDat.theData.Test.miniblock(idx.trialNum);

idx.study1.trialOnset = behavDat.theData.Study1.trialOnset(idx.trialNum);
idx.study2.trialOnset = behavDat.theData.Study2.trialOnset(idx.trialNum);
idx.test.trialOnset = behavDat.theData.Test.trialOnset(idx.trialNum);

idx.test.remember = logical(behavDat.theData.Test.remember);
idx.test.forgot = logical(behavDat.theData.Test.forgot);
idx.test.RT = behavDat.theData.Test.RT;

end





function cond = conID2Cond(imageIDVec)

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







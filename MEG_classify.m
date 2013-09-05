function [ output_args ] = MEG_classify( par )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

MEGDat = load(par.dataConcatNoArt);
behavDat = load(par.behavFile);

idx.study1 = [ones(1, par.miniBlockLength * par.nCycles), zeros(1, 2*par.miniBlockLength * par.nCycles)];
idx.study2 = [zeros(1, par.miniBlockLength * par.nCycles), ones(1, par.miniBlockLength * par.nCycles), zeros(1, par.miniBlockLength * par.nCycles)];
idx.test = [zeros(1, 2*par.miniBlockLength * par.nCycles), ones(1, par.miniBlockLength * par.nCycles)];


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

end

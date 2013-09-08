function [ output_args ] =  MEG_wholeshebang(par, flags)

if ismember(flags, 's'), MEG_splitSqdIntoRuns(par), end;

if ismember(flags, 'p'), MEG_preprocessing_byRun2(par), end;
    
if ismember(flags, 'c'), MEG_classify(par), end;
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


end


function [ output_args ] =  MEG_wholeshebang(par, flags)

if ismember(flags, 's'), MEG_splitSqdIntoRuns(par), end;

if ismember(flags, 'p'), MEG_preprocessing_byRun(par), end;
    
if ismember(flags, 'c'), MEG_concatenatePreprocessedRuns(par), end;
 
if ismember(flags, 'a'), MEG_preprocess_acrossRun(par), end;
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


end


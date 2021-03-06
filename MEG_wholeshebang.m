function [ output_args ] =  MEG_wholeshebang(par, flags)
% perform multiple levels of analysis within the same subject.  

if ismember(flags, 's'), MEG_splitSqdIntoRuns(par), end;

if ismember(flags, 'p'), MEG_preprocessing_byRun2(par), end;
    
if ismember(flags, 'c'), MEG_classify(par), end;

if ismember(flags, 'w'), MEG_preprocessing_FreqDecomp(par), end;

end


function [res idx] = MEG_batchshebang(sa)

if nargin<1
    sa = 1:3;
end

for s=1:length(sa)
    par = MEG_par(sa(s));
    %MEG_preprocessing_byRun2(par);
%     par = MEG_par(sa(s));
    %MEG_preprocessing_FreqDecomp(par);
    %MEG_wholeshebang(par, 'spc');
    [res{s}, idx{s}] = MEG_classify(par);
end
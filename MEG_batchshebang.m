function [res idx] = MEG_batchshebang(sa)
% batch analysis across subjects.

if nargin<1
    sa = 1:3;
end

for s=1:length(sa)
    par = MEG_par(sa(s));
    [res{s}, idx{s}] = MEG_classify(par);
end
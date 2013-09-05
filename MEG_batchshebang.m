function [] = MEG_batchshebang(sa)

if nargin<1
    sa = 1:3;
end

for s=1:length(sa)
    par = MEG_par(sa(s));
    MEG_wholeshebang(par, 'p');
end
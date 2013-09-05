function [hdr] = MEG_make_hdr(par)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

hdr = ft_read_header(par.dataFiles{1});

if ~exist(par.headerDir, 'dir')
    mkdir(par.headerDir)
end

save(fullfile(par.headerDir, 'hdr.mat'), 'hdr')

end


function [hdr] = MEG_make_hdr(par)
% read and save a hdr for a data file.

hdr = ft_read_header(par.dataFiles{1});

if ~exist(par.headerDir, 'dir')
    mkdir(par.headerDir)
end

save(fullfile(par.headerDir, 'hdr.mat'), 'hdr')

end


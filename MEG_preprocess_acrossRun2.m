function [ output_args ] = MEG_preprocess_acrossRun2(par, flags)

%% MEG sample data preprocessing
%   Preprocess MEG using mainly fieldtrip routines
%   Alan Gordon, Stanford University, 09/02/2013
%   Adapted from scripts by Andy Heusser, NYU

if (nargin<2)
    flags = 'tami';
end

%concatData = load(par.dataConcat);

ft_defaults

res.hdr = ft_read_header(par.dataFiles{1});
res.grad = ft_read_sens(par.dataFiles{1},'filename',par.fileformat);

cfg = [];
cfg.layout = par.dataFiles{1};
res.layout = ft_prepare_layout(cfg);

%note: put layout stuff in par structure...
for f = 1:length(par.dataRunsPreproc)
    disp(sprintf('loading run %g', f))
    data_h = load(par.dataRunsPreproc{f});

    R{f} = data_h.res.data_epochs.trial;    
end

allTrials = horzcat(R{:});



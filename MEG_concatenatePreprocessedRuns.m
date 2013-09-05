function [] = MEG_concatenatePreprocessedRuns(par)

% Given a set of filtered miniblocks, concatenate the data and events back
% into a single timeseries.


nPreviousSamples = 0;
ix = 0;

for f = 1:length(par.dataRunsPreproc)
    
    disp(sprintf('loading run %g', f))
    
    data_h = ft_read_data(par.dataRunsPreproc{f});
    hdr_h = ft_read_header(par.dataRunsPreproc{f});
    events_h = ft_read_event(par.dataRunsPreproc{f});
    
    %R_h = load(par.dataRunsPreproc{f});
    
    origSampleRate = par.acquisitionRate;
    curSampleRate = hdr_h.Fs;
    
    %resampled rate / original rate of acquisition
    compressionFactor = curSampleRate/origSampleRate;
    
    %store data from each miniblock
    R{f} = data_h;
    
    %concatenate events
    for k=1:length(events_h)
        ix = ix + 1;
        event(ix) = events_h(k);
        event(ix).sample = round([events_h(k).sample]*compressionFactor) + nPreviousSamples;
    end
    
    %number of samples already processed earlier in the loop
    nPreviousSamples = nPreviousSamples + hdr_h.nSamples;
end


    res.continuous_data_bpf_rs = (horzcat(R{:}))'; %concatenate data across miniblocks
    res.event = event;
    res.hdr = hdr_h;
    
    res.hdr.nSamples = size(res.continuous_data_bpf_rs,1);
    
    if ~exist(par.concatRunsDir, 'dir')
        mkdir(par.concatRunsDir)
    end
    
    %save concatenated data.
    %save(fullfile(par.concatRunsDir, ['MEGDatFiltRSConcat', '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '.mat']))
    thisFile = fullfile(par.concatRunsDir,['MEGDatFiltRSConcat', '_HP' num2str(par.HPFreq) '_LP' num2str(par.LPFreq) '.mat']);
    ft_write_data(thisFile, res.continuous_data_bpf_rs, 'header', res.hdr, 'dataformat', 'fcdc_matbin');
    ft_write_event(thisFile, res.event)
end


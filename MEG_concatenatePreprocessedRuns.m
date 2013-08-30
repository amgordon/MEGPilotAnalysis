function [ output_args ] = MEG_concatenatePreprocessedRuns(par)

nPreviousSamples = 0;
ix = 0;
figure;
for f = 1:length(par.dataRunsPreproc)
    
    disp(sprintf('loading run %g', f))
    R_h = load(par.dataRunsPreproc{f});
    
    R{f} = R_h.res.continuous_data_bpf_rs.trial{1};
    
    for k=1:length(R_h.res.event)
        ix = ix + 1;
        event(ix) = R_h.res.event(k);
        event(ix).sample = R_h.res.event(k).sample + nPreviousSamples;
    end
    
    nPreviousSamples = nPreviousSamples + length(R_h.res.continuous_data_bpf_rs.time{1});
    
    q = [R_h.res.event.sample];
    clear t
    t(q) = 1;
    plot(t)
end

    res = R_h.res;
    res.continuous_data_bpf_rs = horzcat(R{:});
    res.event = event;
end


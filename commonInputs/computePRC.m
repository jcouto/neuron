files = listFiles('*prc*.h5');

trials = cell(length(files),1);
for f = 1:length(files)
    [data,~,~] = loadCommonInputFile(files{f});
    pert = data.pert;
    spk = data.spikes;
    trial = {};
    for ii = 1:length(pert)
        idx = find(spk<pert(ii),3,'last');
        idx = vertcat(idx,find(spk>pert(ii),3,'first'));
        if (length(idx)==6)
            trial = [trial,spk(idx)-pert(ii)];
        end
    end
    trials{f} = trial;
end
expName = 'Khaliq-Raman with Ek set to -40mV';
pertArea = cellfun(@(x)1,files,'uniformoutput',0);
folders = cellfun(@(x)x(1:end-3),files,'uniformoutput',0);
save('prc_spiketimes.mat','trials','expName','folders','pertArea')
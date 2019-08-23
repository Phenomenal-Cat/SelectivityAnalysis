function Spikes = SF2_GetSpikesByStim(NeuroStruct, StimTimes, TimeWindow)


for cell = 1:size(NeuroStruct.cells, 1)
    AllSpikes = NeuroStruct.cells{cell, 3}/1000;
    for stim = 1:numel(StimTimes)
        for rep = 1:numel(StimTimes{stim})
            TimeRange = TimeWindow + StimTimes{stim}(rep);
            SelectedSpikes = AllSpikes(AllSpikes >= TimeRange(1) & AllSpikes <= TimeRange(2))- StimTimes{stim}(rep);
            Spikes(cell, stim).SpikeTimes{rep} = SelectedSpikes;
        end
    end
end
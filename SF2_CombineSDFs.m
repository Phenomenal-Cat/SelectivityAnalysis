

SpikeDir    = '/Volumes/Seagate Backup 4/NIH_Neurophys/StereoFaces_2/SpikeTimes/';
Subjects    = {'Spice','StevieRay','Mochi','Wasabi'};
NoChannels  = [64, 128, 128, 64];
ExpNames    = {'FingerPrint','SizeDistance','SizeDistance_Movies','StereoShape'};
HistBins    = -100:1:400;
PlotData    = 1;
KernelWidth = 10;

for S = 1:numel(Subjects)
    for exp = 4:numel(ExpNames)
        SpikeFiles = wildcardsearch(fullfile(SpikeDir, Subjects{S}), ExpNames{exp});
        for f = 1:numel(SpikeFiles)
            load(SpikeFiles{f});
            AllSpikeTimes = cell(size(SpikeTimes,1), size(SpikeTimes,2));
            
            if PlotData == 1
                fh(f)   = figure;
                axh     = tight_subplot(10, 10, 0.02, 0.02, 0.02);
                axIndx  = 1;
            end
            for c = 1:size(SpikeTimes,1)
                for stim = 1:size(SpikeTimes,2)
                    for rep = 1:numel(SpikeTimes(c,stim).SpikeTimes)
                        AllSpikeTimes{c,stim} = [AllSpikeTimes{c,stim}; SpikeTimes(c,stim).SpikeTimes{rep}];
                    end
                    PSTH(c,stim,:)  = hist(AllSpikeTimes{c,stim}*1000, HistBins);
                    SDF(c,stim,:)   = msdf(squeeze(PSTH(c,stim,:)), 'Gauss', KernelWidth); 
                end
                if PlotData == 1
                    axes(axh(axIndx));
                    imagesc(HistBins, 1:size(SpikeTimes,2), squeeze(SDF(c,:,:)));
                    axis tight;
                    box off
                    hold on;
                    plot([0,0],ylim,'-w');
                    set(gca,'xtick',[],'ytick',[]);
                    if mod(axIndx,10) == 1
                        set(gca,'ytick', 10:50:size(SpikeTimes,2));
                    end
                    if axIndx >= 91 
                        set(gca,'xtick', -100:100:400);
                    end
%                     xlabel('Time (ms)','fontsize',18);
%                     ylabel('Stim #','fontsize',18);
%                     title(sprintf('Neuron %d', c), 'fontsize',18);
                    axIndx = axIndx+1;
                    drawnow;
                end
            end
            pause(1);


            
        end
    end
end




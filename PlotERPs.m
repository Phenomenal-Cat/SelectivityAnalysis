
%================================ PlotERPs.m ==============================

Subject     = 'Matcha';
Date        = '20160613';
ExpType     = 'StereoFaces';

TimingDir   = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/Timing';
TimingFile  = fullfile(TimingDir, ExpType, sprintf('StimTimes_%s_%s.mat', Subject, Date));

OutputDir   = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/LFP/';
LFPmatfile  = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/RawLFP/Matcha/20160613/Matcha_Matcha_20160613_StereoFaces_01_Stp1_ch1.mat';
load(LFPmatfile);
load(TimingFile);

Channel = 1;
PreTime     = 0.1;
PostTime    = 0.5;
StimDur     = 0.3;
NoiseThresh = 2000;


%=========== Downsample raw LFP if necessary
if fs >= 24414
    LFPdownsampled = decimate(double(rawdata), 24);
    fs = fs/24;
end
TimeStamps      = linspace(0, numel(LFPdownsampled)/fs, numel(LFPdownsampled));
WinSizeSamples 	= round((PreTime+PostTime)*fs);
WinTimes        = linspace(-PreTime, PostTime, WinSizeSamples);


Fh(1)      	= figure('position',get(0,'ScreenSize'));
axh{1}     	= tight_subplot(8, 8, 0.05, 0.05, 0.05);
FigCount    = 1;

%=========== Get LFP samples for each stimulus presentation
for s = 83:numel(Stim.Onsets) 
    ExcludeTrial = zeros(1, numel(Stim.Onsets{s}));
    for t = 1:numel(Stim.Onsets{s})
        WinStartIndx        = find(TimeStamps >= Stim.Onsets{s}(t)-PreTime);
        SampleIndx          = (1:WinSizeSamples)+WinStartIndx(1);
        LFPtrials{s}(t,:)   = LFPdownsampled(SampleIndx);
        if any(find(LFPtrials{s}(t,:) > NoiseThresh))
            ExcludeTrial(t) = 1;
        end
    end
    IncludeTrials{s}= find(~ExcludeTrial);
    LFPmean{s}      = mean(LFPtrials{s}(IncludeTrials{s},:));
    LFPse{s}        = std(LFPtrials{s}(IncludeTrials{s},:))/sqrt(numel(IncludeTrials{s}));
    
    %============= Plot ERP data
    AxIndx = s-(64*(FigCount-1));
    if AxIndx > numel(axh{FigCount})
        set(axh{FigCount}, 'xlim', [-PreTime, PostTime], 'ylim', [-200, 200]);
        suptitle(sprintf('%s %s %s (%d)', Subject, Date, ExpType, FigCount));
        FigName 
        FigCount        = FigCount+1;
        AxIndx          = s-(64*(FigCount-1));
        Fh(FigCount)    = figure('position',get(0,'ScreenSize'));
        axh{FigCount}   = tight_subplot(8, 8, 0.05, 0.05, 0.05);
    end
	
  	axes(axh{FigCount}(AxIndx));
    ph1         = shadedplot(WinTimes, LFPmean{s}-LFPse{s}, LFPmean{s}+LFPse{s}, [1 0.5 0.5], 'r');
    hold on;
    ph2         = plot(WinTimes, LFPmean{s},'-r','linewidth',2);
    Ylims       = get(gca,'ylim');
    ph3         = patch([0,0,StimDur,StimDur], Ylims([1,2,2,1]), 0, 'facecolor', [0 0 0], 'facealpha', 0.5);
    uistack(ph3, 'bottom');
    grid on;
    title(sprintf('Cond %d (n=%d)' ,s, numel(IncludeTrials{s})));
    drawnow;
   
end
set(axh{FigCount}, 'xlim', [-PreTime, PostTime], 'ylim', [-200, 200]);
suptitle(sprintf('%s %s %s (%d)', Subject, Date, ExpType, FigCount));

Matfile = fullfile(OutputDir, Subject, Date, sprintf('LFPproc_%s_%s_%s_ch%d.mat', Subject, Date, ExpType, Channel));
save(Matfile, 'LFPtrials',  'IncludeTrials');



%============= PLOT RESULTS BY FACTOR
Factors     = {'Elevations','Azimuths','Distances','Scales','Expressions','Identity'};
CondMatCol  = [3, 2, 4, 5, 1, 1];
Fh          = figure('position',get(0,'ScreenSize'));
axh         = tight_subplot(2, 3, 0.05, 0.05, 0.05);

for f = 1:numel(Factors)
    Factor      = eval(sprintf('Params.%s', Factors{f}));
    Colors      = jet(numel(Factor));
    axes(axh(f));

    for el = 1:numel(Factor)
        LegendText{el}  = sprintf('%d deg',Factor(el));
        CondIndx        = find(Params.ConditionMatrix(:,CondMatCol(f))==el);
        ElLFPall{el}    = [];
        for c = 1:numel(CondIndx)
            ElLFPall{el} = [ElLFPall{el}; LFPtrials{CondIndx(c)}(IncludeTrials{CondIndx(c)},:)];
        end
        ElLFPmeans{el}  = mean(ElLFPall{el});
        ElLFPse{el}     = std(ElLFPall{el})/sqrt(numel(CondIndx));
    % 	ph1{el} = shadedplot(WinTimes, ElLFPmeans{el}-ElLFPse{el}, ElLFPmeans{el}+ElLFPse{el}, Colors(el,:), 'color', Colors(el,:));
        hold on;
        ph2{el}	= plot(WinTimes, ElLFPmeans{el},'-r','linewidth',2, 'color', Colors(el,:));
    end
    legend(LegendText, 'location', 'northwest', 'fontsize',18);
    axis tight
    ph3 = patch([0,0,StimDur,StimDur], Ylims([1,2,2,1]), 0, 'facecolor', [0 0 0], 'facealpha', 0.5, 'edgecolor','none');
    uistack(ph3, 'bottom');
    grid on;
    drawnow;
    xlabel('Time (s)', 'fontsize', 16);
    ylabel('Voltage (V)', 'fontsize', 16);
    set(gca,'ylim', [-200, 200]);
    title(sprintf('%s', Factors{f}), 'fontsize', 18);
    
end
suptitle(sprintf('%s %s %s', Subject, Date, ExpType));
export_fig(sprintf('ERPs_%s_%s_%s.png', Subject, Date, ExpType), '-png', '-transparent');



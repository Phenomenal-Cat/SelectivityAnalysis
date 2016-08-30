function PlotERPs(LFPfiles, TimingFile, CondFile)

%================================ PlotERPs.m ==============================
% This function loads raw LFP data form the StereoFaces experiments and
% plots ERP and power bands.
%
% INPUTS:   LFPfiles:   a cell array of full path strings for all channels
%                       from the single session to be analysed.
%           TimingFile: full path of the timing info file (.mat)
%           CondFile:   full path of the condition info file (.mat)
%
% HISTORY:
%   08/29/2016 - Updated to load all channels from a session
%==========================================================================

if nargin == 0
	RootDir     = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/LFP_raw/';
    LFPfiles    = uipickfiles('FilterSpec',RootDir, 'REFilter','\.mat$','Prompt','Select LFP data files to analyse');
    TimingFile  = uipickfiles('FilterSpec',RootDir, 'REFilter','\.mat$','Prompt','Select timing file for this session');
    CondFile    = [];
end

addpath(genpath('/Users/aidanmurphy/Documents/Toolboxes/chronux_2_12'));

% TimingFile  = fullfile(RootDir, 'Timing', ExpType, sprintf('StimTimes_%s_%s.mat', Subject, Date));
% OutputDir   = fullfile(RootDir, 'LFP/');
% % LFPmatfile  = fullfile(RootDir, 'RawLFP/Matcha/20160613/Matcha_Matcha_20160613_StereoFaces_01_Stp1_ch1.mat');
% LFPmatfile = fullfile(RootDir, 'TDT_converted', Subject, Date);
% load(LFPmatfile);
load(TimingFile{1});
OutputDir = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/LFP_proc';
OutputMatfile = fullfile(OutputDir, Subject, Date, sprintf('LFPproc_%s_%s_%s.mat', Subject, Date, ExpType));

%=========== Set parameters
Params.PreTime          = 0.1;           % Time before stimulus onset (seconds)
Params.PostTime         = 0.5;           % Time after stimulus onset (seconds)
Params.StimDur          = 0.3;           % Duration of stimulus presentation (seconds)
Params.NoiseThresh      = 2000;          % Noise threshold (mV)


for Ch = 1:numel(LFPfiles)

    load(LFPfiles{Ch});
    
    %=========== Threshold the raw LFP signal
    Thresh = std(double(rawdata))*2;
    rawdata(rawdata>Thresh & rawdata<-Thresh) = nan;

    %=========== Downsample raw LFP if necessary
    if fs >= 24414
        LFPdownsampled = decimate(double(rawdata), 24);
        fs = fs/24;
    end
    TimeStamps      = linspace(0, numel(LFPdownsampled)/fs, numel(LFPdownsampled));
    WinSizeSamples 	= round((Params.PreTime+Params.PostTime)*fs);
    WinTimes        = linspace(-Params.PreTime, Params.PostTime, WinSizeSamples);
    movingwin       = round([0.05, 0.025]*params.fs);
    params.fs     	= fs;


    Fh(Ch)      	= figure('position',get(0,'ScreenSize'));
    axh{Ch}     	= tight_subplot(8, 8, 0.05, 0.05, 0.05);
    FigCount        = 1;

    %=========== Get LFP samples for each stimulus presentation
    for s = 1:numel(Stim.Onsets)
        ExcludeTrial = zeros(1, numel(Stim.Onsets{s}));
        for t = 1:numel(Stim.Onsets{s})
            WinStartIndx        = find(TimeStamps >= Stim.Onsets{s}(t)-Params.PreTime);
            SampleIndx          = (1:WinSizeSamples)+WinStartIndx(1);
            LFPtrials{Ch,s}(t,:)   = LFPdownsampled(SampleIndx);
            if any(find(LFPtrials{Ch,s}(t,:) > Params.NoiseThresh))
                ExcludeTrial(t) = 1;
            end
        end
        IncludeTrials{Ch,s}= find(~ExcludeTrial);
        LFPmean{Ch,s}      = mean(LFPtrials{Ch,s}(IncludeTrials{Ch,s},:));
        LFPse{Ch,s}        = std(LFPtrials{Ch,s}(IncludeTrials{Ch,s},:))/sqrt(numel(IncludeTrials{Ch,s}));

        %============= Plot ERP data
        AxIndx = s-(64*(FigCount-1));
        if AxIndx > numel(axh{FigCount})
            set(axh{FigCount}, 'xlim', [-Params.PreTime, Params.PostTime], 'ylim', [-200, 200]);
            suptitle(sprintf('%s %s %s (%d)', Subject, Date, ExpType, FigCount));
            FigCount        = FigCount+1;
            AxIndx          = s-(64*(FigCount-1));
            Fh(FigCount)    = figure('position',get(0,'ScreenSize'));
            axh{FigCount}   = tight_subplot(8, 8, 0.05, 0.05, 0.05);
        end

        axes(axh{FigCount}(AxIndx));
        ph1         = shadedplot(WinTimes, LFPmean{Ch,s}-LFPse{Ch,s}, LFPmean{Ch,s}+LFPse{Ch,s}, [1 0.5 0.5], 'r');
        hold on;
        ph2         = plot(WinTimes, LFPmean{Ch,s},'-r','linewidth',2);
        Ylims       = get(gca,'ylim');
        ph3         = patch([0,0,Params.StimDur,Params.StimDur], Ylims([1,2,2,1]), 0, 'facecolor', [0 0 0], 'facealpha', 0.5);
        uistack(ph3, 'bottom');
        grid on;
        title(sprintf('Cond %d (n=%d)' ,s, numel(IncludeTrials{Ch,s})));
        drawnow;

    end
    set(axh{FigCount}, 'xlim', [-Params.PreTime, Params.PostTime], 'ylim', [-200, 200]);
    suptitle(sprintf('%s %s %s (%d)', Subject, Date, ExpType, FigCount));

    if ~exist(OutputMatfile,'file')
        save(OutputMatfile, 'LFPtrials',  'IncludeTrials', 'Params');
    else
        save(OutputMatfile, '-append','LFPtrials', 'IncludeTrials', 'Params');
    end
end
    

%============= PLOT RESULTS BY FACTOR
Params.Factors      = {'Elevations','Azimuths','Distances','Scales','Expressions'};     % All factors tested
Params.CondMatCol   = [3, 2, 4, 5, 1, 1];                                             	% Which column is each factor coded in?
Fhf                 = figure('position',get(0,'ScreenSize'));                                   
axh                 = tight_subplot(2, 3, 0.05, 0.05, 0.05);

for f = 1:numel(Params.Factors)
    Factor{f}  	= eval(sprintf('Params.%s', Params.Factors{f}));
    Colors      = jet(numel(Factor{f}));
    axes(axh(f));

    for el = 1:numel(Factor{f})
        LegendText{el}  = sprintf('%d deg',Params.Factor{f}(el));
        CondIndx        = find(Params.ConditionMatrix(:,CondMatCol(f))==el);
        ElLFPall{f,el} 	= [];
        for c = 1:numel(CondIndx)
            ElLFPall{f,el} = [ElLFPall{f,el}; LFPtrials{CondIndx(c)}(IncludeTrials{CondIndx(c)},:)];
        end
        ElLFPmeans{f,el}  = mean(ElLFPall{f,el});
        ElLFPse{f,el}     = std(ElLFPall{f,el})/sqrt(numel(CondIndx));

    % 	ph1{el} = shadedplot(WinTimes, ElLFPmeans{el}-ElLFPse{el}, ElLFPmeans{el}+ElLFPse{el}, Colors(el,:), 'color', Colors(el,:));
        hold on;
        ph3{el}	= plot(WinTimes, ElLFPmeans{f,el},'-r','linewidth',2, 'color', Colors(el,:));
    end
    legend(LegendText, 'location', 'northwest', 'fontsize',18);
    axis tight
    ph4 = patch([0,0,Params.StimDur,Params.StimDur], Ylims([1,2,2,1]), 0, 'facecolor', [0 0 0], 'facealpha', 0.5, 'edgecolor','none');
    uistack(ph4, 'bottom');
    grid on;
    drawnow;
    xlabel('Time (s)', 'fontsize', 16);
    ylabel('Voltage (V)', 'fontsize', 16);
    set(gca,'ylim', [-200, 200]);
    title(sprintf('%s', Factors{f}), 'fontsize', 18);
    
end
suptitle(sprintf('%s %s %s', Subject, Date, ExpType));
export_fig(sprintf('ERPs_%s_%s_%s.png', Subject, Date, ExpType), '-png', '-transparent');

save(Matfile, '-append', 'ElLFPall', 'Factors', 'Factor', 'Params');
% function MF3D_PlotPSTH(SpikeFiles, TimingFile, CondFile)

%============================= MF3D_PlotPSTHs.m ===========================
% This function loads WaveClus spike sorted data form the StereoFaces 
% experiments and plots post-stimulus time histograms/ spike density
% functions for each independent variable.
%
% INPUTS:   SpikeFiles: a cell array of full path strings for all channels
%                       from the single session to be analysed.
%           TimingFile: full path of the timing info file (.mat)
%           CondFile:   full path of the condition info file (.mat)
%
% HISTORY:
%   03/21/2017 - Written by APM
%==========================================================================


%============= LOAD DATA
% if nargin == 0
    Subject = 'Spice';
    Date    = '20160621';
    Append  = [];
    if ismac, Append = '/Volumes'; end
    addpath(genpath('/projects/murphya/APMSubfunctions'))
    %SpikeFiles 	= fullfile(Append, '/NIF/procdata/murphya/Physio/Matcha/20160613/20160613_sorted.mat');
    ProcDataDir = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs/', Subject, Date);
    SpikeFiles 	= fullfile(Append, '/procdata/murphya/Physio/WaveClusSorted/', Subject, sprintf('%s_sorted.mat', Date));
    TimingFile  = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/', sprintf('StimTimes_%s_%s.mat', Subject, Date));
    SaveFigDir  = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs',Subject,Date);
% end
if ~exist(ProcDataDir, 'dir')
    mkdir(ProcDataDir);
end
load(SpikeFiles)
load(TimingFile)

ExpNameString = 'StereoFaces';
if isempty(strfind(NeuroStruct(1).block, ExpNameString))
    error('The data file was not labelled as a ''%s'' experiment block!', ExpNameString);
end

%============ Concatenate multi-block sessions
if numel(NeuroStruct) > 1
    fprintf('Concatenating NeuroStruct blocks:%s\n', NeuroStruct.block);
    AllNeuroStruct = NeuroStruct(1);
    for n = 2:numel(NeuroStruct)
        for c = 1:size(NeuroStruct(n).cells, 1)
            AllNeuroStruct.cells{c, 3} = [AllNeuroStruct.cells{c, 3}; NeuroStruct(n).cells{c , 3}+AllNeuroStruct.blocklength];
        end
        AllNeuroStruct.blocklength = AllNeuroStruct.blocklength + NeuroStruct(n).blocklength;
    end
    NeuroStruct = AllNeuroStruct;
end
NoChannels = NeuroStruct.cells{end,1};

%============ Plot session summary
figure('units','normalized','position', [0,0,1,0.5]);
SpsThresh   = 2;
for c = 1:size(NeuroStruct.cells, 1)
    NoSpikes(c) = numel(NeuroStruct.cells{c,3});
end
SpikesPerSecond = NoSpikes/(NeuroStruct.blocklength/10^3);
bar(SpikesPerSecond);
grid on
set(gca,'tickdir','out', 'fontsize', 18, 'xtick', 1:5:numel(SpikesPerSecond));
hold on;
plot(xlim, [SpsThresh, SpsThresh], '--r');
xlabel('Cell number','fontsize', 18);
ylabel('Mean firing rate (Hz)','fontsize', 18);
title(sprintf('Summary for %s %s %s', Subject, Date, ExpNameString), 'fontsize', 20);
% saveas(gcf, fullfile(ProcDataDir, sprintf('Summary_%s_%s_%s.fig', Subject, Date, ExpNameString)));
% export_fig(fullfile(ProcDataDir, sprintf('Summary_%s_%s_%s.png', Subject, Date, ExpNameString)), '-png');
[~,CellOrder] = sort(SpikesPerSecond,'descend');

%============= 
SaveFigs        = 1;
OverWrite       = 1;        % 0 = skip plots that already exist; 1 = overwrite existing plots
Params.PreTime  = 0.1;      % Time prior to stimulus onset to include (seconds)
Params.PostTime = 0.4;      % Time after stimulus onset to include (seconds)
AxW             = 10;
AxH             = 10;
PlotColors      = [1 0 0; 0 0 1; 0 1 0; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 0.5 0.5 0.5];
SDFline         = 1;
Xlims       = [-100, 400]; 
Xticks      = [-100:100:400];
Ylims       = [0, 30];
StimOn      = [0, 300];
StimWindow  = [-100, 400];                         % Response window (milliseconds)
BinWidth    = 10;                                 % Histogram bin width (miliseconds)
HistBins    = linspace(StimWindow(1), StimWindow(2), diff(StimWindow)/BinWidth);


%% ====================== Loop through all cells
Figname     = sprintf('%s_%s', Subject, Date);
FigTitle    = Figname;
FigTitle(strfind(Figname,'_')) = ' ';
Fig.H       = figure('position',get(0,'ScreenSize')./[1 1 1 1], 'name',sprintf('%s', ExpNameString),'renderer','painters');
Fig.axh     = tight_subplot(AxW, AxH, 0.02, 0.04, 0.04);
                    
wbh = waitbar(0,'');
for ch = 1:NoChannels
    ChIndx = find(cell2mat(NeuroStruct.cells(:,1))==ch);
    
    line = 0;
    for cellno = 1:numel(ChIndx)
        cell = ChIndx(cellno);
        
        if ishandle(wbh)
            waitbar(cellno/size(NeuroStruct.cells,1), wbh, sprintf('Averaging spikes for channel %d of %d...', ch, size(NeuroStruct.cells,1)));
        end


        %========== Loop through all stimuli and repetitions
        AllSpikes = {};
        for n = 1:numel(PD.OnTimes)
            WinStart        = (PD.OnTimes(n)-Params.PreTime)*10^3;
            WinEnd          = (PD.OnTimes(n)+Params.PostTime)*10^3;
            SpikeIndx       = find(NeuroStruct.cells{cell,3}>WinStart & NeuroStruct.cells{cell,3}<WinEnd);
            if ~isempty(SpikeIndx)
                AllSpikes{cell,n}	= NeuroStruct.cells{cell,3}(SpikeIndx)-PD.OnTimes(n)*10^3;
            else
                AllSpikes{cell,n}	= NaN;
            end
        end

        %============= Plot raster
        AxIndx = ch;
        axes(Fig.axh(AxIndx));                             
        for n = 1:size(AllSpikes,2)                                                                 % For each repetition/ trial...
            line = line+1;
            for sp = 1:numel(AllSpikes{cell, n})                                                 % For each spike...
                ph(n,sp) = plot(repmat(AllSpikes{cell, n}(sp), [1,2]), [line-1, line], '-k', 'linewidth', 2,'color', PlotColors(cellno,:));       % Draw a vertical line
                hold on;
            end
        end
        if cellno == numel(ChIndx)
            axis tight off
            mkh     = plot([0 0], ylim, '-b', 'linewidth', 2);
            AxesPos = get(Fig.axh(AxIndx+AxW), 'position');
            set(Fig.axh(AxIndx), 'position', [AxesPos(1), sum(AxesPos([2,4])), AxesPos(3), AxesPos(4)], 'xlim', Xlims);
            title(sprintf('Ch %d', ch), 'fontsize', 12);
        end

        %============= Plot PSTH / SDF
        axes(Fig.axh(AxIndx+AxW));
        for n = 1:size(AllSpikes,2)
            BinData{cell}(t,:) = (hist(AllSpikes{cell, n}, HistBins))*(1000/BinWidth);
        end
        BinMeans{cell}	= mean(BinData{cell});
        BinSEM{cell}    = std(BinData{cell})/sqrt(size(BinData{cell},1));
    %             BinData{cell, s}  = (hist(SpikeTimes{cell, s}(:), HistBins))/BinWidth/size(SpikeTimes{cell, s},1);
        if SDFline == 1
            [ha, hb, hc] = shadedplot(HistBins, BinMeans{cell}-BinSEM{cell}, BinMeans{cell}+BinSEM{cell}, PlotColors(end,:));
            hold on;
            delete([hb, hc]);
            plot(HistBins, BinMeans{cell}, '-b', 'color', PlotColors(cellno,:), 'linewidth', 2);
        else
            bar(HistBins, BinMeans{cell});
            hold on;
            errorbar(HistBins, BinMeans{cell}, BinSEM{cell}, '.k');
        end

    end
  	ph = patch(StimOn([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'none', 'facealpha', 0.5);
    uistack(ph, 'bottom')
    axis tight
    set(gca, 'tickdir', 'out', 'xlim', Xlims, 'xtick', Xticks, 'ylim', Ylims)
    box off
    grid on
    if mod(ch,AxW) == 1
        ylabel('Firing rate (Hz)');
    else
        set(gca, 'yticklabel', []);
    end
    if ch > AxW*(AxH-1)
       xlabel('Time (ms)'); 
    end
    drawnow
    
    
end

suptitle(FigTitle);
% export_fig([fullfile(SaveFigDir, Figname),'.png'],'-png'); 

% suptitle(sprintf('%s %s %s channel %d cell %d', Subject, Date, ExpNameString, NeuroStruct.cells{cell,1}, NeuroStruct.cells{cell,2}));
% saveas(Fhf, fullfile(ProcDataDir, sprintf('ERPs_%s_%s_%s_ch%d_cell%d.fig', Subject, Date, ExpNameString, NeuroStruct.cells{cell,1}, NeuroStruct.cells{cell,2})));
% export_fig(fullfile(ProcDataDir, sprintf('ERPs_%s_%s_%s_ch%d_cell%d.png', Subject, Date, ExpNameString, NeuroStruct.cells{cell,1}, NeuroStruct.cells{cell,2})), '-png', '-transparent');

delete(wbh)
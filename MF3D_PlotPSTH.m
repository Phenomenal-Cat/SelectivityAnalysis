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
    Date    = '20160620';
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
            AllNeuroStruct.cells{c, 3} = [AllNeuroStruct.cells{c, 3}, NeuroStruct(n).cells{c , 3}+AllNeuroStruct.blocklength];
        end
        AllNeuroStruct.blocklength = AllNeuroStruct.blocklength + NeuroStruct(n).blocklength;
    end
    NeuroStruct = AllNeuroStruct;
end

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
saveas(gcf, fullfile(ProcDataDir, sprintf('Summary_%s_%s_%s.fig', Subject, Date, ExpNameString)));
[~,CellOrder] = sort(SpikesPerSecond,'descend');

%============= 
SaveFigs        = 1;
Params.PreTime  = 0.1;      % Time prior to stimulus onset to include (seconds)
Params.PostTime = 0.4;      % Time after stimulus onset to include (seconds)
AxW             = 10;
AxH             = 10;
PlotColors      = [1 0 0; 1 0.5,0.5];
SDFline         = 1;
Xlims       = [-100, 400]; 
Xticks      = [-100:100:400];
Ylims       = [0, 100];
StimOn      = [0, 300];
StimWindow  = [-100, 400];                         % Response window (milliseconds)
BinWidth    = 10;                                 % Histogram bin width (miliseconds)
HistBins    = linspace(StimWindow(1), StimWindow(2), diff(StimWindow)/BinWidth);


%============= Loop through all cells
wbh = waitbar(0,'');
for cell = 1:size(NeuroStruct.cells,1)
    waitbar(cell/size(NeuroStruct.cells,1), wbh, sprintf('Arranging spikes for cell %d of %d...', cell, size(NeuroStruct.cells,1)));
    Fig.Current     = 0;
    
    %=========== Loop through all stimuli
    for s = 1:numel(Stim.Onsets)
        
        if mod(s, AxW*AxH/2) == 1
           	if Fig.Current > 0
                Figname = sprintf('%s_%s_cell%d_plot%d', Subject, Date, cell, Fig.Current);
                FigTitle = Figname;
                FigTitle(strfind(Figname,'_')) = ' ';
                suptitle(FigTitle);
                %saveas(Fig.H(cell, Fig.Current), [fullfile(SaveFigDir, Figname),'.fig'],'fig');
                export_fig([fullfile(SaveFigDir, Figname),'.png'],'-png');
            end
            Fig.Current                 = Fig.Current+1;
          	Fig.H(cell, Fig.Current)  	= figure('position',get(0,'ScreenSize')./[1 1 1 1], 'name',sprintf('%s cell %d', ExpNameString, cell),'renderer','painters');
            Fig.axh{cell, Fig.Current}	= tight_subplot(AxW, AxH, 0.02, 0.04, 0.04);
        end
        
        %========== Loop through all stimulus repetitions
        for t = 1:numel(Stim.Onsets{s})
            WinStart        = (Stim.Onsets{s}(t)-Params.PreTime)*10^3;
            WinEnd          = (Stim.Onsets{s}(t)+Params.PostTime)*10^3;
            SpikeIndx       = find(NeuroStruct.cells{cell,3}>WinStart & NeuroStruct.cells{cell,3}<WinEnd);
            if ~isempty(SpikeIndx)
                AllSpikes{cell,s,t}	= NeuroStruct.cells{cell,3}(SpikeIndx)-Stim.Onsets{s}(t)*10^3;
            else
                AllSpikes{cell,s,t}	= NaN;
            end
        end

        %============= Plot raster
        AxPerFig = s-(Fig.Current-1)*AxW*AxH/2;
        AxIndx = floor((AxPerFig-1)/AxW)*AxW*2 +mod(AxPerFig-1, AxW)+1;
        Fig.RasterAx{cell}(s) = Fig.axh{cell, Fig.Current}(AxIndx);
        axes(Fig.axh{cell, Fig.Current}(AxIndx));                             
        line = 1;
        for t = 1:size(AllSpikes,3)                                                                 % For each repetition/ trial...
            line = line+1;
            for sp = 1:numel(AllSpikes{cell, s, t})                                                 % For each spike...
                ph(t,sp) = plot(repmat(AllSpikes{cell, s, t}(sp), [1,2]), [line-1, line], '-k');  	% Draw a vertical line
                hold on;
            end
        end
        axis tight off
        mkh(s)  = plot([0 0], ylim, '-b', 'linewidth', 2);
        AxesPos = get(Fig.axh{cell, Fig.Current}(AxIndx+AxW), 'position');
        set(Fig.axh{cell, Fig.Current}(AxIndx), 'position', [AxesPos(1), sum(AxesPos([2,4])), AxesPos(3), AxesPos(4)/2], 'xlim', Xlims);
        title(sprintf('Stim %d', s), 'fontsize', 12);

        %============= Plot PSTH / SDF
        Fig.SDFAx{cell}(s) = Fig.axh{cell}(AxIndx+AxW);
        axes(Fig.axh{cell, Fig.Current}(AxIndx+AxW));
        for t = 1:size(AllSpikes, 3)
            BinData{cell, s}(t,:) = (hist(AllSpikes{cell, s, t}, HistBins))*(1000/BinWidth);
        end
        BinMeans{cell, s}	= mean(BinData{cell, s});
        BinSEM{cell, s}     = std(BinData{cell, s})/sqrt(size(BinData{cell, s},1));
%             BinData{cell, s}  = (hist(SpikeTimes{cell, s}(:), HistBins))/BinWidth/size(SpikeTimes{cell, s},1);
        if SDFline == 1
            [ha, hb, hc] = shadedplot(HistBins, BinMeans{cell, s}-BinSEM{cell, s}, BinMeans{cell, s}+BinSEM{cell, s}, PlotColors(2,:));
            hold on;
            delete([hb, hc]);
            plot(HistBins, BinMeans{cell, s}, '-b', 'color', PlotColors(1,:), 'linewidth', 2);
        else
            bar(HistBins, BinMeans{cell, s});
            hold on;
            errorbar(HistBins, BinMeans{cell, s}, BinSEM{cell, s, sz}, '.k');
        end
        ph = patch(StimOn([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'none', 'facealpha', 0.5);
        uistack(ph, 'bottom')
        axis tight
%         set(gca, 'tickdir', 'out', 'xlim', Xlims, 'xtick', Xticks, 'ylim', Ylims)
        box off

        if mod(s,AxW) == 1
            ylabel('Firing rate (Hz)');
        else
            set(gca, 'yticklabel', []);
        end
        if s > AxW*(AxH-1)
           xlabel('Time (ms)'); 
        end
        drawnow
        
    end
    close all
    


    %% ===================== PLOT RESULTS BY FACTOR =======================
    Params.Factors      = {'Elevations','Azimuths','Distances','Scales','Expressions'};     % All factors tested
    Params.CondMatCol   = [3, 2, 4, 5, 1, 1];                                             	% Which column is each factor coded in?
    Fhf                 = figure('position',get(0,'ScreenSize'),'renderer','painters');                                   
    axh                 = tight_subplot(2, 3, 0.05, 0.05, 0.05);

    for f = 1:numel(Params.Factors)-1
        Factor{f}  	= eval(sprintf('Params.%s', Params.Factors{f}));
        Colors      = jet(numel(Factor{f}));
        axes(axh(f));

        for el = 1:numel(Factor{f})
            LegendText{el}  = sprintf('%d deg', Factor{f}(el));
            CondIndx        = find(Params.ConditionMatrix(:,Params.CondMatCol(f))==el);
            ElLFPall{f,el} 	= [];
            for c = 1:numel(CondIndx)
                ElLFPall{f,el} = [ElLFPall{f,el}; BinData{cell, CondIndx(c)}];
            end
            ElLFPmeans{f,el}  = mean(ElLFPall{f,el});
            ElLFPse{f,el}     = std(ElLFPall{f,el})/sqrt(numel(CondIndx));

        % 	ph1{el} = shadedplot(WinTimes, ElLFPmeans{el}-ElLFPse{el}, ElLFPmeans{el}+ElLFPse{el}, Colors(el,:), 'color', Colors(el,:));
            hold on;
            ph3{el}	= plot(HistBins, ElLFPmeans{f,el},'-r','linewidth',2, 'color', Colors(el,:));
        end
        legend(LegendText, 'location', 'northwest', 'fontsize',18);
        axis tight
        ph4 = patch([0,0,StimOn([2,2])], Ylims([1,2,2,1]), 0, 'facecolor', [0.5 0.5 0.5], 'facealpha', 1, 'edgecolor','none');
        uistack(ph4, 'bottom');
        grid on;
        xlabel('Time (s)', 'fontsize', 16);
        ylabel('Firing rate (Hz)', 'fontsize', 16);
        title(sprintf('%s', Factor{f}), 'fontsize', 18);
        set(gca,'ytick',Ylims(1):20:Ylims(2),'tickdir','out')
    	drawnow;
    end
    suptitle(sprintf('%s %s %s cell %d', Subject, Date, ExpNameString, cell));
    saveas(Fhf, fullfile(ProcDataDir, sprintf('ERPs_%s_%s_%s_cell%d.fig', Subject, Date, ExpNameString, cell)));
    export_fig(fullfile(ProcDataDir, sprintf('ERPs_%s_%s_%s_cell%d.png', Subject, Date, ExpNameString, cell)), '-png', '-transparent');

end

%================= Save processed data to .mat files
SortedDataFilename = fullfile(ProcDataDir, sprintf('%s_%s.mat', Subject, Date));
save(SortedDataFilename, 'AllSpikes','BinData','BinMeans','BinSEM','HistBins');

save(SortedDataFilename, '-append', 'ElLFPall', 'Factor', 'Params');

%========================= MF3D_PlotWaveforms.m ===========================
% Plot the spike waveforms for all cells on all channels across sessions.
%

Subject     = 'Avalanche';
NoChannels  = 128;

Append = [];
if ismac, Append = '/Volumes'; end
WaveClusDir = fullfile(Append, '/NIF/procdata/murphya/Physio/StereoFaces/',Subject);
Sessions = dir(WaveClusDir);
AllDates = {Sessions.name};
AllDates(~cellfun(@isempty, strfind(AllDates, '.'))) = [];

SampleRate 	= 24000;
ChPerFig    = 12;
FRthresh    = 2000;                             % Minimum number of spikes for a cell to be plotted
AmpThresh   = 100;                              % Maximum spike amplitude (uV) for inclusion as a real neuron
NoFigs      = ceil(NoChannels/ChPerFig);
FigIndx     = 0;
LineColors  = [0 0 1; 1,0,0;0 1 0; 1 0 1; 1 1 0; 0 1 1];
IgnoreColor = [0.5,0.5,0.5];

for ch = 1:NoChannels
    if mod(ch-1, ChPerFig) == 0
        FigIndx         = FigIndx +1;
        fh(FigIndx)     = figure('position',get(0,'screensize'));
        axh{FigIndx}	= tight_subplot(numel(AllDates), ChPerFig, 0.02, 0.05, 0.05);
    end
    
    for d = 1:numel(AllDates)
        ChannelDir      = fullfile(WaveClusDir, AllDates{d}, sprintf('%d',ch));
        ChannelFile     = wildcardsearch(ChannelDir, 'times*_concat_spikes.mat');
        load(ChannelFile{1});
        Cells           = unique(cluster_class(:,1));
        TimePoints      = linspace(0, size(spikes,2)/SampleRate*1000, size(spikes,2));
        AxIndx          = ch+(d-1)*ChPerFig;
        axes(axh{FigIndx}(AxIndx))
        for c = 1:numel(Cells)
            SpikeIndx{c}        = find(cluster_class(:,1)==Cells(c));
            AllSpikes{ch,d,c}   = spikes(SpikeIndx{c},:);
            WaveMean{ch,d,c}    = mean(AllSpikes{ch,d,c});
            WaveSEM{ch,d,c}     = std(AllSpikes{ch,d,c})/sqrt(size(AllSpikes{ch,d,c},1));
            if max(WaveMean{ch,d,c}) > AmpThresh || numel(SpikeIndx{c}) < FRthresh
                PlotColor = IgnoreColor;
            else
                PlotColor = LineColors(c,:);
            end
            Ph{ch,d}(c)       	= plot(TimePoints, WaveMean{ch,d,c},'-k','color',PlotColor); 
            hold on;
        end
        grid on
        axis tight
        if mod(ch-1, ChPerFig) == 0
            if d == 1
                title(sprintf('%s\tCh %d', AllDates{d}, ch), 'fontsize', 12, 'HorizontalAlignment', 'left');
            else
                title(AllDates{d}, 'fontsize', 12, 'HorizontalAlignment', 'left');
            end
            ylabel('(\muV)', 'fontsize', 12);
        else
            set(gca,'yticklabel',[]);
        end
        if d == 1 && mod(ch-1, ChPerFig)~=0
            title(sprintf('Ch %d', ch), 'fontsize', 12);
        end
        if d < numel(AllDates)
            set(gca,'xticklabel',[]);
        else
            xlabel('Time (ms)', 'fontsize', 12);
        end
        drawnow
    end
    
end
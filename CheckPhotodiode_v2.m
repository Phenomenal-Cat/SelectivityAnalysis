function PD = CheckPhotodiode_v2(Signal, SampleRate, PD, BlockTimes)

%========================== CheckPhotodiode.m =============================
% This function analyses the signal from a photodiode attached to the
% visual stimulus presentation display and extracts the time stamps of
% state changes in photodiode input luminance. Data is plotted to allow the
% user to interactively adjust filtering and threshold parameters.
%
% INPUTS:   Signal:     
%           SampleRate: 
%           PD:         (optional) structure returned from previous call of
%                       CheckPhotodiode.m, with user accepted parameters.
% REVISIONS:
%   01/07/2016 - Written by APM
%     ___  ______  __   __
%    /   ||  __  \|  \ |  \    APM SUBFUNCTIONS
%   / /| || |__/ /|   \|   \   Aidan P. Murphy - murphyap@mail.nih.gov
%  / __  ||  ___/ | |\   |\ \  Section on Cognitive Neurophysiology and Imaging
% /_/  |_||_|     |_| \__| \_\ National Institute of Mental Health
%==========================================================================

global PD Fig

PD.Signal       = Signal;
PD.SampleRate   = SampleRate;
PD.TimeStamps   = linspace(0,numel(PD.Signal)/PD.SampleRate,numel(PD.Signal));
PD.FilterOn     = 0; 
PD.Inverted     = 0;
PD.PlotOn       = 1;

%============== High-pass filter signal
PD.Filt.Cutoff    	= 0.001;  
if PD.FilterOn == 1                                          
    [b,a]               = butter(2, PD.Filt.Cutoff/PD.SampleRate*2,'high');
    PD.SignalFilt       = filtfilt(b,a,PD.Signal);
else
    PD.SignalFilt       = PD.Signal;
    PD.ThreshFilt       = mean(PD.SignalFilt);
end

%============== Threshold signal
if ~isfield(PD,'Thresh')
    PD.Thresh       = mean(PD.Signal);
    PD.ThreshFilt   = mean(PD.SignalFilt);
end
PD.ThreshSignal     = PD.SignalFilt <= PD.ThreshFilt;                                   
PD.OnSamples        = find(diff(PD.ThreshSignal)==1);                        	% Get photodiode onset samples
PD.OffSamples       = find(diff(PD.ThreshSignal)==-1);                       	% Get photodiode offset samples
PD.OnTimes          = PD.TimeStamps(PD.OnSamples);
PD.OffTimes         = PD.TimeStamps(PD.OffSamples);                         	% Get time stamp for samples
if numel(PD.OnTimes) ~= numel(PD.OffTimes) 
    PD.OnOffMatch   = 0;
    PD.Durations   	= NaN;                                                      % Get duration of photodiode responses     
else
    PD.OnOffMatch   = 1;
    PD.Durations   	= PD.OffTimes-PD.OnTimes;                                   % Get duration of photodiode responses     
end
if PD.OnTimes(1) < PD.OffTimes(1)
    PD.OnOffMatch   = 0;   
end
        

%============== Plot data
if PD.PlotOn == 1
    Fig.Pos = get(0, 'ScreenSize');
    Fig.H   = figure('Name','Photodiode signal analysis','position',Fig.Pos);

    Fig.axh(1)      = subplot(3,3,[1,2]);
    Fig.TCraw       = plot(PD.TimeStamps, PD.Signal-PD.Thresh,'-b');                    % Plot raw signal
    hold on;
    Fig.TCfilt      = plot(PD.TimeStamps, PD.SignalFilt,'-r');                          % Plot filtered signal
    Fig.TCthresh(1) = plot(PD.TimeStamps, PD.ThreshSignal*range(PD.SignalFilt),'-g'); 	% plot thresholded signal
    set(Fig.axh(1), 'xlim', [0, PD.TimeStamps(end)]);
    ylabel('Signal (V)','fontsize',16)
    xlabel('Time (seconds)','fontsize',16);
    legend({'Raw','Filtered','Thresholded'},'Location','EastOutside','fontsize',16);
    title('Photodiode timecourse','fontsize',18)
    if exist('BlockTimes','var')
        for b = 1:numel(BlockTimes)
            plot(repmat(BlockTimes(b),[1,2]), ylim, '-m', 'linewidth',2);
        end
    end
    
  	Fig.axh(2) = subplot(3,3,4);
    Fig.HistH = histogram(PD.Signal, 100);
    hold on;
    Fig.ThreshH = plot(repmat(PD.Thresh,[1,2]),ylim, '-r', 'linewidth',2);
    xlabel('Photodiode voltage','fontsize',16);
    ylabel('# samples','fontsize',16);
    title('Original signal','fontsize',18)
    
    Fig.axh(3) = subplot(3,3,5);
  	Fig.HistFiltH = histogram(PD.SignalFilt, 100);
    hold on;
    Fig.ThreshFiltH = plot(repmat(PD.ThreshFilt,[1,2]),ylim, '-r', 'linewidth',2);
    xlabel('Photodiode voltage','fontsize',16);
    ylabel('# samples','fontsize',16);
    title('Filtered signal','fontsize',18)
       
    Fig.axh(4) = subplot(3,3,[7,8]);
    Fig.TCthresh(2) = plot(PD.TimeStamps, PD.ThreshSignal,'-b');
    hold on;
    Fig.PDonsets = plot(PD.OnTimes, ones(size(PD.OnTimes)), '*r');
%     plot(BlockStartTimes, ones(size(BlockStartTimes)), '.g', 'markersize',40);
    set(Fig.axh(4), 'xlim', [0, PD.TimeStamps(end)],'ylim',[-0.1,1.1]);
    xlabel('Time (seconds)','fontsize',16);
    ylabel('State (on/off)','fontsize',16);
    title('Photodiode timecourse','fontsize',18);

    linkaxes(Fig.axh([1,4]),'x');
    set(zoom(Fig.axh(1)),'Motion','horizontal');
    
    Fig.axh(5) = subplot(3,3,9);
    if PD.OnOffMatch == 1
        PD.DurHist = histogram(PD.Durations);
        xlabel('State duration (seconds)');
        ylabel('Frequency');
    else
        set(Fig.axh(5), 'color', [0.5,0.5,0.5]);
    end
    
    
    %==================== Add UI menus
    Fig.UI.PanelPos = [Fig.Pos(3)-300, Fig.Pos(4)-500, 250, 300];
    Fig.UI.PanelH  	= uipanel('Title','Parameters','FontSize', 16,'Units','pixels','Position',Fig.UI.PanelPos);
    Fig.UI.Labels   = {'High pass filter?','Cutoff (Hz)','Raw thresh (V)','Filt Thresh (V)','Accept?','Invert signal?'};
    Fig.UI.Style    = {'Checkbox','Edit', 'Edit','Edit','PushButton','Checkbox'};
    Fig.UI.List     = {'on/off', num2str(PD.Filt.Cutoff), num2str(PD.Thresh),num2str(PD.ThreshFilt),'Accept','on/off'};
    Fig.UI.InputDim = [100, 20];
    for i = 1:numel(Fig.UI.Labels)
        Pos = numel(Fig.UI.Labels)-i;
        Fig.UI.LabelPos{i}      = [10, 10+Pos*(Fig.UI.InputDim(2)+5),Fig.UI.InputDim];
        Fig.UI.LabelHandle(i)   = uicontrol('Style','Text',...
                                            'String',Fig.UI.Labels{i},...
                                            'HorizontalAlignment','Left',...
                                            'pos',Fig.UI.LabelPos{i},...
                                            'parent',Fig.UI.PanelH);
        Fig.UI.InputHandle(i) = uicontrol(  'Style',Fig.UI.Style{i},...
                                            'String',Fig.UI.List{i},...
                                            'HorizontalAlignment','Left',...
                                            'pos',[Fig.UI.InputDim(1)+10,15+Pos*(Fig.UI.InputDim(2)+5),Fig.UI.InputDim],...
                                            'Callback',{@UpdateParams,i},...
                                            'parent',Fig.UI.PanelH);
    end
    set(Fig.UI.InputHandle(1),'value',PD.FilterOn);
    set(Fig.UI.InputHandle(6),'value',PD.Inverted);
end

uiwait(Fig.H);


end

%=================== USER UPDATED PARAMETER
function UpdateParams(hObj, evnt, indx)
global Fig PD

switch indx
    case 1
        PD.FilterOn = get(hObj,'value');
        
    case 2
        PD.Filt.Cutoff = str2double(get(hObj,'string'));
        
    case 3
        PD.Thresh = str2double(get(hObj,'string'));
        
    case 4
        PD.ThreshFilt = str2double(get(hObj,'string'));
        
    case 5
        close(Fig.H);
        return
    case 6
        PD.Inverted = get(hObj,'value');
end
UpdateFig(indx);

end

%=================== UPDATE PLOTS
function UpdateFig(indx)
global Fig PD

switch indx
    case 1
        if PD.FilterOn == 1
            set(Fig.axh(3), 'Color',[1 1 1]);
            set(Fig.TCfilt, 'Color',[1 0 0]);
        elseif PD.FilterOn == 0
            set(Fig.axh(3), 'Color',[0.5 0.5 0.5]);
            set(Fig.TCfilt, 'Color',[0.5 0.5 0.5]);
        end
    case 2

    case 3
        set(Fig.ThreshH, 'xdata', repmat(PD.Thresh,[1,2]));
    case 4
        set(Fig.ThreshFiltH, 'xdata', repmat(PD.ThreshFilt,[1,2]));
    case 6
        
end
UpdateFilt;

end

%% ==================== Update filter
function UpdateFilt
global Fig PD

%========== Filter signal
if PD.FilterOn == 1 
    [b,a]               = butter(2, PD.Filt.Cutoff/PD.SampleRate*2,'high');
    PD.SignalFilt       = filtfilt(b,a,PD.Signal);
    set(Fig.TCfilt, 'ydata', PD.SignalFilt);
end

%========== Use filtered signal?
if PD.FilterOn == 0
    PD.ThreshSignal     = double(PD.Signal <= PD.Thresh);                                        
  	
elseif PD.FilterOn == 1
    PD.ThreshSignal     = double(PD.SignalFilt <= PD.ThreshFilt);                                   
 	set(Fig.TCfilt, 'ydata', PD.SignalFilt);
    delete(Fig.HistFiltH)
    axes(Fig.axh(3));
    Fig.HistFiltH       = histogram(PD.SignalFilt, 100);
end
set(Fig.TCthresh(1), 'ydata', PD.ThreshSignal*range(PD.SignalFilt));
set(Fig.TCthresh(2), 'ydata', PD.ThreshSignal);

%========== Update plots
if PD.Inverted == 0
    PD.OnSamples        = find(diff(PD.ThreshSignal)==1);                        	% Get photodiode onset samples
    PD.OffSamples       = find(diff(PD.ThreshSignal)==-1);                       	% Get photodiode offset samples
else
  	PD.OnSamples        = find(diff(PD.ThreshSignal)==-1);                        	% Get photodiode onset samples
    PD.OffSamples       = find(diff(PD.ThreshSignal)==1);                       	% Get photodiode offset samples
end
PD.OnTimes          = PD.TimeStamps(PD.OnSamples);
PD.OffTimes         = PD.TimeStamps(PD.OffSamples);       
if numel(PD.OnTimes) ~= numel(PD.OffTimes) 
    PD.OnOffMatch   = 0;
    PD.Durations   	= NaN;                                                      % Get duration of photodiode responses     
else
    PD.OnOffMatch   = 1;
    PD.Durations   	= PD.OffTimes-PD.OnTimes;                                   % Get duration of photodiode responses     
end
if isfield(PD,'DurHist')
    delete(PD.DurHist);
end
if PD.OnOffMatch   == 1
    axes(Fig.axh(5));
    set(Fig.axh(5), 'color', [1 1 1]);
    PD.DurHist = histogram(PD.Durations);
else
    set(Fig.axh(5), 'color', [0.5,0.5,0.5]);
end

delete(Fig.PDonsets);
axes(Fig.axh(4));
Fig.PDonsets = plot(PD.OnTimes, ones(size(PD.OnTimes)), '*r');

end

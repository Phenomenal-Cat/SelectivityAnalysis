
%===================== MF3D_PlotResponses_Expressions.m ===================
% Plot facial epxression tuning for STS cells tested in 'StereoFace' pilot
% experiments.
%
%==========================================================================

Subject     = 'Avalanche';
Date        = '20160630';

switch Subject
    case 'Avalanche'
        if ~any(~cellfun(@isempty, strfind({'20160630','20160701'}, Date)))
            error('Invalid session for expression analysis!')
        end
    case 'Matcha'
        if ~any(~cellfun(@isempty, strfind({'20160616','20160617'}, Date)))
            error('Invalid session for expression analysis!')
        end
    case 'Spice'
        if ~any(~cellfun(@isempty, strfind({'20160623','20160624'}, Date)))
            error('Invalid session for expression analysis!')
        end
end

%============= Set directories
if ~exist('AllSpikes','var')
    Append  = [];
    if ismac, Append = '/Volumes'; end
    StimDir                 = fullfile(Append, '/projects/murphya/MacaqueFace3D/BlenderFiles/Renders/Monkey_1/');
    TimingData              = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/',sprintf('StimTimes_%s_%s.mat', Subject, Date));
    ProcessedSessionData    = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs/',Subject,Date,sprintf('%s_%s.mat', Subject, Date));
    HeadOrientationDir      = fullfile(Append, '/projects/murphya/MacaqueFace3D/PilotData/PNGs/');
    load(TimingData)
    load(ProcessedSessionData);
end
SaveDir = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/FacialExpressionTuning', Subject);
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

%=========== Load expression condition images
for exp = 1:numel(Params.Expressions)
    HeadExpImage = fullfile('/Volumes/projects/murphya/MacaqueFace3D/PilotData/PNGs/',sprintf('OrientStim3x3_%s.png',Params.Expressions{exp}));
    [HeadIm,cm,HeadAlpha]   = imread(HeadExpImage);                                   % Load head orientation overlay
    HeadImSize              = size(HeadIm);
    HeadImages{exp}         = HeadIm;
    HeadImAlpha{exp}        = HeadAlpha;
end

%============= Settings
Channel     = 117;
CellIndx   	= 2;
Cell        = find(ismember(ChIndx(:,[2,3]), [Channel, CellIndx], 'rows'));
Twin        = [80, 150];                                                    % Specify time window to calculate mean response from (ms)
BaseWin     = [-50 50];                                                     % Specify time window to calculate baseline response from (ms)
Params.Depths = {'Concave','Flat','Convex'};
Params.DepthIndx    = [-1, 0, 1];
if ~isempty(find(Params.ConditionMatrix(:,6)==-1))
    for p = numel(Params.DepthIndx):-1:1
        Params.ConditionMatrix(Params.ConditionMatrix(:,6)==Params.DepthIndx(p), 6) = p;
    end
end
        
%========== Pool neural data across scales and depths
figure('position', [-1911, 168, 1799, 840]);
Axh     = tight_subplot(6,6,0.02, 0.05,0.05);
AxIndx  = 1;
Ylims   = [0 100];                                                  % Specify 
Tindx   = find(HistBins>Twin(1) & HistBins<Twin(2));
BaseIndx = find(HistBins>BaseWin(1) & HistBins<BaseWin(2));
WinColor = [0.5, 0.5, 0.5];
WinAlpha = 0.5;

ExpColors = [0.75,0.75,1; 1,1,0.75; 0.75,1,0.75; 1, 0.75, 0.75];
ExpAxIndx = [1,2,3,7,8,9,13,14,15];

for exp = 1:numel(Params.Expressions)
    AxIndx = 1;
    if exp == 2 || exp == 4
        ExpAxIndx = ExpAxIndx+3;
    elseif exp == 3
        ExpAxIndx = ExpAxIndx+15;
    end
    for el = 1:numel(Params.Elevations)
        for az = 1:numel(Params.Azimuths)
            CondIndx = find(ismember(Params.ConditionMatrix(:,[1,2,3]), [exp,az,el],'rows'));
            OrientationSDF{exp,az,el} = [];
            for c = 1:numel(CondIndx)
                for t = 1:size(AllSpikes, 3)
                     if ~isnan(AllSpikes{Cell, CondIndx(c), t})
                         OrientationSDF{exp,az,el}(end+1,:) = hist(AllSpikes{Cell, CondIndx(c), t}, HistBins)*10^3/diff(HistBins([1,2]));
                     end
                end
            end
            OrientRawMean(exp,az,el)   = mean(mean(OrientationSDF{exp,az,el}(:,Tindx)));
            OrientBaseMean(exp,az,el)  = mean(mean(OrientationSDF{exp,az,el}(:,BaseIndx)));
            OrientRawSEM{exp}(az,el)    = std(std(OrientationSDF{exp,az,el}(:,Tindx)))/sqrt(size(OrientationSDF{exp,az,el},1));
            OrientDiffMat{exp}(az,el)  = OrientRawMean(exp,az,el)-OrientBaseMean(exp,az,el);
        end
    end
    
    %========= Plot depth data
    for dist = 1:numel(Params.Distances)
        for d = 1:numel(Params.Depths)
            CondIndx = find(ismember(Params.ConditionMatrix(:,[1,4,6]), [exp,dist,d],'rows'));
            StereoSDF{exp,dist,d} = [];
            for c = 1:numel(CondIndx)
                for t = 1:size(AllSpikes, 3)
                     if ~isnan(AllSpikes{Cell, CondIndx(c), t})
                         StereoSDF{exp,dist,d}(end+1,:) = hist(AllSpikes{Cell, CondIndx(c), t}, HistBins)*10^3/diff(HistBins([1,2]));
                     end
                end
            end
        	StereoRawMean(exp,dist,d)   = mean(mean(StereoSDF{exp,dist,d}(:,Tindx)));
            StereoBaseMean(exp,dist,d)  = mean(mean(StereoSDF{exp,dist,d}(:,BaseIndx)));
            StereoRawSEM{exp}(dist,d) 	= std(std(StereoSDF{exp,dist,d}(:,Tindx)))/sqrt(size(StereoSDF{exp,dist,d},1));
            StereoDiffMat{exp}(dist,d)  = StereoRawMean(exp,dist,d)-StereoBaseMean(exp,dist,d);
            
            %========= Plot depth data
            axes(Axh(ExpAxIndx(AxIndx)));
            BinSEM{exp,dist,d} = std(StereoSDF{exp,dist,d})/sqrt(size(StereoSDF{exp,dist,d}, 1));
            [ha, hb, hc] = shadedplot(HistBins, mean(StereoSDF{exp,dist,d})-BinSEM{exp,dist,d}, mean(StereoSDF{exp,dist,d})+BinSEM{exp,dist,d}, [1,0.5,0.5]);
            hold on;
            delete([hb, hc]);
            plot(HistBins, mean(StereoSDF{exp,dist,d}), '-r');
            ph = patch(Twin([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
            uistack(ph, 'bottom')
            
            if ismember(ExpAxIndx(AxIndx), [1:6:(6*6)])
                ylabel(sprintf('Dist = %d cm', Params.Distances(dist)), 'fontsize', 16);
            end
            if ExpAxIndx(AxIndx) >= (5*6)+1
                xlabel(sprintf('%s', Params.Depths{d}), 'fontsize', 16);
            end
            set(gca, 'color', ExpColors(exp,:));
            AxIndx = AxIndx+1;
            grid on
            %axis off
        end
    end
end
set(Axh, 'ylim', Ylims);
suptitle(sprintf('%s %s channel %d cell %d', Subject, Date, Channel, CellIndx));

twin = inputdlg({'Time window start (ms)','Time window end (ms)'},'Time window', 1, {sprintf('%0.2f', Twin(1)), sprintf('%0.2f', Twin(2))});
if isempty(twin)
    return
end
Twin = [str2num(twin{1}), str2num(twin{2})];
export_fig(fullfile(SaveDir, sprintf('AllExpressions_%s_%s_ch%03d_cell%d.png',Subject, Date, Channel, CellIndx)), '-png');


%% ================ PLOT EXPRESSION SUMMARY FIGURE =======================
mfh                     = figure('name',sprintf('%s %s - cell %d orientation analysis', Subject, Date, Cell),'position',[1,1,1223, 912]);
Colors                  = jet(size(StereoSDF,2));
for exp = 1:size(StereoSDF,1)
    Axh(exp)            = subplot(3,4,exp);
    for dist = 1:size(StereoSDF,2)
        AllEl{dist} = [];
        for d = 1:size(StereoSDF,3)
            AllEl{dist} = [AllEl{dist}; StereoSDF{exp,dist,d}];
        end
        phdist(dist) = plot(HistBins, mean(AllEl{dist}),'linewidth',2,'color',Colors(dist,:));
        hold on;
    end
    StimH = plot([0, 300], [1, 1], '-k', 'linewidth', 10);
    YlimsC = get(gca,'ylim');
    ph = patch([Twin([1,1]),Twin([2,2])], [0,YlimsC([2,2]),0], 0, 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
    uistack(ph, 'bottom');
    set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(StimH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    if exp == size(StereoSDF,1)
        for dist = 1:numel(Params.Distances)
            LegendTextDist{dist}  = sprintf('%d cm',Params.Distances(dist));
        end
        legend(phdist, LegendTextDist, 'location', 'northeast', 'fontsize',14);
    end

    xlabel('Time (ms)', 'fontsize', 18);
    ylabel('Firing rate (Hz)',  'fontsize', 18);
    set(gca,'fontsize',14,'xlim',[HistBins(1), HistBins(end)],'tickdir','out');
    grid on
    box off
 	Params.Expressions{exp}(1) = upper(Params.Expressions{exp}(1));
    title(Params.Expressions{exp},'fontsize', 18)
end

%============== Plot depth profile x expression
for exp = 1:size(StereoSDF,1)
    Axh(4+exp)            = subplot(3,4,exp+4);
    for d = 1:size(StereoSDF,3)
        AllEl{d} = [];
        for dist = 1:size(StereoSDF,2)
            AllEl{d} = [AllEl{d}; StereoSDF{exp,dist,d}];
        end
        phdepth(d) = plot(HistBins, mean(AllEl{d}),'linewidth',2,'color',Colors(d,:));
        hold on;
    end
    StimH = plot([0, 300], [1, 1], '-k', 'linewidth', 10);
    YlimsC = get(gca,'ylim');
    ph = patch([Twin([1,1]),Twin([2,2])], [0,YlimsC([2,2]),0], 0, 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
    uistack(ph, 'bottom');
    set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(StimH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    if exp == size(StereoSDF,1)
        legend(phdepth, Params.Depths, 'location', 'northeast', 'fontsize',14);
    end
    xlabel('Time (ms)', 'fontsize', 18);
    if exp == 1
        ylabel('Firing rate (Hz)',  'fontsize', 18);
    end
    set(gca,'fontsize',14,'xlim',[HistBins(1), HistBins(end)],'tickdir','out');
    grid on
    box off
 	Params.Expressions{exp}(1) = upper(Params.Expressions{exp}(1));
    title(Params.Expressions{exp},'fontsize', 18)
end
linkaxes(Axh(1:8))

%============== Plot expression x orientation firing rate matrices
for exp = 1:size(OrientationSDF,1)
    Axh(8+exp)             	= subplot(3,4,8+exp);
    PixelOffset             = round([HeadImSize(1)/size(OrientDiffMat{exp},2), HeadImSize(2)/size(OrientDiffMat{exp},1)]/2);
    imh(exp,1)            	= imagesc([PixelOffset(2), HeadImSize(2)-PixelOffset(2)],[PixelOffset(1), HeadImSize(1)-PixelOffset(1)], OrientDiffMat{exp}');                                                       
    hold on
    imh(exp,2)            	= image(HeadImages{exp});                                                       % Draw head image overlay
    alpha(imh(exp,2), HeadImAlpha{exp});                                                                    % Set alpha transparency
    axis equal tight
    colormap hot
    box off
    set(gca,'xtick',linspace(PixelOffset(2), HeadImSize(2)-PixelOffset(2), size(OrientDiffMat{exp},1)),...
            'xticklabel',Params.Azimuths, ...
            'ytick', linspace(PixelOffset(1), HeadImSize(1)-PixelOffset(1), size(OrientDiffMat{exp},2)),...
            'yticklabel',Params.Elevations,...
            'fontsize',16,...
            'tickdir', 'out');
 	
   	xlabel('Azimuth (°)', 'fontsize', 18);
	if exp==1  
        ylabel('Elevation (°)',  'fontsize', 18);
    end
    AxColLims(exp,:) = get(gca,'clim');
end
cbh = colorbar;                                                             % Add a color bar
set(cbh.Label, 'String', '\Delta Firing rate (Hz)', 'FontSize', 16);        % Give the color bar a title
cbh.Position = cbh.Position+[0.05, 0,0,0];                                  % Adjust colorbar position
Ax3Pos = get(Axh(11),'position');                                           % Get position of penultimate axis
Ax4Pos = get(Axh(12),'position');                                           % Get position of last axis
set(Axh(12),'position', [Ax4Pos([1,2]), Ax3Pos(3), Ax4Pos(4)]);             % Adjust width of last axis
set(Axh(9:12),'clim', [min(AxColLims(:,1)), max(AxColLims(:,2))]);          % Link colormap scales for all axes


% %============== Plot orientation tuning curves for azimuth angle
% Colors    	= jet(size(OrientationSDF,2)+1);
% for exp = 1:numel(OrientDiffMat)
%     Axh(8+exp) 	= subplot(3,4,8+exp);
%     for el = 1:size(OrientDiffMat{exp},2)
%     %     [ha, hb, hc] = shadedplot(1:size(OrientDiffMat{exp},1), [OrientDiffMat{exp}(:,el)-OrientRawSEM{exp}(:, el)]', [OrientDiffMat{exp}(:,el)+OrientRawSEM{exp}(:, el)]', [0.75, 0.75, 0.75]);
%     %     delete([hb, hc]);
%     %     hold on
%         plh(el) = plot(OrientDiffMat{exp}(:,el),'linewidth',2);
%         hold on;
%         ebh(el) = errorbar(1:size(OrientDiffMat{exp}, 1), OrientDiffMat{exp}(:,el), OrientRawSEM{exp}(:, el), OrientRawSEM{exp}(:, el));
%         set(ebh(el), 'color', get(plh(el), 'color'));
%     end
%     plh(el+1) = plot(mean(OrientDiffMat{exp}'),'--k','linewidth',3);
%     grid on
%     box off
% 
%     LegendTextEl{end+1} = 'Mean';
%     legend(plh, LegendTextEl, 'fontsize',18);
%     set(gca,'xlim',[0.5, 7.5],'xtick',1:1:7,'xticklabel', Params.Azimuths,'tickdir','out','fontsize',16);
%     xlabel('Azimuth (°)', 'fontsize', 16);
%     if exp == 1
%         ylabel('\Delta Firing rate (Hz)',  'fontsize', 18);
%     end
% end

suptitle(sprintf('Facial expression summary - %s %s channel %d cell %d', Subject, Date, Channel, CellIndx));

%============== Save figure
saveas(mfh, fullfile(Append, '/PROCDATA/murphya/Physio/StereoFaces/FacialExpressionTuning', Subject,sprintf('ExpressionTuning_%s_%s_ch%03d_cell%d.eps',Subject, Date, Channel, CellIndx)),'epsc');
%  export_fig(fullfile(Append, '/projects/murphya/MacaqueFace3D/PilotData/PNG/',sprintf('OrientationMatrix_%s_%s_cell%d.png',Subject, Date, Cell)), '-png');

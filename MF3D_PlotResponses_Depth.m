function MF3D_PlotResponses_Depth(Subject, Date, Channel, CellIndx, Output)

%====================== MF3D_PlotResponses_Depth.m ========================
% Plot facial epxression tuning for STS cells tested in 'StereoFace' pilot
% experiments.
%
%==========================================================================

if nargin == 0
    Subject     = 'Matcha';
    Date        = '20160616';
    Channel     = 5;
    CellIndx   	= 1;
    Output      = 'gif';
end

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
    HeadExpImage = fullfile('/Volumes/projects/murphya/MacaqueFace3D/PilotData/PNGs/DepthStim.png');
    [HeadIm,cm,HeadAlpha]   = imread(HeadExpImage);                                   % Load head depth overlay
    HeadImSize              = size(HeadIm);
    HeadImages{exp}         = HeadIm;
    HeadImAlpha{exp}        = HeadAlpha;
end

%============= Settings
Cell                = find(ismember(ChIndx(:,[2,3]), [Channel, CellIndx], 'rows'));
Params.Depths       = {'Concave','Flat','Convex'};
Params.DepthIndx    = [-1, 0, 1];
RenameExp   = {'fear','Bared teeth'; 'lipsmack', 'Pursed lips';'threat', 'Open mouth'};
for n = 1:size(RenameExp,1)
    Params.Expressions{~cellfun(@isempty, strfind(Params.Expressions, RenameExp{n,1}))} = RenameExp{n,2};
end
if ~isempty(find(Params.ConditionMatrix(:,6)==-1))
    for p = numel(Params.DepthIndx):-1:1
        Params.ConditionMatrix(Params.ConditionMatrix(:,6)==Params.DepthIndx(p), 6) = p;
    end
end
        

%========== Modified plot
fh          = figure('position', get(0,'screensize'));              % Open fullscreen figure window
Axh         = tight_subplot(6,12,0, 0.05,0.05);                     % Generate 72 axes
RastAxIndx  = [1,2,3,25,26,27,49,50,51];                            % Specify which axes are for raster plots
SDFAxIndx   = RastAxIndx+12;                                        % Specify axes for SDF plots
Ylims       = [0 100];                                            	% Specify y-axis limits for SDFs (spikes per second)
Xlims       = [-100, 400];                                          % Specify x-axis limits for all axes (ms)          
WinColor    = [0.5, 0.5, 0.5];
WinAlpha    = 0.5;
AxesXshift  = [-0.01, 0.01, 0.03,0.05];
ExpColors   = [0.75,0.75,1; 1,1,0.75; 0.75,1,0.75; 1, 0.75, 0.75];
SDFdims     = [0.07,0.08];
Rastdims    = [0.07,0.05];
PlotXpos    = 0.04+[0,1,2]*SDFdims(1);
PlotXpos    = [PlotXpos, PlotXpos+0.23, PlotXpos+0.46, PlotXpos+0.69];
PlotYpos    = 0.02+cumsum([0, SDFdims(2), Rastdims(2),SDFdims(2), Rastdims(2),SDFdims(2)]);
PlotYpos    = PlotYpos(end:-1:1)+0.5;
ph          = [];
lh          = [];
AllSDFAxIndx = [];
Twin        = [80, 150];                                            % Specify time window to calculate mean response from (ms) *if output type is 'svg'
BaseWin     = [-50 50];                                         	% Specify time window to calculate baseline response from (ms)
BaseIndx    = find(HistBins>BaseWin(1) & HistBins<BaseWin(2));      % Calculate which bins to include as baseline measure
TwinPos 	= 0:10:300;                                             % Specify timepoints of window start position (ms)
TwinWidth   = 50;                                                   % Specify window width (ms)

%========== Loop through condiions
for exp = 1:numel(Params.Expressions)
    AxIndx = 1;
    if exp > 1
        RastAxIndx  = RastAxIndx+3;
        SDFAxIndx   = SDFAxIndx+3;
    end
    AllSDFAxIndx = [AllSDFAxIndx,SDFAxIndx];
    
    
    %========= Plot depth data
    for dist = 1:numel(Params.Distances)
        for d = 1:numel(Params.Depths)
            CondIndx = find(ismember(Params.ConditionMatrix(:,[1,4,6]), [exp,dist,d],'rows'));
            StereoSDF{exp,d,dist} = [];
            
            %========= Plot depth rasters
            axes(Axh(RastAxIndx(AxIndx)));
            line = 1;
            for c = 1:numel(CondIndx)
                for t = 1:size(AllSpikes, 3)
                     if ~isnan(AllSpikes{Cell, CondIndx(c), t})
                        StereoSDF{exp,d,dist}(end+1,:) = hist(AllSpikes{Cell, CondIndx(c), t}, HistBins)*10^3/diff(HistBins([1,2]));
                        for sp = 1:numel(AllSpikes{Cell, CondIndx(c), t})                                                          % For each spike...
                            rph = plot(repmat(AllSpikes{Cell, CondIndx(c), t}(sp), [1,2]), [line-1, line], '-k');                   % Draw a vertical line
                         	hold on;
                     	end
                        line = line+1;
                     end
                end
            end
            axis tight
            box off
            set(gca,'yticklabels',[], 'xticklabels',[]);
            nx = rem(RastAxIndx(AxIndx),12);
            ny = ceil(RastAxIndx(AxIndx)/12);
            if nx == 0
                nx = 12;
            end
            set(gca, 'Position', [PlotXpos(nx), PlotYpos(ny),Rastdims]);
            if AxIndx == 2
                Params.Expressions{exp}(1) = upper(Params.Expressions{exp}(1));
                title(Params.Expressions{exp}, 'fontsize', 18);
            end
            
            %========= Calculate time window matrices
            for t = 1:numel(TwinPos)
                Twin                            = [TwinPos(t)-TwinWidth/2, TwinPos(t)+TwinWidth/2];                                                    % Specify time window to calculate mean response from (ms)
                Tindx                           = find(HistBins>Twin(1) & HistBins<Twin(2));     
                StereoRawMean(exp,d,dist,t)     = mean(mean(StereoSDF{exp,d,dist}(:,Tindx)));
                StereoBaseMean(exp,d,dist,t)    = mean(mean(StereoSDF{exp,d,dist}(:,BaseIndx)));
                StereoRawSEM{exp,t}(d,dist) 	= std(std(StereoSDF{exp,d,dist}(:,Tindx)))/sqrt(size(StereoSDF{exp,d,dist},1));
                StereoDiffMat{exp,t}(d,dist)    = StereoRawMean(exp,d,dist,t)-StereoBaseMean(exp,d,dist,t);
              	MaxDiffMat(exp,t)               = max(StereoDiffMat{exp,t}(:));
                MinDiffMat(exp,t)               = min(StereoDiffMat{exp,t}(:));
            end

            %========= Plot depth data SDFs
            axes(Axh(SDFAxIndx(AxIndx)));
            BinSEM{exp,d,dist} = std(StereoSDF{exp,d,dist})/sqrt(size(StereoSDF{exp,d,dist}, 1));
            [ha, hb, hc] = shadedplot(HistBins, mean(StereoSDF{exp,d,dist})-BinSEM{exp,d,dist}, mean(StereoSDF{exp,d,dist})+BinSEM{exp,d,dist}, [1,0.5,0.5]);
            hold on;
            delete([hb, hc]);
            plot(HistBins, mean(StereoSDF{exp,d,dist}), '-r');
            ph(end+1) = patch(Twin([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
            lh(end+1) = plot([0,0],Ylims, '-k','linewidth',2);
            uistack(ph(end), 'bottom')
          	set(gca,'xlim', Xlims, 'ylim', Ylims);
            if ismember(SDFAxIndx(AxIndx), [1:12:(6*12)])
                ylabel(sprintf('Dist = %d cm', Params.Distances(dist)), 'fontsize', 16);
            else
                set(gca,'yticklabels',[]);
            end
            if SDFAxIndx(AxIndx) >= (5*12)+1
                xlabel(sprintf('%s', Params.Depths{d}), 'fontsize', 16);
            else
                set(gca,'xticklabels',[]);
            end
%             set(gca, 'color', ExpColors(exp,:));
          	nx = rem(SDFAxIndx(AxIndx),12);
            ny = ceil(SDFAxIndx(AxIndx)/12);
         	if nx == 0
                nx = 12;
            end
            set(gca, 'Position', [PlotXpos(nx), PlotYpos(ny),SDFdims]);
            AxIndx = AxIndx+1;
            grid on
            box off
            %axis off
            drawnow
        end
    end
    
end

%========= Adjust y-axis limits
Ylims = [0, ceil(max(MaxDiffMat(:))/50)*50];
set(Axh(AllSDFAxIndx),'ylim', Ylims);
set(ph, 'Ydata', Ylims([1,2,2,1]));
set(lh, 'Ydata', Ylims);
    
%============== Plot expression x orientation firing rate matrices
t = 1;
MatClims = [min(MinDiffMat(:)), max(MaxDiffMat(:))];
for exp = 1:numel(Params.Expressions)
    AxMat(exp)           	= axes('position', [PlotXpos(1+(exp-1)*3)-0.04, 0.1, 0.32 0.32]);
    PixelOffset             = round([HeadImSize(1)/size(StereoDiffMat{exp, t},2), HeadImSize(2)/size(StereoDiffMat{exp, t},1)]/2);
    imh(exp,1)            	= imagesc([PixelOffset(2), HeadImSize(2)-PixelOffset(2)],[PixelOffset(1), HeadImSize(1)-PixelOffset(1)], StereoDiffMat{exp, t}');                                                       
    hold on
    imh(exp,2)            	= image(HeadImages{exp});                                                       % Draw head image overlay
    alpha(imh(exp,2), HeadImAlpha{exp});                                                                    % Set alpha transparency
    axis equal tight
    colormap hot
    box off
    set(gca,'xtick',linspace(PixelOffset(2), HeadImSize(2)-PixelOffset(2), size(StereoDiffMat{exp, t},1)),...
            'xticklabel',Params.Depths, ...
            'ytick', linspace(PixelOffset(1), HeadImSize(1)-PixelOffset(1), size(StereoDiffMat{exp, t},2)),...
            'yticklabel',Params.Distances,...
            'fontsize',16,...
            'tickdir', 'out');

    xlabel('Depth profile', 'fontsize', 18);
    if exp==1  
        ylabel('Position in depth (cm)',  'fontsize', 18);
    end
    AxColLims(exp,:) = get(gca,'clim');
end
    
cbh             = colorbar;                                               	% Add a color bar
set(cbh.Label, 'String', '\Delta Firing rate (Hz)', 'FontSize', 16);        % Give the color bar a title
cbh.Position    = cbh.Position+[0.03, 0,0,0];                            	% Adjust colorbar position
Ax3Pos          = get(AxMat(3),'position');                                	% Get position of penultimate axis
Ax4Pos          = get(AxMat(4),'position');                               	% Get position of last axis
set(AxMat(4),'position', [Ax4Pos([1,2]), Ax3Pos(3), Ax4Pos(4)]);          	% Adjust width of last axis
set(AxMat,'clim', [min(AxColLims(:,1)), max(AxColLims(:,2))]);              % Link colormap scales for all axes
suptitle(sprintf('%s %s channel %d cell %d', Subject, Date, Channel, CellIndx), 20);


%============= Animate moving window
if strcmpi(Output, 'gif')
    giffilename = fullfile(SaveDir,sprintf('StereoTimeline_%s_%s_ch%03d_cell%d.gif',Subject, Date, Channel, CellIndx));
    set(AxMat,'clim', MatClims);                                                % Adjust colormap scales for all axes at all timepoints
    ClockH = axes('position', [0.05,0.95,0.1,0.05]);
    ClockTextH = text( 0.5, 0, '0 ms', 'FontSize', 30, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom' ) ;
    axis off
    for t = 1:numel(TwinPos)
        Twin = [TwinPos(t)-TwinWidth/2, TwinPos(t)+TwinWidth/2];    
        set(ClockTextH, 'string', sprintf('%d ms',TwinPos(t)));
        set(ph, 'Xdata', Twin([1,1,2,2]));
        set(lh, 'Xdata', repmat(TwinPos(t),[1,2]));
        for exp = 1:numel(Params.Expressions)
            set(imh(exp,1), 'cdata', StereoDiffMat{exp, t}');
        end
        drawnow
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if t == 1
          imwrite(imind,cm,giffilename,'gif', 'Loopcount',inf,'DelayTime', 0.3);
        else
          imwrite(imind,cm,giffilename,'gif','WriteMode','append');
        end
    end
elseif strcmpi(Output, 'svg')
    print('-painters', '-dsvg', fullfile(SaveDir, sprintf('StereoTuning_%s_%s_ch%03d_cell%d.svg',Subject, Date, Channel, CellIndx)))
    export_fig(fullfile(SaveDir, sprintf('StereoTuning_%s_%s_ch%03d_cell%d.png',Subject, Date, Channel, CellIndx)), '-png');
end


%% ============= Pool neural data across scales and depths
% FigVar = {'Orient','Depth'};
% for f = 1:2
%     fh{f}      = figure('position', [-1911, 168, 1799, 840]);
%     Axh{f}     = tight_subplot(6,12,0, 0.05,0.05);
% end

% 
%     
%     %========= Plot depth data
%     for dist = 1:numel(Params.Distances)
%         for d = 1:numel(Params.Depths)
%             CondIndx = find(ismember(Params.ConditionMatrix(:,[1,4,6]), [exp,d,dist],'rows'));
%             StereoSDF{exp,d,dist} = [];
%             for c = 1:numel(CondIndx)
%                 for t = 1:size(AllSpikes, 3)
%                      if ~isnan(AllSpikes{Cell, CondIndx(c), t})
%                          StereoSDF{exp,d,dist}(end+1,:) = hist(AllSpikes{Cell, CondIndx(c), t}, HistBins)*10^3/diff(HistBins([1,2]));
%                      end
%                 end
%             end
%         	StereoRawMean(exp,dist,d)   = mean(mean(StereoSDF{exp,dist,d}(:,Tindx)));
%             StereoBaseMean(exp,dist,d)  = mean(mean(StereoSDF{exp,dist,d}(:,BaseIndx)));
%             StereoRawSEM{exp}(dist,d) 	= std(std(StereoSDF{exp,dist,d}(:,Tindx)))/sqrt(size(StereoSDF{exp,dist,d},1));
%             StereoDiffMat{exp}(dist,d)  = StereoRawMean(exp,dist,d)-StereoBaseMean(exp,dist,d);
%             
%             %========= Plot depth data
%             axes(Axh{2}(ExpAxIndx(AxIndx(2))));
%             BinSEM{exp,dist,d} = std(StereoSDF{exp,dist,d})/sqrt(size(StereoSDF{exp,dist,d}, 1));
%             [ha, hb, hc] = shadedplot(HistBins, mean(StereoSDF{exp,dist,d})-BinSEM{exp,dist,d}, mean(StereoSDF{exp,dist,d})+BinSEM{exp,dist,d}, [1,0.5,0.5]);
%             hold on;
%             delete([hb, hc]);
%             plot(HistBins, mean(StereoSDF{exp,dist,d}), '-r');
%             ph = patch(Twin([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
%             uistack(ph, 'bottom')
%             
%             if ismember(ExpAxIndx(AxIndx(2)), [1:6:(6*6)])
%                 ylabel(sprintf('Dist = %d cm', Params.Distances(dist)), 'fontsize', 16);
%             end
%             if ExpAxIndx(AxIndx(2)) >= (5*6)+1
%                 xlabel(sprintf('%s', Params.Depths{d}), 'fontsize', 16);
%             end
%             set(gca, 'color', ExpColors(exp,:));
%             AxIndx(2) = AxIndx(2)+1;
%             grid on
%             %axis off
%         end
%     end
% end
% for f = 1:2
%     set(Axh{f}, 'ylim', Ylims);
%     figure(fh{f})
%     suptitle(sprintf('%s %s channel %d cell %d', Subject, Date, Channel, CellIndx));
% end
%     
% twin = inputdlg({'Time window start (ms)','Time window end (ms)'},'Time window', 1, {sprintf('%0.2f', Twin(1)), sprintf('%0.2f', Twin(2))});
% if isempty(twin)
% %     for f = 1:numel(fh)
% %         close(fh{f});
% %     end
%     return
% end
% Twin = [str2num(twin{1}), str2num(twin{2})];
% for f = 1:2
%     figure(fh{f});
%     export_fig(fullfile(SaveDir, sprintf('Exp_x_%s_%s_%s_ch%03d_cell%d.png',FigVar{f}, Subject, Date, Channel, CellIndx)), '-png');
% end
% 
% %% ================ PLOT EXPRESSION SUMMARY FIGURE =======================
% mfh                     = figure('name',sprintf('%s %s - cell %d orientation analysis', Subject, Date, Cell),'position',[1,1,1223, 912]);
% Colors                  = jet(size(StereoSDF,2));
% for exp = 1:size(StereoSDF,1)
%     Axhd(exp)            = subplot(3,4,exp);
%     for dist = 1:size(StereoSDF,2)
%         AllEl{dist} = [];
%         for d = 1:size(StereoSDF,3)
%             AllEl{dist} = [AllEl{dist}; StereoSDF{exp,dist,d}];
%         end
%         phdist(dist) = plot(HistBins, mean(AllEl{dist}),'linewidth',2,'color',Colors(dist,:));
%         hold on;
%     end
%     StimH = plot([0, 300], [1, 1], '-k', 'linewidth', 10);
%     YlimsC = get(gca,'ylim');
%     ph = patch([Twin([1,1]),Twin([2,2])], [0,YlimsC([2,2]),0], 0, 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
%     uistack(ph, 'bottom');
%     set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(StimH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
%     if exp == size(StereoSDF,1)
%         for dist = 1:numel(Params.Distances)
%             LegendTextDist{dist}  = sprintf('%d cm',Params.Distances(dist));
%         end
%         legend(phdist, LegendTextDist, 'location', 'northeast', 'fontsize',14);
%     end
% 
%     xlabel('Time (ms)', 'fontsize', 18);
%     ylabel('Firing rate (Hz)',  'fontsize', 18);
%     set(gca,'fontsize',14,'xlim',[HistBins(1), HistBins(end)],'tickdir','out');
%     grid on
%     box off
%  	Params.Expressions{exp}(1) = upper(Params.Expressions{exp}(1));
%     title(Params.Expressions{exp},'fontsize', 18)
% end
% 
% %============== Plot depth profile x expression
% for exp = 1:size(StereoSDF,1)
%     Axhd(4+exp)            = subplot(3,4,exp+4);
%     for d = 1:size(StereoSDF,3)
%         AllEl{d} = [];
%         for dist = 1:size(StereoSDF,2)
%             AllEl{d} = [AllEl{d}; StereoSDF{exp,dist,d}];
%         end
%         phdepth(d) = plot(HistBins, mean(AllEl{d}),'linewidth',2,'color',Colors(d,:));
%         hold on;
%     end
%     StimH = plot([0, 300], [1, 1], '-k', 'linewidth', 10);
%     YlimsC = get(gca,'ylim');
%     ph = patch([Twin([1,1]),Twin([2,2])], [0,YlimsC([2,2]),0], 0, 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
%     uistack(ph, 'bottom');
%     set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(StimH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     if exp == size(StereoSDF,1)
%         legend(phdepth, Params.Depths, 'location', 'northeast', 'fontsize',14);
%     end
%     xlabel('Time (ms)', 'fontsize', 18);
%     if exp == 1
%         ylabel('Firing rate (Hz)',  'fontsize', 18);
%     end
%     set(gca,'fontsize',14,'xlim',[HistBins(1), HistBins(end)],'tickdir','out');
%     grid on
%     box off
%  	Params.Expressions{exp}(1) = upper(Params.Expressions{exp}(1));
%     title(Params.Expressions{exp},'fontsize', 18)
% end
% linkaxes(Axhd(1:8))
% 
% %============== Plot expression x orientation firing rate matrices
% for exp = 1:size(OrientationSDF,1)
%     Axhd(8+exp)             	= subplot(3,4,8+exp);
%     PixelOffset             = round([HeadImSize(1)/size(OrientDiffMat{exp},2), HeadImSize(2)/size(OrientDiffMat{exp},1)]/2);
%     imh(exp,1)            	= imagesc([PixelOffset(2), HeadImSize(2)-PixelOffset(2)],[PixelOffset(1), HeadImSize(1)-PixelOffset(1)], OrientDiffMat{exp}');                                                       
%     hold on
%     imh(exp,2)            	= image(HeadImages{exp});                                                       % Draw head image overlay
%     alpha(imh(exp,2), HeadImAlpha{exp});                                                                    % Set alpha transparency
%     axis equal tight
%     colormap hot
%     box off
%     set(gca,'xtick',linspace(PixelOffset(2), HeadImSize(2)-PixelOffset(2), size(OrientDiffMat{exp},1)),...
%             'xticklabel',Params.Azimuths, ...
%             'ytick', linspace(PixelOffset(1), HeadImSize(1)-PixelOffset(1), size(OrientDiffMat{exp},2)),...
%             'yticklabel',Params.Elevations,...
%             'fontsize',16,...
%             'tickdir', 'out');
%  	
%    	xlabel('Azimuth (°)', 'fontsize', 18);
% 	if exp==1  
%         ylabel('Elevation (°)',  'fontsize', 18);
%     end
%     AxColLims(exp,:) = get(gca,'clim');
% end
% cbh = colorbar;                                                             % Add a color bar
% set(cbh.Label, 'String', '\Delta Firing rate (Hz)', 'FontSize', 16);        % Give the color bar a title
% cbh.Position = cbh.Position+[0.05, 0,0,0];                                  % Adjust colorbar position
% Ax3Pos = get(Axhd(11),'position');                                           % Get position of penultimate axis
% Ax4Pos = get(Axhd(12),'position');                                           % Get position of last axis
% set(Axhd(12),'position', [Ax4Pos([1,2]), Ax3Pos(3), Ax4Pos(4)]);             % Adjust width of last axis
% set(Axhd(9:12),'clim', [min(AxColLims(:,1)), max(AxColLims(:,2))]);          % Link colormap scales for all axes
% 
% 
% % %============== Plot orientation tuning curves for azimuth angle
% % Colors    	= jet(size(OrientationSDF,2)+1);
% % for exp = 1:numel(OrientDiffMat)
% %     Axhd(8+exp) 	= subplot(3,4,8+exp);
% %     for el = 1:size(OrientDiffMat{exp},2)
% %     %     [ha, hb, hc] = shadedplot(1:size(OrientDiffMat{exp},1), [OrientDiffMat{exp}(:,el)-OrientRawSEM{exp}(:, el)]', [OrientDiffMat{exp}(:,el)+OrientRawSEM{exp}(:, el)]', [0.75, 0.75, 0.75]);
% %     %     delete([hb, hc]);
% %     %     hold on
% %         plh(el) = plot(OrientDiffMat{exp}(:,el),'linewidth',2);
% %         hold on;
% %         ebh(el) = errorbar(1:size(OrientDiffMat{exp}, 1), OrientDiffMat{exp}(:,el), OrientRawSEM{exp}(:, el), OrientRawSEM{exp}(:, el));
% %         set(ebh(el), 'color', get(plh(el), 'color'));
% %     end
% %     plh(el+1) = plot(mean(OrientDiffMat{exp}'),'--k','linewidth',3);
% %     grid on
% %     box off
% % 
% %     LegendTextEl{end+1} = 'Mean';
% %     legend(plh, LegendTextEl, 'fontsize',18);
% %     set(gca,'xlim',[0.5, 7.5],'xtick',1:1:7,'xticklabel', Params.Azimuths,'tickdir','out','fontsize',16);
% %     xlabel('Azimuth (°)', 'fontsize', 16);
% %     if exp == 1
% %         ylabel('\Delta Firing rate (Hz)',  'fontsize', 18);
% %     end
% % end
% 
% suptitle(sprintf('Facial expression summary - %s %s channel %d cell %d', Subject, Date, Channel, CellIndx));
% 
% %============== Save figure
% print('-painters', '-dsvg', fullfile(Append, '/PROCDATA/murphya/Physio/StereoFaces/FacialExpressionTuning', Subject,sprintf('ExpressionTuning_%s_%s_ch%03d_cell%d.svg',Subject, Date, Channel, CellIndx)))
% % saveas(mfh, fullfile(Append, '/PROCDATA/murphya/Physio/StereoFaces/FacialExpressionTuning', Subject,sprintf('ExpressionTuning_%s_%s_ch%03d_cell%d.eps',Subject, Date, Channel, CellIndx)),'epsc');
% export_fig(fullfile(Append, '/PROCDATA/murphya/Physio/StereoFaces/FacialExpressionTuning', Subject,sprintf('ExpressionTuning_%s_%s_ch%03d_cell%d.png',Subject, Date, Channel, CellIndx)), '-png');

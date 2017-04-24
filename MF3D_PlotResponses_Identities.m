function MF3D_PlotResponses_Identities(Subject, Date, Channel, CellIndx, Output)

%===================== MF3D_PlotResponses_Identities.m ====================
% Plot facial epxression tuning for STS cells tested in 'StereoFace' pilot
% experiments.
%
%==========================================================================

if nargin == 0
    Subject     = 'Matcha';
    Date        = '20160720';
    Channel     = 25;
    CellIndx   	= 1;
    Output      = 'gif';
end

switch Subject
    case 'Avalanche'
        if ~any(~cellfun(@isempty, strfind({'20160712','20160713'}, Date)))
            error('Invalid session for expression analysis!')
        end
    case 'Matcha'
        if ~any(~cellfun(@isempty, strfind({'20160719','20160720','20160721'}, Date)))
            error('Invalid session for expression analysis!')
        end
    case 'Spice'
%         if ~any(~cellfun(@isempty, strfind({'20160623','20160624'}, Date)))
            error('Invalid session for expression analysis!')
%         end
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
Cell 	= find(ismember(ChIndx(:,[2,3]), [Channel, CellIndx], 'rows'));
SaveDir = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/FacialIdentityTuning', Subject);
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

%=========== Load expression condition images
for id = 1:numel(Params.MonkeyIDs)
    HeadExpImage = fullfile('/Volumes/projects/murphya/MacaqueFace3D/PilotData/PNGs/',sprintf('MonkeyID_%d.png',id));
    [HeadIm,cm,HeadAlpha]   = imread(HeadExpImage);                                   % Load head orientation overlay
    HeadImSize{id}        	= size(HeadIm);
    HeadImages{id}          = HeadIm;
    HeadImAlpha{id}         = HeadAlpha;
end


%========== Plot settings
fh          = figure('position', get(0,'screensize'));              % Open fullscreen figure window
NoAxY       = numel(Params.Elevations)*2;
NoAxX       = numel(Params.Azimuths)*numel(Params.MonkeyIDs);       
Axh         = tight_subplot(NoAxY, NoAxX,0, 0.05,0.05);           	% Generate axes
RastAxIndx  = [];
for x = 1:3
    RastAxIndx(end+1:end+3) = (1:3) + (x-1)*(NoAxX*2);            	% Specify which axes are for raster plots
end
SDFAxIndx   = RastAxIndx+NoAxX;                                   	% Specify axes for SDF plots
Ylims       = [0 100];                                            	% Specify y-axis limits for SDFs (spikes per second)
Xlims       = [-100, 400];                                          % Specify x-axis limits for all axes (ms)       
WinColor    = [0.5, 0.5, 0.5];
WinAlpha    = 0.5;
AxesXshift  = [-0.01, 0.01, 0.03,0.05];
ExpColors   = [0.75,0.75,1; 1,1,0.75; 0.75,1,0.75; 1, 0.75, 0.75];
SDFdims     = [0.05,0.08];
Rastdims    = [0.05,0.05];
PlotXpos    = 0.04+[0,1,2]*SDFdims(1);
PlotXpos    = [PlotXpos, PlotXpos+0.19, PlotXpos+0.38, PlotXpos+0.57, PlotXpos+0.76];
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
for id = 1:numel(Params.MonkeyIDs)
    AxIndx = 1;
    if id > 1
        RastAxIndx  = RastAxIndx+3;
        SDFAxIndx   = SDFAxIndx+3;
    end
    AllSDFAxIndx = [AllSDFAxIndx,SDFAxIndx];
    for el = 1:numel(Params.Elevations)
        for az = 1:numel(Params.Azimuths)
            CondIndx = find(ismember(Params.ConditionMatrix(:,[1,2,3]), [id,az,el],'rows'));
            OrientationSDF{id,az,el} = [];
            
         	%========= Plot expression x head oreintation rasters
            axes(Axh(RastAxIndx(AxIndx)));
            line = 1;
            for c = 1:numel(CondIndx)
                for t = 1:size(AllSpikes, 3)
                     if ~isnan(AllSpikes{Cell, CondIndx(c), t})
                         OrientationSDF{id,az,el}(end+1,:) = hist(AllSpikes{Cell, CondIndx(c), t}, HistBins)*10^3/diff(HistBins([1,2]));
                         for sp = 1:numel(AllSpikes{Cell, CondIndx(c), t})                                                  % For each spike...
                            rph = plot(repmat(AllSpikes{Cell, CondIndx(c), t}(sp), [1,2]), [line-1, line], '-k');  	% Draw a vertical line
                            hold on;
                         end
                        line = line+1;
                     end
                end
            end
            axis tight
            box off
            set(gca,'yticklabels',[], 'xticklabels',[]);
            nx = rem(RastAxIndx(AxIndx),NoAxX);
            ny = ceil(RastAxIndx(AxIndx)/NoAxX);
            if nx == 0
                nx = NoAxX;
            end
            set(gca, 'Position', [PlotXpos(nx), PlotYpos(ny),Rastdims]);
            if AxIndx == 2
                title(sprintf('ID %d', id), 'fontsize', 18);
            end
            
            %========= Calculate time window matrices
            for t = 1:numel(TwinPos)
                Twin                            = [TwinPos(t)-TwinWidth/2, TwinPos(t)+TwinWidth/2];                                   	% Specify time window to calculate mean response from (ms)
                Tindx                           = find(HistBins>Twin(1) & HistBins<Twin(2));     
                OrientRawMean(id,az,el,t)      = mean(mean(OrientationSDF{id,az,el}(:,Tindx)));
                OrientBaseMean(id,az,el,t)     = mean(mean(OrientationSDF{id,az,el}(:,BaseIndx)));
                OrientRawSEM{id,t}(az,el)      = std(std(OrientationSDF{id,az,el}(:,Tindx)))/sqrt(size(OrientationSDF{id,az,el},1));
                OrientDiffMat{id,t}(az,el)     = OrientRawMean(id,az,el,t)-OrientBaseMean(id,az,el,t);
                MaxDiffMat(id,t)               = max(OrientDiffMat{id,t}(:));
                MinDiffMat(id,t)               = min(OrientDiffMat{id,t}(:));
            end

            %========= Plot expression x head oreintation SDF
            axes(Axh(SDFAxIndx(AxIndx)));
            BinSEM{id,az,el} = std(OrientationSDF{id,az,el})/sqrt(size(OrientationSDF{id,az,el}, 1));
            [ha, hb, hc] = shadedplot(HistBins, mean(OrientationSDF{id,az,el})-BinSEM{id,az,el}, mean(OrientationSDF{id,az,el})+BinSEM{id,az,el}, [1,0.5,0.5]);
            hold on;
            delete([hb, hc]);
            plot(HistBins, mean(OrientationSDF{id,az,el}), '-r');
            ph(end+1) = patch(Twin([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
            lh(end+1) = plot([0,0],Ylims, '-k','linewidth',2);
            uistack(ph(end), 'bottom');
            
            set(gca,'xlim', Xlims, 'ylim', Ylims,'xtick', Xlims(1):100:Xlims(2));
            if ismember(SDFAxIndx(AxIndx), [1:NoAxX:(6*NoAxX)])
                ylabel(sprintf('%d °', Params.Elevations(el)), 'fontsize', 16);
            else
                set(gca,'yticklabels',[]);
            end
            if SDFAxIndx(AxIndx) >= (5*NoAxX)+1
                xlabel(sprintf('%d °', Params.Azimuths(az)), 'fontsize', 16);
            else
                set(gca,'xticklabels',[]);
            end
%             set(gca, 'color', ExpColors(id,:));
          	nx = rem(SDFAxIndx(AxIndx),NoAxX);
            ny = ceil(SDFAxIndx(AxIndx)/NoAxX);
         	if nx == 0
                nx = NoAxX;
            end
            set(gca, 'Position', [PlotXpos(nx), PlotYpos(ny), SDFdims]);
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
for id = 1:numel(Params.MonkeyIDs)
 	AxMat(id)           	= axes('position', [PlotXpos(1+(id-1)*3)-0.08, 0.1, 0.32 0.32]);
    PixelOffset             = round([HeadImSize{id}(1)/size(OrientDiffMat{id, t},2), HeadImSize{id}(2)/size(OrientDiffMat{id, t},1)]/2);
    imh(id,1)            	= imagesc([PixelOffset(2), HeadImSize{id}(2)-PixelOffset(2)],[PixelOffset(1), HeadImSize{id}(1)-PixelOffset(1)], OrientDiffMat{id, t}');                                                       
    hold on
    imh(id,2)            	= image(HeadImages{id});                                                       % Draw head image overlay
    alpha(imh(id,2), HeadImAlpha{id});                                                                    % Set alpha transparency
    axis equal tight
    colormap hot
    box off
    set(gca,'xtick',linspace(PixelOffset(2), HeadImSize{id}(2)-PixelOffset(2), size(OrientDiffMat{id, t},1)),...
            'xticklabel',Params.Azimuths, ...
            'ytick', linspace(PixelOffset(1), HeadImSize{id}(1)-PixelOffset(1), size(OrientDiffMat{id, t},2)),...
            'yticklabel',Params.Elevations,...
            'fontsize',16,...
            'tickdir', 'out');

    xlabel('Azimuth (°)', 'fontsize', 18);
    if id==1  
        ylabel('Elevation (°)',  'fontsize', 18);
    end
    AxColLims(id,:) = get(gca,'clim');
end
    
cbh             = colorbar;                                               	% Add a color bar
set(cbh.Label, 'String', '\Delta Firing rate (Hz)', 'FontSize', 16);        % Give the color bar a title
cbh.Position    = cbh.Position+[0.02, 0,0,0];                            	% Adjust colorbar position
Ax3Pos          = get(AxMat(3),'position');                                	% Get position of penultimate axis
Ax4Pos          = get(AxMat(4),'position');                               	% Get position of last axis
set(AxMat(4),'position', [Ax4Pos([1,2]), Ax3Pos(3), Ax4Pos(4)]);          	% Adjust width of last axis
set(AxMat,'clim', [min(AxColLims(:,1)), max(AxColLims(:,2))]);              % Link colormap scales for all axes
suptitle(sprintf('%s %s channel %d cell %d', Subject, Date, Channel, CellIndx), 20);


%============= Animate moving window
if strcmpi(Output, 'gif')
    giffilename = fullfile(SaveDir, sprintf('IdentityTimeline_%s_%s_ch%03d_cell%d.gif',Subject, Date, Channel, CellIndx));
    set(AxMat,'clim', MatClims);                                                % Adjust colormap scales for all axes at all timepoints
  	ClockH = axes('position', [0.05,0.95,0.1,0.05]);
    ClockTextH = text( 0.5, 0, '0 ms', 'FontSize', 30, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom' ) ;
    axis off
    for t = 1:numel(TwinPos)
        Twin = [TwinPos(t)-TwinWidth/2, TwinPos(t)+TwinWidth/2];
        set(ClockTextH, 'string', sprintf('%d ms',TwinPos(t)));
        set(ph, 'Xdata', Twin([1,1,2,2]));
        set(lh, 'Xdata', repmat(TwinPos(t),[1,2]));
        for id = 1:numel(Params.MonkeyIDs)
            set(imh(id,1), 'cdata', OrientDiffMat{id, t}');
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
    print('-painters', '-dsvg', fullfile(SaveDir, sprintf('IdentityTuning_%s_%s_ch%03d_cell%d.svg',Subject, Date, Channel, CellIndx)))
    export_fig(fullfile(SaveDir, sprintf('IdentityTuning_%s_%s_ch%03d_cell%d.png',Subject, Date, Channel, CellIndx)), '-png');
end

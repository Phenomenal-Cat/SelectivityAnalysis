% function MF3D_PlotResponses_PCAIdentity(Subject, Date, Channel, CellIndx, Output)

%==================== MF3D_PlotResponses_PCAIdentity.m ====================
% Plot face-space object tuning for STS cells tested in 'StereoFace' pilot
% experiments.
%
%==========================================================================


Subject        = 'Spice';
SessionDate     = '170928';
ExpType         = 'GF';
Cell            = 73;
Channel         = 1;
headang       	= 5;
Output          = 'gif';
InterpOn        = 0;

%================= SET PATHS AND LOAD DATA
if ismac
    Prefix = '/Volumes'; 
else
    Prefix = [];
end
DataDir         = fullfile(Prefix,'/procdata/leathersml/physiology/raster/',[Subject, SessionDate],ExpType);
MatfileName    	= fullfile(DataDir,'StimAlignedSpikes.mat');
SaveDir         = fullfile(Prefix,'/procdata/murphya/Physio/StereoFaces/PCAFaceIdentityTuning/',Subject);
switch ExpType
    case 'GF'
        StimFile        = fullfile(Prefix, '/projects/leathersml/Stimuli/3Davatar/GrayFace/pilot170928/PCA_Grayscale_Processed/PCA_greyscale_processed.mat');
    case 'SF'
        StimFile        = fullfile(Prefix, '/projects/leathersml/Stimuli/3Davatar/SkinFace/pilot170929/PCA_Textured_Processed/PCA_Textured_processed.mat');
end

%============= Load image and data
HeadImFile  = fullfile(Prefix, '/projects/murphya/MorphBlender/StimArray_TexPCA.png');
HeadIm      = imread(HeadImFile);
[~,~,HeadImAlpha] = imread(HeadImFile);
HeadImSize  = size(HeadIm);

load(StimFile);
load(MatfileName);
AllSpikes = StimAlignedSpikes;
if ~exist(SaveDir,'dir')
    mkdir(SaveDir);
end

%================= SET PARAMETERS
PC1        	= unique([Design.PC1],'stable');
PC2       	= unique([Design.PC2],'stable');
HeadAngle  	= unique({Design.HeadAngle},'stable');
fh          = figure('position', get(0,'screensize')./[1 1 1 1]);   	% Open half-screen figure window
NoAxY       = numel(PC1)*2;
NoAxX       = numel(PC2);       
Axh         = tight_subplot(NoAxY, NoAxX,0, 0.05,0.05);           	% Generate axes
RastAxIndx  = [];
for r = 1:numel(PC1)
    RastAxIndx(end+1:end+NoAxX) = (1:NoAxX) + (r-1)*(NoAxX*2);    	% Specify which axes are for raster plots
end
SDFAxIndx   = RastAxIndx+NoAxX;                                   	% Specify axes for SDF plots
Ylims       = [0 100];                                            	% Specify y-axis limits for SDFs (spikes per second)
Xlims       = [-100, 400];                                          % Specify x-axis limits for all axes (ms)       
WinColor    = [0.5, 0.5, 0.5];                                      % Color of sliding window 
WinAlpha    = 0.5;                                                  % Transparency of sliding window
SDFdims     = [0.05,0.05];
Rastdims    = [0.05,0.04];
PlotXpos    = 0.05+(0:(NoAxX-1))*(SDFdims(1)+0.01);
PlotYpos    = 0.02+cumsum([0, SDFdims(2), repmat([Rastdims(2),SDFdims(2)],[1,(NoAxY/2)-1])]);
PlotYpos    = PlotYpos(end:-1:1)+0.20;
ph          = [];
lh          = [];
Twin        = [80, 150];                                            % Specify time window to calculate mean response from (ms) *if output type is 'svg'
BaseWin     = [-50 50];                                         	% Specify time window to calculate baseline response from (ms)
BaseIndx    = find(HistBins>BaseWin(1) & HistBins<BaseWin(2));      % Calculate which bins to include as baseline measure
TwinPos 	= 0:10:300;                                             % Specify timepoints of window start position (ms)
TwinWidth   = 50;                                                   % Specify window width (ms)

%===================

AxIndx  = 1;
for pc1 = 1:numel(PC1)
    for pc2 = 1:numel(PC2)
     	OrientationSDF{pc2,pc1} = [];
        
      	StimIndices = find(ismember([[Design.PC1]',[Design.PC2]'], [PC1(pc1),PC2(pc2)],'rows'));
        HAindx      = find(~cellfun(@isempty, strfind({Design(StimIndices).HeadAngle}, HeadAngle{headang})));
        CondIndx    = StimIndices(HAindx);
%         CondIndx    = StimIndices;
        
        %========= Plot rasters
    	axes(Axh(RastAxIndx(AxIndx)));
        line = 1;
        for c = 1:numel(CondIndx)
            %OrientationSDF{pc2,pc1}(end+1,:) = hist(StimAlignedSpikesAll{Cell, CondIndx(c)}, HistBins)*10^3/diff(HistBins([1,2]));
            OrientationSDF{pc2,pc1}(end+1,:) = PSTH{Cell}(CondIndx(c),:);
            for t = 1:size(AllSpikes, 3)
                 if ~isnan(AllSpikes{Cell, CondIndx(c), t})
                     for sp = 1:numel(AllSpikes{Cell, CondIndx(c), t})                                            	% For each spike...
                        rph = plot(repmat(AllSpikes{Cell, CondIndx(c), t}(sp)*1000, [1,2]), [line-1, line], '-k');       % Draw a vertical line
                        hold on;
                     end
                    line = line+1;
                 end
            end
        end
        plot([0,0],[0,line],'-k');
        axis tight
        set(gca,'yticklabels',[], 'xticklabels',[], 'xlim', Xlims);
        box off
        nx = rem(RastAxIndx(AxIndx),NoAxX);
        ny = ceil(RastAxIndx(AxIndx)/NoAxX);
        if nx == 0
            nx = NoAxX;
        end
        set(gca, 'Position', [PlotXpos(nx), PlotYpos(ny), Rastdims]);


        %========= Calculate all time window matrices
        for t = 1:numel(TwinPos)
            Twin                    	= [TwinPos(t)-TwinWidth/2, TwinPos(t)+TwinWidth/2];                                   	% Specify time window to calculate mean response from (ms)
            Tindx                   	= find(HistBins>Twin(1) & HistBins<Twin(2));    
            OrientRawMean(pc2, pc1, t)    = mean(mean(OrientationSDF{pc2,pc1}(:,Tindx)));
            OrientBaseMean(pc2, pc1, t)   = mean(mean(OrientationSDF{pc2,pc1}(:,BaseIndx)));
            OrientRawSEM{t}(pc2, pc1)     = std(std(OrientationSDF{pc2,pc1}(:,Tindx)))/sqrt(size(OrientationSDF{pc2,pc1},1));
            OrientDiffMat{t}(pc2, pc1)    = OrientRawMean(pc2, pc1, t)-OrientBaseMean(pc2, pc1, t);
            MaxDiffMat(t)               = max(OrientDiffMat{t}(:));
            MinDiffMat(t)               = min(OrientDiffMat{t}(:));
        end
        
        %========= Plot head oreintation SDFs
        [SDF{pc2,pc1},~] = msdf(OrientationSDF{pc2,pc1}','Gauss',10);
        axes(Axh(SDFAxIndx(AxIndx)));
        BinSEM{pc2, pc1} = std(OrientationSDF{pc2,pc1})/sqrt(size(OrientationSDF{pc2,pc1}, 1));
%         [ha, hb, hc] = shadedplot(HistBins, [SDF{pc2,pc1}-BinSEM{pc2, pc1}]', [SDF{pc2,pc1} +BinSEM{pc2, pc1}]', [1,0.5,0.5]);
        hold on;
%         delete([hb, hc]);
        plot(HistBins, SDF{pc2,pc1}', '-r');
        ph(end+1) = patch(Twin([1,1,2,2]), Ylims([1,2,2,1]), Ylims([1,2,2,1]), 'facecolor', WinColor, 'edgecolor', 'none', 'facealpha', WinAlpha);
        lh(end+1) = plot([0,0],Ylims, '-k','linewidth',2);
        uistack(ph(end), 'bottom')
        plot([0,0],ylim,'-k');
        
        set(gca,'xlim', Xlims, 'ylim', Ylims,'xtick', Xlims(1):100:Xlims(2));
        if ismember(SDFAxIndx(AxIndx), [1:NoAxX:(NoAxY*NoAxX)])
            ylabel(sprintf('%d', PC1(pc1)), 'fontsize', 16);
        else
            set(gca,'yticklabels',[]);
        end
        if SDFAxIndx(AxIndx) >= ((7-1)*NoAxX*2)+1
            xlabel(sprintf('%d', PC2(pc2)), 'fontsize', 16);
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

%========= Adjust y-axis limits
Ylims = [0, ceil(max(MaxDiffMat(:))/50)*50];
set(Axh(SDFAxIndx),'ylim', Ylims);
set(ph, 'Ydata', Ylims([1,2,2,1]));
set(lh, 'Ydata', Ylims);
    
%============== Plot orientation firing rate matrix
t = 1;
MatClims = [min(MinDiffMat(:)), max(MaxDiffMat(:))];

AxMat                   = axes('position', [0.55, 0.3, 0.45 0.5]);
PixelOffset             = round([HeadImSize(1)/size(OrientDiffMat{t},2), HeadImSize(2)/size(OrientDiffMat{t},1)]/2);
if InterpOn == 0
   
    imh(1)                  = imagesc([PixelOffset(2), HeadImSize(2)-PixelOffset(2)],[PixelOffset(1), HeadImSize(1)-PixelOffset(1)], OrientDiffMat{t}');                                                       
    hold on
    imh(2)                  = image(HeadIm);                                                        % Draw head image overlay
    alpha(imh(2), HeadImAlpha);   

elseif InterpOn == 1
    A       = OrientDiffMat{t};
    [X,Y]   = meshgrid(1:size(A,1));
    [Xi,Yi] = meshgrid(1:0.1:size(A,1));
    A       = interp2(X,Y,A,Xi,Yi,'linear');
    imh(1)                  = imagesc([PixelOffset(2), HeadImSize(2)-PixelOffset(2)],[PixelOffset(1), HeadImSize(1)-PixelOffset(1)], A');                                                       
    hold on
    imh(2)                  = image(HeadIm);                                                        % Draw head image overlay
    alpha(imh(2), HeadImAlpha);   
    axis off;
    
%     [X,Y]   = meshgrid(1:size(A,2), 1:size(A,1));
%     F       = scatteredInterpolant(X(:), Y(:), A(:), 'linear');
%     [U,V]   = meshgrid(linspace(1,size(A,2),50), linspace(1,size(A,1),50));
%    	imh(1)  = surf(U,V, F(U,V), 'edgecolor','none');
%     view([90,90]);
%     hold on;
%     patch([1,1,7,7],[1,7,7,1],[],'cdata', HeadIm);
    
end
                                                                      % Set alpha transparency
axis equal tight
colormap hot
box off
set(gca,'xtick',linspace(PixelOffset(2), HeadImSize(2)-PixelOffset(2), size(OrientDiffMat{t},1)),...
        'xticklabel',PC2, ...
        'ytick', linspace(PixelOffset(1), HeadImSize(1)-PixelOffset(1), size(OrientDiffMat{t},2)),...
        'yticklabel',PC1,...
        'fontsize',16,...
        'tickdir', 'out');
set(gca, 'xticklabel', []);                                                 % Turn off x-tick labels
% xlabel('PC 2 (SD)', 'fontsize', 18);
ylabel('PC 1 (SD)',  'fontsize', 18);
AxColLims       = get(gca,'clim');
cbh             = colorbar;                                               	% Add a color bar
set(cbh.Label, 'String', '\Delta Firing rate (Hz)', 'FontSize', 18);        % Give the color bar a title
cbh.Position    = cbh.Position+[0.03, 0,0,0];                            	% Adjust colorbar position
MatrixPos       = get(AxMat, 'position');

%============== Plot orientation tuning curves for azimuth angle
AxTuneX      = axes('position', [MatrixPos(1)+0.065, 0.07, MatrixPos(3)-0.1, 0.2]);
Colors    	= jet(size(OrientationSDF,2)+1);
if InterpOn == 0
    for pc1 = 1:size(OrientDiffMat{1},2)
    %     [ha, hb, hc] = shadedplot(1:size(OrientDiffMat,1), [OrientDiffMat(:,pc1)-OrientRawSEM(:, pc1)]', [OrientDiffMat(:,pc1)+OrientRawSEM(:, pc1)]', [0.75, 0.75, 0.75]);
    %     delete([hb, hc]);
    %     hold on
        plh{1}(pc1) = plot(OrientDiffMat{1}(:,pc1),'linewidth',2);
        hold on;
        ebh{1}(pc1) = errorbar(1:size(OrientDiffMat{1}, 1), OrientDiffMat{1}(:,pc1), OrientRawSEM{1}(:, pc1), OrientRawSEM{1}(:, pc1));
        set(ebh{1}(pc1), 'color', get(plh{1}(pc1), 'color'));
        LegendTextEl{pc1} = sprintf('%d', PC1(pc1));
        plh{1}(pc1+1) = plot(mean(OrientDiffMat{1}'),'--k','linewidth',3);
    end
else
% 	[ha, hb, hc] = shadedplot(linspace(1, numel(PC2), size(A,2)), [mean(A,2)-std(A)']', [mean(A,2)+std(A)']', [0.5,0.5,1]);
%  	hold on;
% 	set([hb, hc], 'visible','off');
%     plh{1}(2) = ha;
    plh{1}(1) = plot(linspace(1, numel(PC2), size(A,1)), mean(A,2),'linewidth',2);
end
grid on
box off
if InterpOn == 0
    LegendTextEl{end+1} = 'Mean';
    legend(plh{1}, LegendTextEl, 'location', 'EastOutside', 'fontsize',14);
end
set(gca,'xlim',[0.5, numel(PC2)+0.5],'xtick',1:1:numel(PC2),'xticklabel', PC2,'tickdir','out','fontsize',16);
set(gca,'position', [MatrixPos(1)+0.065, 0.07, MatrixPos(3)-0.1, 0.2]);
xlabel('PC 2 (SD)', 'fontsize', 18);
ylabel('\Delta Firing rate (Hz)',  'fontsize', 18);

%====== Y axis
AxTuneY      = axes('position', [MatrixPos(1)-0.035, MatrixPos(2), 0.1, MatrixPos(4)]);
if InterpOn == 0
    for pc2 = 1:size(OrientDiffMat{1},1)
    %     [ha, hb, hc] = shadedplot(1:size(OrientDiffMat,1), [OrientDiffMat(:,pc1)-OrientRawSEM(:, pc1)]', [OrientDiffMat(:,pc1)+OrientRawSEM(:, pc1)]', [0.75, 0.75, 0.75]);
    %     delete([hb, hc]);
    %     hold on
        plh{2}(pc2) = plot(OrientDiffMat{1}(pc2,:),'linewidth',2);
        hold on;
        ebh{2}(pc2) = errorbar(1:size(OrientDiffMat{1}, 1), OrientDiffMat{1}(pc2,:), OrientRawSEM{1}(pc2, :), OrientRawSEM{1}(pc2, :));
        set(ebh{2}(pc2), 'color', get(plh{2}(pc2), 'color'));
        LegendTextEl{pc2} = sprintf('%d', PC2(pc2));
    end
    plh{2}(pc2+1) = plot(mean(OrientDiffMat{1}),'--k','linewidth',3);
elseif InterpOn == 1
%     [ha, hb, hc] = shadedplot(linspace(1, numel(PC1), size(A,2)), mean(A,1)-std(A'), mean(A,1)+std(A'), [0.5,0.5,1]);
%  	hold on;
% 	delete([hb, hc]);
%     plh{2}(2) = ha;
    plh{2}(1) = plot(linspace(1, numel(PC1), size(A,2)), mean(A,1),'linewidth',2);
end
grid on
box off
view([90,90]);
set(gca,'xlim',[0.5, numel(PC1)+0.5],'xtick',1:1:numel(PC1),'xticklabel', PC1,'tickdir','out','fontsize',16);
xlabel('PC 1 (SD)', 'fontsize', 18);
% xlabel('\Delta Firing rate (Hz)',  'fontsize', 18);



suptitle(sprintf('%s %s cell %d', Subject, SessionDate, Cell), 20);


%============= Animate moving window
if strcmpi(Output, 'gif')
    giffilename = fullfile(SaveDir, sprintf('IdentityPCATimeline_%s_%s_ch%03d_cell%d.gif',Subject, SessionDate, Channel, Cell));
    if InterpOn == 1
        giffilename = fullfile(SaveDir, sprintf('IdentityPCATimeline_%s_%s_ch%03d_cell%d_interp.gif',Subject, SessionDate, Channel, Cell));
    end
    set(AxMat,'clim', MatClims);                                                % Adjust colormap scales for all axes at all timepoints
    set(AxTuneX,'ylim', MatClims);  
    set(AxTuneY,'ylim', MatClims); 
  	ClockH = axes('position', [0.1,0.9,0.1,0.05]);
    ClockTextH = text( 0.5, 0, '0 ms', 'FontSize', 24, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom' ) ;
    axis off
    for t = 1:numel(TwinPos)
        Twin = [TwinPos(t)-TwinWidth/2, TwinPos(t)+TwinWidth/2];
        set(ClockTextH, 'string', sprintf('%d ms',TwinPos(t)));
        set(ph, 'Xdata', Twin([1,1,2,2]));
        set(lh, 'Xdata', repmat(TwinPos(t),[1,2]));
        
        if InterpOn == 1
            A       = OrientDiffMat{t};
            [X,Y]   = meshgrid(1:size(A,1));
            [Xi,Yi] = meshgrid(1:0.1:size(A,1));
            A       = interp2(X,Y,A,Xi,Yi,'linear');
            set(imh(1), 'cdata', A');
            
%             A = OrientDiffMat{t};
%             [X,Y] = meshgrid(1:size(A,2), 1:size(A,1));
%             F = scatteredInterpolant(X(:), Y(:), A(:), 'linear');
%             [U,V] = meshgrid(linspace(1,size(A,2),50), linspace(1,size(A,1),50));
%             set(imh(1), 'cdata', F(U,V), 'zdata', F(U,V));
        else
            set(imh(1), 'cdata', OrientDiffMat{t}');
        end
        if InterpOn == 0
            for pc1 = 1:numel(PC1)
             	set(plh{1}(pc1), 'Ydata', OrientDiffMat{t}(:,pc1));
                if exist('ebh', 'var')
                    set(ebh{1}(pc1), 'Ydata', OrientDiffMat{t}(:,pc1), 'LData', OrientRawSEM{t}(:,pc1), 'Udata', OrientRawSEM{t}(:,pc1));
                end
            end
            for pc2 = 1:numel(PC2)
                set(plh{2}(pc2), 'Ydata', OrientDiffMat{t}(:,pc2));
                if exist('ebh', 'var')
                    set(ebh{2}(pc2), 'Ydata', OrientDiffMat{t}(:,pc2), 'LData', OrientRawSEM{t}(:,pc2), 'Udata', OrientRawSEM{t}(:,pc2));
                end
            end
                
        elseif InterpOn == 1
            set(plh{1}(1), 'Ydata', mean(A,1));
            set(plh{2}(1), 'Ydata', mean(A,2));
%             axes(AxTuneX)
%             delete(plh{1}(2));
%           	[ha, hb, hc] = shadedplot(linspace(1, numel(PC1), size(A,2)), mean(A,1)-std(A)', mean(A,1)+std(A)', [0.5,0.5,1]);
%             delete([hb, hc]);
%             plh{1}(2) = ha;
%          	axes(AxTuneY)
%             delete(plh{2}(2));
%           	[ha, hb, hc] = shadedplot(linspace(1, numel(PC2), size(A,1)), mean(A,2)-std(A'), mean(A,2)+std(A'), [0.5,0.5,1]);
%             delete([hb, hc]);
%             plh{2}(2) = ha;
%             uistack(plh{1}(1),'top');
%             uistack(plh{2}(1),'top');
        end
      	if InterpOn == 0
            set(plh{1}(pc1+1), 'Ydata', mean(OrientDiffMat{t}'));
            set(plh{2}(pc2+1), 'Ydata', mean(OrientDiffMat{t}'));
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
    print('-painters', '-dsvg', fullfile(SaveDir, sprintf('OrientationTuning_%s_%s_ch%03d_cell%d.svg',Subject, SessionDate, Channel, Cell)))
    export_fig(fullfile(SaveDir, sprintf('OrientationTuning_%s_%s_ch%03d_cell%d.png',Subject, SessionDate, Channel, Cell)), '-png');
end
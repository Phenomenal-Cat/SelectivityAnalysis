
%============================= MF3D_RunAnovas.m ===========================
% This script loads NeuroMat structures (created by MF3D_GenerateSDFs.m)
% and performs N-way ANOVAs per cell, in order to determine which factors
% were most strongly encoded by the neural population.
%
%==========================================================================

Subject     = 'Matcha'; %

Append      = [];
if ismac, Append = '/Volumes'; end
SaveDir     = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs/StereoFaces/', Subject);
SavePopDir  = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Population/');
load(fullfile(SaveDir, sprintf('SDF_%s.mat',Subject)));

RespWindow  = [80, 150];
BaselineWin = [-80, 20];

%% ================== PLOT AVERAGE VISUAL RESPONSE FOR EVERY CELL
fh  = figure('position',get(0,'screensize'));
axh = tight_subplot(1, numel(NeuroMat), 0.04,0.04,0.04);
for D = 1:numel(NeuroMat)
    axes(axh(D));
    imagesc(NeuroMat(D).HistBins, 1:NeuroMat(D).NoCells, squeeze(mean(NeuroMat(D).SDF,2)));
    hold on;
    plot([0,0],ylim,'-w','linewidth',2);
    if D == 1
        ylabel('Cell #');
    else
        set(axh(D),'yticklabel',[]);
    end
    title(NeuroMat(D).Session);
end
set(axh,'clim', [0 50]);
% colorbar;

%% ================== PLOT AVERAGE RESPONSE FOR EVERY STIMULUS
fh  = figure('position',get(0,'screensize'));
axh = tight_subplot(1, numel(NeuroMat), 0.04,0.04,0.04);
for D = 4:numel(NeuroMat)
    ResponseBins        = find(NeuroMat(D).HistBins > RespWindow(1) & NeuroMat(D).HistBins < RespWindow(2));
    BaselineBins        = find(NeuroMat(D).HistBins > BaselineWin(1) & NeuroMat(D).HistBins < BaselineWin(2));
    BaselineResponse    = squeeze(mean(NeuroMat(D).SDF(:,:,BaselineBins),3));
    MeanResponse        = squeeze(mean(NeuroMat(D).SDF(:,:,ResponseBins),3));
    MeanMean          	= mean(MeanResponse');
    [~,CellOrder]     	= sort(MeanMean,'descend');
    BaselineSubMean     = MeanResponse(CellOrder,:)-BaselineResponse(CellOrder,:);
    
%     for f = 1:numel(NeuroMat(D).Params.Factors)
%         [~,StimOrder] = sort(NeuroMat(D).Params.ConditionMatrix(:,f),'descend');
%         
%     end
    
%     axh = tight_subplot(numel(NeuroMat(D).Params.Scales), numel(NeuroMat(D).Params.Distances), 0.01,0.01,0.01);
%     i = 1;
%     for F1 = 1:numel(NeuroMat(D).Params.Scales)
%         for F2 = 1:numel(NeuroMat(D).Params.Distances)
%             StimIndices{F1,F2} = find(ismember(NeuroMat(D).Params.ConditionMatrix(:,[4,5]), [F1, F2],'rows'));
%             FactorMat{F1,F2}    = MeanResponse(CellOrder,StimIndices{F1,F2})-BaselineResponse(CellOrder,StimIndices{F1,F2});
%             axes(axh(i));
%             imagesc((0:1)+(F1-1), (0:1)+(F2-1), FactorMat{F1,F2});
%             hold on;
%             i = i+1;
%         end
%     end
%     axes(axh(D))
%     imagesc();
%     xlabel('Stim #');
%   	if D == 1
%         ylabel('Cell #');
%     else
%         set(axh(D),'yticklabel',[]);
%     end
%     title(NeuroMat(D).Session);
%     
    
    %================ PLOT POPULATION FIGURE
  	PopMean = mean(BaselineSubMean(1:60,:));
    SEM     = std(BaselineSubMean(1:60,:))/sqrt(60);
    
    fhpa = figure('position',get(0,'screensize')); 
    axh(1) = subplot(2,1,1);
    bar(PopMean);
    hold on
    errorbar(1:numel(PopMean), PopMean, SEM, '.b');
    grid on;
    box off;
    axis tight
	title(sprintf('Population average %s %s (%d-%dms)', NeuroMat(D).Subject, NeuroMat(D).Session, RespWindow),'fontsize',16);
    axh = subplot(2,1,2);
    imagesc(NeuroMat(D).Params.ConditionMatrix(:,2:end)');
    set(gca,'ytick',1:numel(NeuroMat(D).Params.Factors)-1, 'yticklabel', NeuroMat(D).Params.Factors(1:end-1));
    export_fig(fullfile(SavePopDir, sprintf('%s_%s_%d-%dms.png', NeuroMat(D).Subject, NeuroMat(D).Session, RespWindow)),'-png');
    close(fhpa);

    %================ GENERATE ORIENTATION TUNING MATRICES
    for cell = 1:size(NeuroMat(D).PSTH, 1)
        for az = 1:numel(NeuroMat(D).Params.Azimuths)
            for el = 1:numel(NeuroMat(D).Params.Elevations)
                StimIndx = find(ismember(NeuroMat(D).Params.ConditionMatrix(:,[2,3]), [az,el], 'rows'));
                Data = NeuroMat(D).PSTH(cell, StimIndx, ResponseBins);
                NeuroMat(D).OrientationMat(cell,el,az) = mean(Data(:));
            end
        end
    end
    
    %================== Run N-way ANOVA on all cells
%  	AnovaData = [];
%     for cell = 1:CellsToInclude
%         AnovaData = [AnovaData, BaselineSubMean(cell,:)];
%     end
%     if ~isfield(NeuroMat(D).Params, 'CondMatCol')
%         NeuroMat(D).Params.Factors = NeuroMat(D).Params.Factors([5,2,1,3,4]);
%     else
%         NeuroMat(D).Params.Factors = NeuroMat(D).Params.Factors([5,2,1,3,4]);
%     end
    Groups = {};
    NeuroMat(D).Stats.InlcudedFactors = {};
    DepthIndx = ~cellfun(@isempty, strfind(NeuroMat(D).Params.Factors, 'Depth'));
    if any(DepthIndx) && ~isfield(NeuroMat(D).Params, 'Depth')                                                              % If depth profile was a variable...
        NeuroMat(D).Params.Depth = {'Concave','Flat','Convex'};                                                             % Specify level labels
        if min(NeuroMat(D).Params.ConditionMatrix(:,find(DepthIndx))) == -1
            NeuroMat(D).Params.ConditionMatrix(:,find(DepthIndx)) = NeuroMat(D).Params.ConditionMatrix(:,find(DepthIndx))+2;    % Adjust range from -1:1 to 1:3
        end
    end
    for F = 1:numel(NeuroMat(D).Params.Factors)
        if numel(unique(NeuroMat(D).Params.ConditionMatrix(:,F))) > 1
            NeuroMat(D).Stats.InlcudedFactors{end+1} = NeuroMat(D).Params.Factors{F};
            Groups{end+1}   = eval(sprintf('NeuroMat(D).Params.%s(NeuroMat(D).Params.ConditionMatrix(:,F));', NeuroMat(D).Params.Factors{F}));
        end
    end
% 	NeuroMat(D).Stats.InlcudedFactors{end+1} = 'Cell';
%     [p,tbl,stats,terms] = anovan(AnovaData, Groups, 'random', numel(NeuroMat(D).Stats.InlcudedFactors),'model',2,'varnames', NeuroMat(D).Stats.InlcudedFactors);
%     
    
    %================== Run N-way ANOVA on EACH cell
    wbh = waitbar(0);
    for cell = 1:size(BaselineSubMean,1)
        waitbar(cell/size(BaselineSubMean,1), wbh, sprintf('Running ANOVA for cell %d of %d...', cell, size(BaselineSubMean,1)));
        [p,tbl,stats,terms] = anovan(BaselineSubMean(cell,:), Groups,'model',2,'varnames', NeuroMat(D).Stats.InlcudedFactors, 'display','off');
        NeuroStats(cell).p      = p;
        NeuroStats(cell).stats  = stats;
        NeuroStats(cell).table  = tbl;
        NeuroStats(cell).groups = Groups;
        PMatrix(:, cell)        = p;
        NeuroStats(cell).Labels	= NeuroStats(cell).table(2:end-2,1);
        
        FactorIndx = [1,2,5];
        if any(p(FactorIndx)<0.05)
            NeuroStats(cell).OrientTuned = 1;
        else
            NeuroStats(cell).OrientTuned = 0;
        end
            
        
        close all;
    end
    delete(wbh)
    
    TunedCellIndx   = find([NeuroStats.OrientTuned]==1);
    OrientMeanMat   = squeeze(mean(NeuroMat(D).OrientationMat(TunedCellIndx,:,:)));
    OrientSDMat     = squeeze(std(NeuroMat(D).OrientationMat(TunedCellIndx,:,:)));
    
    %============== Plot stats matrix
    ScreenRes = get(0,'screensize');
    fh = figure('position',ScreenRes./[1,1,1,2]);
    imagesc(ones(size(PMatrix))-PMatrix);
    set(gca,'yticklabel', NeuroStats(1).table(2:end-2,1), 'clim', [0.95, 1]);
    cbh = colorbar;
    xlabel('Cell #');
    ylabel('Variable');
    colormap hot;
    title(sprintf('%s %s - ANOVA p-values', NeuroMat(D).Subject, NeuroMat(D).Session),'fontsize', 16);
	export_fig(fullfile(SavePopDir, sprintf('%s_%s_ANOVA_pvals.png', NeuroMat(D).Subject, NeuroMat(D).Session)),'-png');
    saveas(fh, fullfile(SavePopDir, sprintf('%s_%s_ANOVA_pvals.fig', NeuroMat(D).Subject, NeuroMat(D).Session)), 'fig');
    close(fh);
    
    clear PMatrix
end

save(fullfile(SaveDir, sprintf('AnovaStats_%s.mat',Subject)), 'NeuroStats','-v7.3');




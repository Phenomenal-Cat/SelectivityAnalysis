
%========================== MF3D_GenerateSDFs.m ===========================

Subject     = 'Spice';
Dates       = {'20160620','20160621','20160622','20160623','20160624'};

% Subject     = 'Matcha';
% Dates       = {'20160613','20160614','20160615','20160616','20160617','20160719','20160720','20160721','20160722'};
% ExpType     = [1,1,2,3,3,4,4,4,5];

% Subject     = 'Avalanche';
% Dates   	= {'20160627','20160628','20160629','20160630','20160701','20160712','20160713','20160714','20160715'};
% ExpType     = [1,1,2,3,3,4,4,5,5];


RespWindow  = [80, 150];
BaselineWin = [-80, 20];

Append  = [];
if ismac, Append = '/Volumes'; end
SaveDir = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs/StereoFaces/', Subject);
SavePopDir = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Population/');
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end
    
for D = 1:numel(Dates)
    
    %============= Set directories
    TimingData              = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/',sprintf('StimTimes_%s_%s.mat', Subject, Dates{D}));
    ProcessedSessionData    = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs/StereoFaces/',Subject,Dates{D},sprintf('%s_%s.mat', Subject, Dates{D}));
    load(TimingData)
    load(ProcessedSessionData);
    
    %============= Add data to structure
  	NeuroMat(D).Subject     = Subject;
    NeuroMat(D).Session     = Dates{D};
    NeuroMat(D).NoCells     = size(AllSpikes,1);
    NeuroMat(D).Params      = Params;
    NeuroMat(D).HistWidth   = 1;
    NeuroMat(D).ChIndx      = ChIndx;
    if isfield(Params, 'PreTime')
        NeuroMat(D).HistBins    = (-Params.PreTime*10^3):NeuroMat(D).HistWidth:(Params.PostTime*10^3);
    else
        NeuroMat(D).HistBins    = -100:1:400;
    end

    %============= Make PSTH and SDF matrices
    for cell = 1:size(AllSpikes,1)
        for stim = 1:size(AllSpikes,2)
           	AllSpikeTimes = [];
            for t = 1:size(AllSpikes,3)
                if ~isnan(AllSpikes{cell,stim,t})
                    AllSpikeTimes = [AllSpikeTimes; AllSpikes{cell,stim,t}];
                end
            end
            NeuroMat(D).PSTH(cell, stim, :)    = hist(AllSpikeTimes, NeuroMat(D).HistBins)*(10^3/NeuroMat(D).HistWidth)/Stim.Repetitions(stim);
            NeuroMat(D).SDF(cell, stim, :)     = msdf(squeeze(NeuroMat(D).PSTH(cell, stim, :)),'Gauss',5);
        end
        
    end
    
end

save(fullfile(SaveDir, sprintf('SDF_%s.mat',Subject)), 'NeuroMat','-v7.3');


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
for D = 1:numel(NeuroMat)
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

    
    %================== Run N-way ANOVA on all cells
 	AnovaData = [];
    for cell = 1:CellsToInclude
        AnovaData = [AnovaData, BaselineSubMean(cell,:)];
    end
    if ~isfield(NeuroMat(D).Params, 'CondMatCol')
        NeuroMat(D).Params.Factors = NeuroMat(D).Params.Factors([5,2,1,3,4]);
    else
        NeuroMat(D).Params.Factors = NeuroMat(D).Params.Factors([5,2,1,3,4]);
    end
    Groups = {};
    NeuroMat(D).Stats.InlcudedFactors = {};
    for F = 1:numel(NeuroMat(D).Params.Factors)
        if numel(unique(NeuroMat(D).Params.ConditionMatrix(:,F))) > 1
            NeuroMat(D).Stats.InlcudedFactors{end+1} = NeuroMat(D).Params.Factors{F};
            Groups{end+1}   = eval(sprintf('NeuroMat(D).Params.%s(NeuroMat(D).Params.ConditionMatrix(:,F));', NeuroMat(D).Params.Factors{F}));
        end
    end
	NeuroMat(D).Stats.InlcudedFactors{end+1} = 'Cell';
    [p,tbl,stats,terms] = anovan(AnovaData, Groups, 'random', numel(NeuroMat(D).Stats.InlcudedFactors),'model',2,'varnames', NeuroMat(D).Stats.InlcudedFactors);
    
end
set(axh,'clim', [0 50]);
colormap hot


%================== Run N-way ANOVA on EACH cell
for cell = 1:size(BaselineSubMean,1)
    [p,tbl,stats,terms] = anovan(BaselineSubMean(cell,:), Groups,'model',2,'varnames', NeuroMat(D).Stats.InlcudedFactors);
    NeuroStats(cell).p      = p;
    NeuroStats(cell).stats  = stats;
    NeuroStats(cell).table  = tbl;
end


%========================= MF3D_PlotAllTuning.m ===========================
% Plot single cell tuning to head orientation, stereoscopic depth and 
% facial expression based on the first week of pilot data for all monkeys.
%
%==========================================================================

Subject     = 'Avalanche';
Output      = 'gif';

switch Subject
    case 'Avalanche'
        Dates      	= {'20160627','20160628','20160629','20160630','20160701'};
        Channels  	= {[100, 89, 125, 123],...          % Specify Session 1 channels of interest
                        };
        Cell        = 1;

    case 'Matcha';
        Dates       = {'20160613','20160614','20160615','20160616','20160617',};
        Channels    = {};
        
    case 'Spice'
      	Dates       = {'20160620','20160621','20160622','20160623','20160624'};
        Channels    = {};
        
    otherwise
        error('unknown subject ''%s''!\n', Subject)
end

%=================== Lopp through cells
for d = 1:numel(Dates)
    for ch = 1:numel(Channels)
        if d <= 3
            MF3D_PlotResponses_Orientations(Subject, Dates{d}, Channels{d}(ch), Cell, Output);
        elseif ismember(d,[4,5])
            MF3D_PlotResponses_Expressions(Subject, Dates{d}, Channels{d}(ch), Cell, Output);
            MF3D_PlotResponses_Depth(Subject, Dates{d}, Channels{d}(ch), Cell, Output);
        elseif ismember(d,[6,7,8]);
            
        end
    end
end
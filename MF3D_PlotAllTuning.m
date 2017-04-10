

%================== Plot head orientation, depth, and expression tuning
Subject     = 'Avalanche';
Dates      	= {'20160701'};
Channels  	= [100, 89, 125, 123];
Cell        = 1;
Output      = 'gif';


for d = 1:numel(Dates)
    for ch = 1:numel(Channels)
        MF3D_PlotResponses_Expressions(Subject, Dates{d}, Channels(ch), Cell, Output);
        MF3D_PlotResponses_Depth(Subject, Dates{d}, Channels(ch), Cell, Output);
        
    end
end
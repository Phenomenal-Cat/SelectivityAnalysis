
%======================== MF3D_PlotCategoryResponses.m ====================



Subject = 'Matcha';
Date    = '20160613';

Append  = [];
if ismac, Append = '/Volumes'; end
TimingData              = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/',sprintf('StimTimes_%s_%s.mat', Subject, Date));
ProcessedSessionData    = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/PSTHs/',Subject,sprintf('%s_%s.mat', Subject, Date));
load(TimingData)
load(ProcessedSessionData);

Cell = 7;
Cond = 1;

BinData{Cell, Cond}
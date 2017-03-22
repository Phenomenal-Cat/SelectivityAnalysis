
% PlotAllERPs.m

Subject     = 'Matcha';
Date        = '20160613';
ExpType     = 'StereoFaces';
NoChannel 	= 128;

for Channel = 1:NoChannels
    PlotERPs(Subject, Date, ExpType, Channel);
end
function PlotSpect(Subject, Date, ExpType, Channel)

%================================ PlotSpect.m ==============================
% This function loads raw LFP data form the StereoFaces experiments and
% plots spectrograms and power bands.
%
%
%==========================================================================


RootDir     = '/Volumes/APM_128GB/NIH_Postdoc/PilotPhysiology/LFP/';
MatDir      = fullfile(RootDir, Subject, Date);
MatFiles    = wildcardsearch(MatDir, sprintf('LFPproc_*ch%d.mat', Channel));

addpath(genpath('/Users/aidanmurphy/Documents/Toolboxes/chronux_2_12'));


params.tapers   = [5, 9];
params.Fs       = 1000;
params.fpass    = [0, 200];
params.trialave = 1;
params.err      = 0;
params.pad      = 0;
movingwin       = [0.1, 0.025];



Params.PreTime         = 0.1;
Params.PostTime        = 0.5;
Params.StimDur         = 0.3;
Params.NoiseThresh     = 2000;

%================ Plot spectrograms
for f = 1:4%numel(Factors)
    Fhs(f)     	= figure('position',get(0,'ScreenSize'));                                   
    axh         = tight_subplot(4, 4, 0.05, 0.05, 0.05);

    for el = 1:numel(Factor{f})
        axes(axh(el));
        [S,t,fq]    = mtspecgramc(ElLFPall{f,el}', movingwin, params);
        FqTicks   	= linspace(params.fpass(1), params.fpass(2), numel(fq));
        TimeTicks  	= linspace(-Params.PreTime, Params.PostTime, numel(t));
        Data{f,el}  = log(S);
        
        imagesc(TimeTicks, FqTicks, Data{f,el}');
        axis xy;
        hold on;
        plot([0 0], get(gca, 'ylim'), '--r');
        xlabel('Time (s)','fontsize', 18);
        ylabel('Frequency (Hz)','fontsize', 18);
        title(sprintf('%s = %d', Factors{f}, Factor{f}(el)), 'fontsize', 18);
    end
%     suptitle(sprintf('%s %s %s %s', Subject, Date, ExpType, Factors{f}));
end

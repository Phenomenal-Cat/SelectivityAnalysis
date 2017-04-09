

%============================= MF3D_PlotGaze.m ============================
% This function loads the gaze data for a single session date and performs
% some basic analysis to check whether fixation behaviour differed across
% experimental conditions.
%
%==========================================================================

Subject = 'Matcha';
Date    = '20160613';

TimingDir   = fullfile('/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/Timing/StereoFaces');
TimingFile  = fullfile(TimingDir, sprintf('StimTimes_%s_%s.mat', Subject, Date));
GazeDir     = fullfile('/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/TDT_converted/',Subject,Date);
GazeFiles  	= wildcardsearch(GazeDir, '*eyeSignal.mat');

%=============== Load gaze data from all blocks
Eye.Signal  = [];
Eye.Times   = [];
for f = 1:numel(GazeFiles)
    load(GazeFile{f})
    
    [eyeCh,eyeSig,eyeSigPerSample] = size(eyeCodesAll);                                                     % Check matrix dimensions
    for n = 1:eyeCh                                                                                         % For each channel...
        AllEyeCodes(n,:) = reshape(permute(eyeCodesAll(n,:,:),[3,2,1]),[1,numel(eyeCodesAll(n,:,:))]);      % Reshape eye signal
    end
    eyeTimes            = linspace(0, eyeTimesAll(end), length(AllEyeCodes));
    Eye.Signal        	= [Eye.Signal, AllEyeCodes];
    if f > 1
        Eye.Times   	= [Eye.Times, eyeTimes+Eye.Times(end)+eyeTimes(2)];
    else
        Eye.Times     	= eyeTimes;
    end
end

%=============== Gather gaze data by stimulus ID
TrialPeriod     = [-0.1, 0.4]; 
StimSamples     = round(diff(TrialPeriod)*eyeSampleRate);
load(TimingFile)
for S = 1:numel(Stim.Onsets)
    EyeSigX{S} = [];
    EyeSigY{S} = [];
    EyeSigP{S} = [];
    for t = 1:numel(Stim.Onsets{S})
        StimTime    = Stim.Onsets{S}(t)+TrialPeriod;
        StimSample1 = find(Eye.Times > Stim.Onsets{S}(t)+TrialPeriod(1));
        EyeSigX{S}     = [EyeSigX{S}; Eye.Signal(1, StimSample1(1)+(0:StimSamples-1))];
        EyeSigY{S}     = [EyeSigY{S}; Eye.Signal(2, StimSample1(1)+(0:StimSamples-1))];
        EyeSigP{S}     = [EyeSigP{S}; Eye.Signal(3, StimSample1(1)+(0:StimSamples-1))];
    end
    
    %============= Plot gaze distribution over stimuli
    axes(AxH(S));
    [Im,Cm,ImAlpha]	= imread(Params.FullFilenames{S});
    ImH(S,1)        = imagesc(imresize(Im(:,1:size(Im,2)/2, :), size(Im)));
    alpha(ImH(S,1), imresize(ImAlpha(:,1:size(ImAlpha,2)/2, :), size(ImAlpha)));
    
    HistH = ndhist(reshape(EyeSigX{S}',[1,numel(EyeSigX{S})]), reshape(EyeSigY{S}',[1,numel(EyeSigY{S})]));
    HistH = hist3([reshape(EyeSigX{S}',[numel(EyeSigX{S}),1]), reshape(EyeSigY{S}',[numel(EyeSigY{S}),1])]);
    ImH(S,2)        = imagesc(HistH);
    alpha(ImH(S,2), 0.5);
    axis equal tight
    box off
    set(gca,'color', [0.5, 0.5, 0.5]);
    
end



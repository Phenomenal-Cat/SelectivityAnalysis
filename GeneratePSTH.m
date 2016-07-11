

%============================= GeneratePSTH.m =============================

ExpName     = 'StereoFaces';
Subject     = 'Spice';
Dates       = {'20160620','20160621'};
TimingDir   = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/Timing';
LFPDir      = '/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/TDT_converted';

PreStimulus     = 0.1;          % Start time window (seconds before stim on)
PostStimulus    = 0.2;          % End time window (seconds after stim off)

%================== GET ALL SPIKE TIMES IN WINDOW
for D = 1:numel(Dates)
    TimingFile  = fullfile(TimingDir, sprintf('StimTimes_%s_%s.mat', Subject, Dates{D}));
    load(TimingFile);

   	StimDuration    = unique([ExpParam.Trial_Time]);
    ISI             = unique([ExpParam.Blank_btw_stim]);
    InitialFix      = unique([ExpParam.Initial_Fix_time]);
    
    LFPfiles     	= wildcardsearch(fullfile(LFPDir,Subject,Dates{D}), ['*',ExpName,'*LFP*.mat']);    
    for Ch = 1:numel(LFPfiles)
        ChFile      = wildcardsearch(fullfile(LFPDir,Subject,Dates{D}), ['*',ExpName,'*LFP-ch',num2str(Ch),'.mat']); 
        if numel(ChFile) > 1
            disp(ChFile);
            error('Multiple files match LFP channel number!');
        end
        load(ChFile{1});
        if ~exist('TrialCount','var')
            TrialCount = ones(numel(AllChannels), numel(Stim.Onsets));      % Start channel count
        end
        for S = 1:numel(Stim.Onsets)                                        % For each unique stimulus...
            for T = 1:numel(Stim.Onsets{S});                                % For each repetition...
                TimeWindow = [-PreStimulus, StimDuration+PostStimulus]+Stim.Onsets{S}(T);
                % <<<<< X = all spike times 
                PSTHSpikeTimes{S}{TrialCount(S)} = X(Ch, find(X >= TimeWindow(1) && X <= TimeWindow(2)));
                TrialCount(Ch, S) = TrialCount(Ch, S)+1;
            end
        end
    end
end

%================== BIN SPIKES AND SMOOTH




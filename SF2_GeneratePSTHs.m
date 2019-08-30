
%======================== SF2_GeneratePSTHs.m =============================
% Generate PSTHs for all data collected for StereoFaces2 experiments.

RootDataDir = '/Volumes/Seagate Backup 4/NIH_Neurophys/StereoFaces_2';
Subjects    = {'Spice','StevieRay','Mochi','Wasabi'};
NoChannels  = [64, 128, 128, 64];
ExpNames    = {'FingerPrint','SizeDistance','SizeDistance_Movies','StereoShape'};
SortSDs     = [27,30,33,36,40,50,55];
SDIndx      = 5;
StimParamsDir = fullfile(RootDataDir, 'StimParams');
TimeWindow  = [-0.1, 0.4];

if ~exist('wildcardseacrh.m','file')
    addpath(fullfile(RootDataDir, 'SF2_Matlab/Subfunctions'));
end

for S = 1:numel(Subjects)
    
    %=========== Find all data files for current subject
    SpikesDir   = fullfile(RootDataDir, 'WaveClusSorted', Subjects{S});
    Temp        = dir(SpikesDir);
    Temp        = {Temp.name};
    SubjectDates{S} = Temp(cellfun(@isempty, strfind(Temp, '.')));
    EventsDir   = fullfile(RootDataDir, 'TDT_converted', Subjects{S});
    Temp        = dir(EventsDir);
    Temp        = {Temp.name};
    BlockEvents{S} = Temp(cellfun(@isempty, strfind(Temp, '.')));
    
    for d = 6:numel(SubjectDates{S})
        fprintf('Processing session: %s %s (%d/%d)...\n', Subjects{S}, SubjectDates{S}{d}, d, numel(SubjectDates{S}));
        
        %============ Find necessary .mat files
        SpikeFile   = fullfile(SpikesDir, SubjectDates{S}{d}, sprintf('SD%d_%s_sorted.mat', SortSDs(SDIndx), SubjectDates{S}{d}));
        EventFiles  = wildcardsearch(EventsDir, sprintf('TDTconv_%s_%s*.mat', Subjects{S}, SubjectDates{S}{d}));
        if isempty(EventFiles)
            error('ERROR: No event .mat files found for %s on %s!', Subjects{S}, SubjectDates{S}{d});
        end
        if ~exist(SpikeFile,'file')
            error('ERROR: No spike sorted .mat file found for %s on %s!', Subjects{S}, SubjectDates{S}{d});
        end
        
        %============ Load data and check
        load(SpikeFile);
        NoBlocks    = numel(NeuroStruct);
        NoCells     = size(NeuroStruct(1).cells,1);
        BlockDurs   = [NeuroStruct.blocklength]/1000;
        
        %========= Find which experiment type was run
        for b = 1:numel(EventFiles)
            [~,EvntFile]    = fileparts(EventFiles{b});
            EvntFiles{b}     = EvntFile;
        end
        for exp = 1:numel(ExpNames)
            ExpType(~cellfun(@isempty, strfind(EvntFiles, ExpNames{exp}))) = exp;
        end
        AllExpTypes = unique(ExpType);
        
        %========= Load data for same experiment type
        for b = 3:numel(ExpType)
            ExpName = ExpNames{ExpType(b)};
            StimParamsMat = wildcardsearch(StimParamsDir, sprintf('ImParams_%s*.mat', ExpName));
            load(StimParamsMat{numel(StimParamsMat)});                     % Load stimulus image parameters data
            NoStim = ImParams.TotalImages;
            fprintf('Experiment: %s (%d stim)...\n', ExpName, NoStim);
            
            %========= Find matching block for event times and neural data
            %ClockError  = -0.0040;
            ClockError  = 0;
            load(EventFiles{b});
            RecordTime      = datenum(TDTdata.info.starttime, 'HH:MM:SS');
            RecordDuration  = duration(TDTdata.info.duration,'InputFormat','hh:mm:ss');
            if ~isempty(strfind(NeuroStruct(1).block, '_'))
                TimeFormat = 'HH_MM_SS';
                TimeIndx   = 8;
            else
             	TimeFormat = 'HHMMSS';
                TimeIndx   = 5;
            end
            for ns = 1:numel(NeuroStruct)
                NeuroBlockTimes(ns)     = datenum(NeuroStruct(ns).block((end-TimeIndx):end), TimeFormat)+ClockError;
                NeuroBlockDurs(ns)      = seconds(NeuroStruct(ns).blocklength/1000);
            end
            TimeDiffMin    	= (NeuroBlockTimes - RecordTime)*24*60;
            DurDiffMin      = (NeuroBlockDurs-RecordDuration);
            [MatchDiff1, MatchIndx1] 	= min(abs(TimeDiffMin));
            [MatchDiff2, MatchIndx2] 	= min(abs(DurDiffMin));
            if MatchIndx1 ~= MatchIndx2
                MatchIndx1      = find(TimeDiffMin < 0, 1, 'last');
            end
            if MatchIndx1 ~= MatchIndx2
                figure; 
                plot(NeuroBlockTimes, NeuroBlockDurs, '.--r','MarkerSize', 40);
                hold on;
                plot(RecordTime, RecordDuration, 'ob','MarkerSize', 30);
                xlabel('Block start time','fontsize', 18)
                ylabel('Block duration','fontsize', 18)
                MatchIndx2 = MatchIndx1;
%                 error('Closest time and duration matches are not the same block!');
            end
            
           	%========= Get stim times and extract spikes
            AllTimes(d,b).Stim  = SF2_GetStimTimes(TDTdata, NoStim);  
            SpikeTimes          = SF2_GetSpikesByStim(NeuroStruct(MatchIndx1), AllTimes(d,b).Stim.Times, TimeWindow);
            
            SpikeDir            = fullfile(RootDataDir, 'SpikeTimes', Subjects{S});
            if ~exist(SpikeDir)
                mkdir(SpikeDir);
            end
            SpikeFilename       = fullfile(SpikeDir, sprintf('%s_%s_%s_%d.mat', Subjects{S}, ExpName, SubjectDates{S}{d}, b));
            if ~exist(SpikeFilename, 'file')
                save(SpikeFilename, 'AllTimes', 'SpikeTimes');
            end
        end
        pause(1);
    end
end
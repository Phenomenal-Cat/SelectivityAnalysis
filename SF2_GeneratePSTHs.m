
%======================== SF2_GeneratePSTHs.m =============================
% Generate PSTHs for all data collected for StereoFaces2 experiments.

RootDataDir = '/Volumes/Seagate Backup 4/NIH_Neurophys/StereoFaces_2';
Subjects    = {'Spice','StevieRay','Mochi','Wasabi'};
NoChannels  = [64, 128, 128, 64];
ExpNames    = {'FingerPrint','SizeDistance','SizeDistance_Movies','StereoShape','Stereo'};
SortSDs     = [27,30,33,36,40,50,55];
SDIndx      = 5;
StimParamsDir = fullfile(RootDataDir, 'StimParams');
TimeWindow  = [0.1, 0.4];

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
    
    for d = 1:numel(SubjectDates{S})
        
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
        
        
        for b = 1:numel(EventFiles)
           
            
            for n = 1:numel(ExpNames)
                [~,file] = fileparts(EventFiles{b});
                if ~isempty(strfind(lower(file), lower(ExpNames{n})))
                    AllTimes(d,b).BlockType = n;
                end
            end
            
            
            %========= Find which experiment type was run
            [~,EvntFile] = fileparts(EventFiles{b});
            for exp = 1:numel(ExpNames)
                if ~isempty(strfind(EvntFile, ExpNames{exp}))
                    ExpName = ExpNames{exp};
                end
            end
            StimParamsMat = wildcardsearch(StimParamsDir, sprintf('ImParams_%s*.mat', ExpName));
            load(StimParamsMat{1});                     % Load stimulus image parameters data
            
            %========= Find matching block for event times and neural data
            load(EventFiles{b});
            RecordTime      = datenum(TDTdata.info.starttime, 'HH:MM:SS');
            RecordDuration  = duration(TDTdata.info.duration,'InputFormat','hh:mm:ss');
            for ns = 1:numel(NeuroStruct)
                NeuroBlockTimes(ns)     = datenum(NeuroStruct(ns).block((end-8):end), 'HH_MM_SS');
                NeuroBlockDurs(ns)      = seconds(NeuroStruct(ns).blocklength/1000);
            end
            TimeDiffMin    	= (NeuroBlockTimes - RecordTime)*24*60;
            DurDiffMin      = (NeuroBlockDurs-RecordDuration);
            [~,MatchIndx1] 	= min(abs(TimeDiffMin));
            [~,MatchIndx2] 	= min(abs(DurDiffMin));
            if MatchIndx1 ~= MatchIndx2
                error('Closest time and duration matches are not the same block!');
            end
            
           	%========= Get stim times and extract spikes
            AllTimes(d,b).Stim  = SF2_GetStimTimes(TDTdata);  
            SpikeTimes          = SF2_GetSpikesByStim(NeuroStruct(MatchIndx1), AllTimes(d,b).Stim.Times, TimeWindow);
            
            pause(1);
        end
        pause(1);
    end
end
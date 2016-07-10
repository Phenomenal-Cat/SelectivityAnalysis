
%========================= AnalyzeStereoFaces.m ===========================


SubjectID       = 'Matcha';
DateStrings     = {'20160613','20160614','20160615','20160616','20160617'};


NeuralDataDir 	= fullfile('/Volumes/PROCDATA/murphya/Physio/WaveClusSorted/',SubjectID);
ProcDataDir     = '/Volumes/PROCDATA/murphya/Physio/StereoFaces';


for d = 1:numel(DateStings)
    ConditionsFile{d}   = fullfile(ProcDataDir, sprintf('StereoFaces_Conditions_%s.mat', DateStrings{d}));
    TimingFile{d}       = fullfile(ProcDataDir, sprintf('StimTimes_%s_%s.mat', SubjectID, DateStrings{d}));
 	load(ConditionsFile{d});
    load(TimingFile{d});
    
%     [QNX, PD] = GetStimulusTimes(TDTdir)


    
    
end


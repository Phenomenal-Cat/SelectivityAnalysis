
%============================ GetAllStimTimes.m ===========================
% This is a wrapper script that allows the user to select which subjects 
% and session dates to get stimulus times for using GetStimulusTimes.m

% ExpName = 'StereoFaces';
ExpName = 'FingerPrint';
ExpType = 1;

[~,CompName] = system('hostname');  
if strcmpi(CompName(1:end-1), 'Aidans-MacBook-Pro.local')
  	DataDir  = fullfile('/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/TDT_converted');
else
    DataDir  = fullfile('/Volumes/rawdata/murphya/Physio/TDT_converted/');
end

Subjects = dir(DataDir);
Subjects = {Subjects(cellfun(@isempty, strfind({Subjects.name},'.'))).name};
[s,v] = listdlg('ListString',Subjects,'PromptString','Select subjects');
for S = s
    Dates = dir(fullfile(DataDir, Subjects{S}));
    Dates = {Dates(cellfun(@isempty, strfind({Dates.name},'.'))).name};
    [d,v] = listdlg('ListString',Dates,'PromptString','Select session dates');
    for D = d
        GetStimulusTimes(ExpName, Subjects{S}, Dates{D}, ExpType, 1);
    end
end

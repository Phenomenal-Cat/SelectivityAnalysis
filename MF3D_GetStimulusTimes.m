% function [Stim, QNX, PD] = GetStimulusTimes(ExpName, SubjectID, DateString, ExpType, Verbose)

%======================== GetStimulusTimes.m ==============================
% This function retreives the timing and stimulus identity information for
% every stimulus presentation that occured within the TDT block specified
% by the input directory.
%
% INPUTS:   ExpName:    string specifying experiment filename
%           SubjectID:  subject identifier string (i.e. monkey name)
%           DateStr:    date string in YYYYMMDD format
%           Verbose:    plot data and print analyses of timing data
%
% OUTPUTS:  Stim:   a 1xN cell array, where N is the maximum stimulus ID 
%                   number. Each cell contains all the stimulus onset times
%                   from successfully completed stimulus presentations.
% HISTORY:
%   13/06/2016 - Written by APM (murphyap@mail.nih.gov)
%==========================================================================

% if nargin == 0
    ExpName     = 'FingerPrint60';
    SubjectID   = 'Spice';
    DateString  = '20170503';
    ExpType     = 6;
    OnlyCompletedObs = 0;           % 1 = only include trials from completes obs blocks
    Verbose     = 1;
% end

%============== SET PATHS FOR CURRENT SYSTEM
[~,CompName] = system('hostname');  
if strcmpi(CompName(1:end-1), 'Aidans-MacBook-Pro.local')
  	OutputFile  = fullfile('/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/Timing/', ExpName, sprintf('StimTimes_%s_%s.mat', SubjectID, DateString));
    TDTdir      = fullfile('/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/TDT_converted/', SubjectID, DateString);
    QNX.dir     = fullfile('/Volumes/Seagate Backup 1/NeuralData/FacePatchPilot/QNX/', SubjectID, DateString);
else
    OutputFile  = fullfile('/Volumes/procdata/murphya/Physio/StereoFaces/Timing', ExpName, sprintf('StimTimes_%s_%s.mat', SubjectID, DateString));
    TDTdir      = fullfile('/Volumes/RAWDATA/murphya/Physio/TDT_converted/', SubjectID, DateString);
    QNX.dir     = fullfile('/Volumes/RAWDATA/murphya/Physio/QNX/', SubjectID, DateString);
end


%======================== READ QNX EVENT DATA ===========================
QNX.codes.Begin         = 1000;     %% the beginning time of a block;
QNX.codes.BeginInBlock  = 2;        %  there are two codeBegin in each block;
QNX.codes.FixOn         = 1001;     %% the time of fix cross on;
QNX.codes.Reward        = 1066;     %% every reward;
QNX.codes.Finish        = 1999;     %% the ending time of block;

PDchannelNo = 7;
if datenum(DateString,'yyyymmdd')>= 736818
    PDchannelNo = 3;
end


%% ======================== READ DGZ FILE DATA ============================
DGZFiles = sort_nat(wildcardsearch(QNX.dir, sprintf('%s_%s_%s*.dgz', SubjectID, DateString, ExpName)));
for i = 1:numel(DGZFiles)                                                       % For each .dgz file...
    fprintf('Loading DGZ %s...\n', DGZFiles{i});
    try
        [ExpParam(i), DGZ(i)] = MF3D_LoadDGZdata(DGZFiles{i});
    catch
        DGZ(i).e_params = [];
    end
end
if isfield(ExpParam, 'Experiment')
    AllBlocks       = {ExpParam.Experiment};
    NoName          = find(cellfun(@isempty,  AllBlocks));
end
TrialsPerObs  	= unique([ExpParam.Trials_Per_Obs]);
if numel(TrialsPerObs) > 1
    fprintf('WARNING: number of trials per obs was not kept constant during this session!\n');
    TrialsPerObs = [ExpParam.Trials_Per_Obs];
elseif isempty(TrialsPerObs)
    if strcmp(ExpName, 'FingerPrint')
        TrialsPerObs = 12;
    elseif strcmp(ExpName, 'StereoFaces')
        TrialsPerObs = 5;
    end
end


%======================== Load all data
PD.Signal   = [];
AllQNXtimes = [];
AllQNXcodes = [];
DataTypes   = {'AnlgSignal','QNX'}; 
for dt = 1:numel(DataTypes)
    Filename = sort_nat(wildcardsearch(TDTdir, ['*',ExpName,'*',DataTypes{dt},'*']));
    if numel(Filename)> 1
        answer = questdlg(sprintf('%d %s files were located. Do you want to merge them?', numel(Filename), DataTypes{dt}), 'Multiple files found','Merge all','Merge some','Abort','Merge some');
        if strcmpi(answer, 'Abort')
            return
        elseif strcmpi(answer, 'Merge some')
            [s,v] = listdlg('ListString',Filename, 'PromptString','Select files to merge','ListSize',[500, 200]);
            Filename = Filename(s);
        end
    elseif numel(Filename)< 1
        error('No files found matching %s!', ['*',ExpName,'*',DataTypes{dt},'*'])
    end

    for f = 1:numel(Filename)
        load(Filename{f});
        
        %=================== READ PHOTODIODE SIGNAL =======================
        if strcmpi(DataTypes{dt}, 'AnlgSignal')
            if isnan(anlgSampleRate)
                anlgSampleRate = 1017.3;                                    % <<<<< Default sample rate used in SCNI rig #2 in 2017
            end
            [anlgCh,anlgSig,anlgSigPerSample] = size(anlgCodesAll);                                                 % Check matrix dimensions
            for n = 1:anlgCh                                                                                        % For each channel...
                AllAnlgCodes(n,:) = reshape(permute(anlgCodesAll(n,:,:),[3,2,1]),[1,numel(anlgCodesAll(n,:,:))]);   % Reshape analog codes
            end
            BlockEndTime(f)     = size(AllAnlgCodes,2)/anlgSampleRate;
            PD.Signal           = [PD.Signal, AllAnlgCodes(PDchannelNo,:)];
            clear AllAnlgCodes;
            
     	%======================= READ QNX STROBES =========================
        elseif strcmpi(DataTypes{dt}, 'QNX')
            if f > 1
                AllQNXtimes = [AllQNXtimes, QNXtimes+sum(BlockEndTime(1:(f-1)))+(1/anlgSampleRate)];
            else
                AllQNXtimes = QNXtimes;
            end
            if strcmp(ExpName, 'StereoFaces') & strcmp(SubjectID, 'Matcha') & strcmp(DateString, '20160615')  	% <<<<<<<<<< TEMPORARY HACK FOR MATCHA FEAR EXPERIMENT
                QNXcodes(QNXcodes>189 & QNXcodes <= 2*189) = QNXcodes(QNXcodes>189 & QNXcodes <= 2*189)-189;    % <<<<<<<<<< TEMPORARY HACK FOR MATCHA FEAR EXPERIMENT
            end
            AllQNXcodes = [AllQNXcodes, QNXcodes];
        end
    end
end


% PD.Signal = CleanDiodeFlicker(PD.Signal);       	% <<<<<<< TEMPORARY FIX FOR PHOTODIODE IN SESSION #1 IN RIG 2


%============== Get experiment specific params
if strcmp(ExpName, 'StereoFaces')
    Params             	= MF3D_GetConditions(ExpType);
else
    Params = [];
end

%===================== ANALYSE PHOTODIODE SIGNAL ==========================
PD = CheckPhotodiode(PD.Signal, anlgSampleRate);                            % Extract photodiode state changes from signal
for p = 1:numel(PD.OnTimes)                                              	% For each photodiode onset...
    QNXindx         = find(AllQNXtimes < PD.OnTimes(p));                  	% Find the preceding code that was received from QNX
    PD.QNXcodes(p)  = AllQNXcodes(QNXindx(end));            
end
NonStimOnsets       = find(ismember(PD.QNXcodes, [QNX.codes.Begin, QNX.codes.FixOn, QNX.codes.Finish]));	% Find all photodiode onsets not related to stimuli
NonStimDiffs        = diff(NonStimOnsets)-1;                            	% Find the number of photodiode onsets between non-stim onsets
if OnlyCompletedObs == 1
    CompletedObs        = find(NonStimDiffs==TrialsPerObs);               	% Find observations for which all trials were completed
    CompletedObsBegin   = NonStimOnsets(CompletedObs);                          
    StimOnCompBlock     = [];
    for t = 1:TrialsPerObs
        StimOnCompBlock	= sort([StimOnCompBlock, CompletedObsBegin+t]);             
    end
elseif OnlyCompletedObs == 0
    StimOnCompBlock = find(PD.QNXcodes>0 & PD.QNXcodes<=numel(Params.Filenames));
end
StimOnTimes         = PD.OnTimes(StimOnCompBlock);                         	% Find the photodiode onset times for all valid stimulus presentations
StimOnIDs           = PD.QNXcodes(StimOnCompBlock);                       	% Find the stimulus ID number for each valid stimulus presentations
Stim.AllStimuli 	= unique(StimOnIDs);                                    % Find all stimulus numbers for which observations were completed
for s = 1:max(Stim.AllStimuli)                                              % For each stimulus that was sucessfully presented
    Stim.Onsets{s} = StimOnTimes(find(StimOnIDs==s));                       % Find all stimulus onset times
    Stim.Repetitions(s) = numel(Stim.Onsets{s});                            % Count how many repetitions of the stimulus were presented
end


%========== Alternative method


%============= Get all stimulus ID numbers from DGZ
for block = 1:numel(DGZ)
    if ~isempty(DGZ(block).e_params)
        for obs = 1:numel(DGZ(block).e_params)
            AllChars = DGZ(block).e_params{obs}(find(cellfun(@ischar, DGZ(block).e_params{obs})));
            AllChars = AllChars(~cellfun(@isempty, strfind(AllChars, 'pic#')));
            for t = 1:numel(AllChars)
                WhtSpcIndx = strfind(AllChars{t},' ');
                QNX.ImageIDs{block}(obs,t) = str2num(AllChars{t}(WhtSpcIndx(2)+1:WhtSpcIndx(3)-1));
            end
        end
        QNX.ValidTrials{block} = find(QNX.ImageIDs{block}~=0);
    end
end

BlockStartIndx          = find(QNXcodes == QNX.codes.Begin);
BlockStartTimes         = QNXtimes(BlockStartIndx);

save(OutputFile, 'Stim', 'ExpParam', 'QNX', 'PD', 'Params');


%======================== PLOT SOME DATA
if Verbose == 1
    figure;
    axh(1) = subplot(3,2,1);
    plot(1:max(Stim.AllStimuli), Stim.Repetitions, '.b');
    xlabel('Stimulus number');
    ylabel('Number of repetitions');
    Ylims = get(gca,'ylim');
    set(gca,'ylim',[0, Ylims(2)]);
    
    axh(2) = subplot(3,2,2);
    hist(PD.Durations);
    xlabel('Photodiode duration (seconds)');
    ylabel('Frequency');
        
    
    suptitle(sprintf('%s %s %s', SubjectID, DateString, 'StereoFaces'));
    
end


% function [ExpParam, DGZ] = LoadDGZdata(DGZfile)
% 
% DGZ                 = dg_read(DGZfile);                                 % Load the .dgz file            
% par_num             = size(DGZ.e_pre,1);                              	% Count how many params
% ExpParam.Experiment = DGZ.e_pre{1}{2};                                  % Get experiment name
% ExpParam.Subject    = DGZ.e_pre{2}{2};                                	% Get subject name
% ExpParam.ADCperDeg  = DGZ.e_params{1}{2};                             	% Get ADC/deg scale
% for j = 3:2:(floor(par_num/2)*2)                                       	% For each parameter...
%     DGZ.e_pre{j}{2}(ismember(DGZ.e_pre{j}{2}, ' ()!?*-:')) = '_';    	% Remove punctuation from parameter names
%     if ~isempty(DGZ.e_pre{j}{2})                                       	% If parameter name cell is not empty...
%         eval(sprintf('ExpParam.%s = %d;', DGZ.e_pre{j}{2}, str2double(DGZ.e_pre{j+1}{2})));
%     else
%         return                                                              
%     end
% end
% %             tmp1 = getEVTParams(DGZ,26,1,0);           	% Use dg utils to get values
% %             tmp2 = getEVTParams(DGZ,26,2,0);
% %             Alt_ADCdegHV{d} = [tmp1, tmp2];
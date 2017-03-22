
%========================= RemoveSEVheader.m ==============================
% This function loads neural data from .sev files and removes header
% information before saving again as in .sev format, ready for spike
% sorting via WaveClus.
%==========================================================================

SubjectID       = 'Matcha';
ExpName         = 'StereoFaces';

DataDir         = fullfile('/spin1/users/murphyap/Rawdata/',SubjectID);
SessionDates    = dir(DataDir);
SessionDates    = {SessionDates(find(cellfun(@isempty, strfind({SessionDates.name}, '.')))).name};

greattic = tic;
for d = 1:numel(SessionDates)
    SEVfiles{d} = wildcardsearch(fullfile(DataDir, SessionDates{d}), sprintf('*%s*.sev', ExpName));
    for ch = 1:numel(SEVfiles{d})
        SEVfiles{d,ch} = wildcardsearch(fullfile(DataDir, SessionDates{d}), sprintf('*%s*_ch%d.sev', ExpName, ch));
        [rawdata, fs, header] = SEV2mat_singleCh(SEVfiles{d,ch});
        MatFilename{d,ch} = [SEVfiles{d,ch}(1:end-2), 'mat'];
        save(MatFilename{d,ch}, 'rawdata','fs','header');
    end
   	ProcDuration = toc(minortic);
    fprintf('Time elapsed = %s\n', datestr(ProcDuration/86400, 'HH:MM:SS'))
end

ProcDuration = toc(greattic);
fprintf('Time elapsed = %s\n', datestr(ProcDuration/86400, 'HH:MM:SS'))



%======================= RenameFiles

DataDir         = '/Volumes/RAWDATA/murphya/Physio/QNX/Dango/20170622/';
AllFiles        = wildcardsearch(DataDir, '*.dgz');

RemoveString    = '170622_';
ReplaceString   = '_20170622_';

for f = 1:numel(AllFiles)
    NewName = AllFiles{f};
    Indx    = strfind(NewName, RemoveString);
    if numel(Indx) ~= 1
        error('String to remove is either not unique or could not be found!');
    end
    if isempty(ReplaceString)
        NewName = NewName([1:Indx-1, (Indx+length(RemoveString)):end]);
    else
        NewName = [NewName(1:Indx-1), ReplaceString, NewName((Indx+length(RemoveString)):end)];
    end
    movefile(AllFiles{f},NewName);
end
%% Data parameters

TDTtankPath = '/data/NIMH_NIF/rawdata/koyanok/tdt/';
SavePath = '/data/NIMH_NIF/scratch/koyanok/';
% tanklist = {'Matcha160524';'Matcha160525';'Matcha160526';'Matcha160527';'Matcha160529';'Matcha160530';'Matcha160531';'Matcha160601';'Matcha160602'};
tanklist = {'Spice160606'; 'Spice160607'; 'Spice160609'};
channels = [1:64];
channel_group = 64;
channel_dead = [];
channel_dead_indx = ismember(channels,channel_dead);

filter_hp = 100;    % high-pass frequency
filter_lp = 5000;   % low-pass frequency
filter_n  = 4;      % filter order

% int16_scale_factor = 3276700;

%% Determine all Tanks in TDT path
if ~exist('tanklist','var')
    tanklist = [];
end
if isempty(tanklist)
    folderlist = dir(TDTtankPath);
    isd = [folderlist(:).isdir]; % only use folders
    tanklist = {folderlist(isd).name}';
    tanklist(ismember(tanklist,{'.','..'})) = []; % remove the . and .. from the list
    clear folderlist isd
end

%% Checking availability of statistics toolbox
disp('    checking statistics toolbox license...')
avail_license = 0;  % flag for statistics toolbox license
while ~avail_license
    msg = evalc('!licenses -stat');
    avail_license = str2double(msg(2));
    if avail_license>0
        disp('    statistics toolbox is available')
    else
        pause(5);
    end
end


%% Run the Preprocessing function for all paths 

greattic = tic;
for ii = 1:length(tanklist)                                                     % loop for tank
    NeuroTankPath = [TDTtankPath tanklist{ii} filesep];
    disp(['Beginning tank ' tanklist{ii}]);
    folderlist = dir(NeuroTankPath);
    isd = [folderlist(:).isdir];
    blocklist = {folderlist(isd).name}';
    blocklist(ismember(blocklist,{'.','..'})) = [];
    clear folderlist isd;
%     blocklist = {'S11_RestingState'};
return
    for jj = 1:length(blocklist)                                                % loop for block
        if strcmp(blocklist{jj},'TempBlk'); continue; end;
        BlockPath   = [NeuroTankPath blocklist{jj} filesep];
%         if ~(exist([SavePath  'proc_mat' filesep tanklist{ii}],'dir'))
%             mkdir([SavePath  'proc_mat' filesep tanklist{ii}]);SEV2mat_singleCh
%         end;
%         if ~(exist([SavePath  'proc_bin' filesep tanklist{ii}],'dir'))
%             mkdir([SavePath  'proc_bin' filesep tanklist{ii}]);
%         end;
%         if ~(exist([SavePath  'proc_bin' filesep tanklist{ii} filesep blocklist{jj}],'dir'))
%             mkdir([SavePath  'proc_bin' filesep tanklist{ii} filesep blocklist{jj}]);
%         end;
        if ~(exist([SavePath  'proc_sev' filesep tanklist{ii}],'dir'))
            mkdir([SavePath  'proc_sev' filesep tanklist{ii}]);
        end;
        if ~(exist([SavePath  'proc_sev' filesep tanklist{ii} filesep blocklist{jj}],'dir'))
            mkdir([SavePath  'proc_sev' filesep tanklist{ii} filesep blocklist{jj}]);
        end;
%         MatSaveFile = [SavePath  'proc_mat' filesep tanklist{ii} filesep blocklist{jj} '.mat'];
%         BinSavePath = [SavePath  'proc_bin' filesep tanklist{ii} filesep blocklist{jj} filesep];
        SevSavePath = [SavePath  'proc_sev' filesep tanklist{ii} filesep blocklist{jj} filesep];
        disp(['  Beginning block ' blocklist{jj}]);
        disp('    reading sev file...')
        for kk = 1:length(channels);                                            % loop for channel
            filelist = dir([BlockPath '*ch' num2str(channels(kk)) '.sev']);
            if ~isempty(filelist)
                disp(['      Ch ' num2str(channels(kk))]);
                SevFile = [BlockPath filelist(1).name];
                [rawdat, fs] = SEV2mat_singleCh(SevFile);                        % load data
                if ~exist('alldat','var')
                    alldat = nan(length(channels),length(rawdat));
                    fslist = nan(length(channels),1);
                end
                if isa(rawdat,'int16')
                    rawdat = single(rawdat);
                    dattype = 'int16';
%                     rawdat = rawdat./int16_scale_factor;
                elseif isa(rawdat,'single')
                    dattype = 'single';
                else
                    dattype = 'other';
                end
                
%                 if length(rawdat) > size(alldat,2)
%                     rawdat = rawdat(1:size(alldat,2));
%                 elseif length(rawdat) < size(alldat,2)
%                     rawdat = [rawdat nan([1 size(alldat,2)-length(rawdat)])];
%                 end
                
                alldat(kk,:) = rawdat;
                fslist(kk,:) = fs;
            else
                disp(['    Ch ' num2str(channels(kk)) ' does not exist']);
            end
            clear rawdat fs filelist SevFile;
        end
        disp('    all loaded');
        process_time = toc(greattic);
        process_hour = floor(process_time/3600);
        process_min  = floor(process_time/60)-process_hour*60;
        process_sec  = round(process_time - process_hour*3600 - process_min*60);
        disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
        disp(' ');
        
        %% Mean signal subtraction
        disp('    subtracting mean...');
        for kk = 1:ceil(length(channels)/channel_group)
            start_dat = channel_group*(kk-1)+1;
            end_dat   = channel_group*kk;
            if end_dat>length(channels)
                end_dat = length(channels);
            end
            disp(['      Ch ' num2str(channels(start_dat)) ' to ' num2str(channels(end_dat))]);
            channel_indx = start_dat:end_dat;
            channel_indx(channel_dead_indx(start_dat:end_dat)) = [];
            mean_subtract = nanmean(alldat(channel_indx,:),1);
            alldat(start_dat:end_dat,:) = alldat(start_dat:end_dat,:) - repmat(mean_subtract,[end_dat-start_dat+1,1]);
            clear start_dat end_dat mean_subtract;
        end
        disp('    subtracted');
        process_time = toc(greattic);
        process_hour = floor(process_time/3600);
        process_min  = floor(process_time/60)-process_hour*60;
        process_sec  = round(process_time - process_hour*3600 - process_min*60);
        disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
        disp(' ');
      

        %% Remove 1st principal component
        disp('    subtracting 1st principal component...');
        for kk = 1:ceil(length(channels)/channel_group)
            start_dat = channel_group*(kk-1)+1;
            end_dat   = channel_group*kk;
            if end_dat>length(channels)
                end_dat = length(channels);
            end
            disp(['      Ch ' num2str(channels(start_dat)) ' to ' num2str(channels(end_dat))]);
            channel_indx = start_dat:end_dat;
            channel_indx(channel_dead_indx(start_dat:end_dat)) = [];
            alldat(channel_indx,:) = (pcares(alldat(channel_indx,:)',1))';
            clear start_dat end_dat;
        end
        disp('    subtracted');
        process_time = toc(greattic);
        process_hour = floor(process_time/3600);
        process_min  = floor(process_time/60)-process_hour*60;
        process_sec  = round(process_time - process_hour*3600 - process_min*60);
        disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
        disp(' ');
          
        %% Filtering
        disp('    filtering...')
        fs = fslist(1);
        WnHP = filter_hp / (fs/2) ;
        WnLP = filter_lp / (fs/2) ;
        [bHP, aHP] = butter(filter_n, WnHP,'high');
        [bLP, aLP] = butter(filter_n, WnLP,'low');        
        tmp_filtered = filtfilt(bHP,aHP,alldat');
        alldat = filtfilt(bLP,aLP,tmp_filtered)';
        clear tmp_filtered;
        disp('    filtered.')
        process_time = toc(greattic);
        process_hour = floor(process_time/3600);
        process_min  = floor(process_time/60)-process_hour*60;
        process_sec  = round(process_time - process_hour*3600 - process_min*60);
        disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
        disp(' ');
        
%         %% Save binary data
%         disp('    savind binary data...');
%         for kk = 1:length(channels)
%             disp(['      Ch ' num2str(channels(kk))]);
%             BinSaveFile = [BinSavePath 'ch' num2str(channels(kk)) '.bin'];
%             savedat = alldat(kk,:);
%             savedat = int16(savedat.*1000.*(2^15));
%             fid=fopen(BinSaveFile,'w');
%             fwrite(fid,savedat,'int16');
%             fclose(fid);
%             clear savedat BinSaveFile;
%         end
%         disp('    saved.')
%         process_time = toc(greattic);
%         process_hour = floor(process_time/3600);
%         process_min  = floor(process_time/60)-process_hour*60;
%         process_sec  = round(process_time - process_hour*3600 - process_min*60);
%         disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
%         disp(' ');
        
        %% Save sev data
        disp('    savind sev data...');
        for kk = 1:length(channels)
            disp(['      Ch ' num2str(channels(kk))]);
            SevSaveFile = [SevSavePath 'ch' num2str(channels(kk)) '.sev'];
            fid=fopen(SevSaveFile,'w');
            switch dattype
                case 'int16'
                    fwrite(fid,alldat(kk,:),'int16');
                case 'single'
                    fwrite(fid,alldat(kk,:),'single');
                otherwise
                    fwrite(fid,alldat(kk,:));
            end
            fclose(fid);
            clear SevSaveFile;
        end
        disp('    saved.')
        process_time = toc(greattic);
        process_hour = floor(process_time/3600);
        process_min  = floor(process_time/60)-process_hour*60;
        process_sec  = round(process_time - process_hour*3600 - process_min*60);
        disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
        disp(' ');
        
%         %% Save mat data
%         disp('    savind mat data...')
%         for kk = 1:length(channels)
%             disp(['      Ch ' num2str(channels(kk))]);
%             if channels(kk)<10
%                 channel_string = ['00' num2str(channels(kk))];
%             elseif channels(kk)<100
%                 channel_string = ['0'  num2str(channels(kk))];
%             else
%                 channel_string =       num2str(channels(kk));
%             end
%             var_name = ['MUA' channel_string];
%             eval([var_name '= alldat(kk,:);']);
%             if exist(MatSaveFile,'file')
%                 save(MatSaveFile,var_name,'-append');
%             else
%                 save(MatSaveFile,var_name,'-v7.3');
%             end
%             clear(var_name,'var_name');
%         end
%         disp('    saved.')

        clear alldat;
    end
    clear blocklist;
    disp(['Finished ' tanklist{ii} ' convert']);
    process_time = toc(greattic);
    process_hour = floor(process_time/3600);
    process_min  = floor(process_time/60)-process_hour*60;
    process_sec  = round(process_time - process_hour*3600 - process_min*60);
    disp([num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
    disp(' ');
end

disp('Finished convertitng');
process_time = toc(greattic);
process_hour = floor(process_time/3600);
process_min  = floor(process_time/60)-process_hour*60;
process_sec  = round(process_time - process_hour*3600 - process_min*60);
disp(['    ' num2str(process_hour) 'h ' num2str(num2str(process_min)) 'min ' num2str(process_sec) 'sec']);
disp(' ');
clear tanklist;
% exit
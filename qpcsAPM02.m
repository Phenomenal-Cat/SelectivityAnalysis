function qpcsAPM02(DebugOn)
%% ============================= qpcsAPM02.m ==============================
% Sets up socket connection between QNX machine and stimulus display PC. 
%
% EXPERIMENTS:
%   apm_attention:      Posner-style cued attention paradigm
%   apm_binorivOKN: 	Binocular rivalry generating optokinetic nystagmus (OKN)
%   apm_binorivPA:      Sham binocular rivalry with physical alternation between eyes
%   apm_binorivSLOF:    Binocular rivalry generating short-latency occular following (SLOF)
%   apm_calib:          9-point calibration for eyetracking
%   apm_CFS:            Continuous flash suppression
%   apm_CFSBehav:       behavioural CFS experiment
%   apm_DelaySac:       Delayed saccade experiment
%   apm_movie:          Movie presentation
%   apm_motion:         Presents short animations (optionally looped) as trials
%   apm_picTuning:      Presents static images from a specified directory containing .pngs
%   apm_picTuningStereo:Presents static images at varying positions in depth
%   apm_Retinotopy:     Presents rotating wedge and expanding or contracting ring stimuli
%   apm_RFMapping:      Receptive field mapping using motion patches
%   apm_SFM:        	Bistable structure-from-motion stimulus
%   apm_SFMPA:          Stereoscopically disambiguated SFM stimulus
%   apm_StereoFaces:    presents stereoscopic images of faces
%   apm_StereoTest:     Random dot stereogram version of 9-point calibration
%   apm_Tonotopy:       Presents audio at a range of frequencies
%   
% KEYBOARD INPUTS:
%   'm':        Toggle menu on/off (icons appear top left screen)
%   's':        Play start tone (if audio is on)
%   'p':        deliver penalty (screen flash)
%
% REQUIREMENTS:
%   APMsubfunctions:    Matlab path set to include subfunction directory 
%   ImageDir:           Image directory containing any image format
%   MovieDir:           Movie directory containing .avi/.mp4/.mov files
%
% REVISIONS:
%   14/11/2012 - Adapted from qpcsBER.m for use with csd_picTuning.c on QNX.
%   14/01/2013 - movie function added
%   01/02/2013 - intermittent presentation, frequency tagging and
%                short-latency occular following options added to binocular rivalry.
%   06/02/2013 - binorivPA, RFMapping & attention experiments added
%   14/02/2013 - SFM and SFMPA experiments added
%   20/02/2013 - Binocular rivalry with motion probe added. Audio added.
%   25/03/2013 - Photodiode target, stereo vision test and delayed saccade experiment added.
%   02/05/2013 - Yoked transition options added to binorivPA
%   14/05/2013 - Continuous flash suppression experiment added
%   26/06/2013 - Retinotopic localizer experiment added + Matlab keyboard control
%   14/01/2014 - Behavioural dCFS experiment added
%   31/01/2013 - MovieClips experiment added
%   22/10/2015 - Scene tranistion markers added via photodiode signal in movie experiment
%   03/12/2015 - Added picTuningStereo experiment
%     ___  ______  __   __
%    /   ||  __  \|  \ |  \    APM SUBFUNCTIONS
%   / /| || |__/ /|   \|   \   Aidan P. Murphy - murphyap@mail.nih.gov
%  / __  ||  ___/ | |\   |\ \  Section on Cognitive Neurophysiology and Imaging
% /_/  |_||_|     |_| \__| \_\ NIMH
%==========================================================================

% Declare persistent variables (available to nested functions)
persistent Display
persistent currentSystem;
persistent currentSystemHandle;
persistent reply;
persistent con;
persistent ImageDir;
persistent MotionDir;
persistent Audio;
persistent Photodiode;
persistent Penalty;
persistent GammaCorrect;
persistent Key;
persistent Debug;
persistent Sleep;
persistent Exit;
persistent Button;
persistent EPI;
persistent Metronome;

if nargin < 1
    Debug.On = 0;                                       	% Set to 1 for half-size screen window and printed info
else
	Debug.On = DebugOn;
end

try
    addpath(genpath('../../APMSubfunctions'));            	% Add APMSubfunctions to MATLAB path
    
 	%================== SET DISPLAY SETTINGS
    Display = DisplaySettings(1);                           % Get display settings 
%     if strcmpi(Display.Settings.Bits,'win32')
%         h = errordlg('You are currently running a 32-bit version of Matlab! Video playback will not work!');
%         uiwait(h);
%         return;
%     end
    Display.MultiFlip       = 1;                                    % Screen('Flip'... to all windows simultaneously
    Display.DontSync        = 1;                                    % If used in the flip command, Matlab does not wait for next refresh before continuing
    Display.Background      = [0.5 0.5 0.5]*255;                    % Set default background color to mid-grey
    Display.FlashBackground = [255 0 0];                            % Set color of 'flash' screen to red
    Penalty.On           	= 0;
    Penalty.NoFlashes       = 6;
    Display.Mirror          = 0;                                	% QNX will calculate mirroring for monocular presentations
    UsePsychImaging         = 1;  
    GammaCorrect.On         = 1;                                    % Load gamma table?
  	if GammaCorrect.On == 1
        Display.CLUT = 'NIH_Setup3_ASUS_V27_reduced.mat';
        Display.OriginalGamma = Screen('ReadNormalizedGammaTable', Display.ScreenID);
        try
            load(Display.CLUT);                            	% Load gamma lookup tables
            try
            Screen('LoadNormalizedGammaTable', Display.ScreenID, inverseCLUT{2});   % Apply gamma table
            catch
                Screen('LoadNormalizedGammaTable', Display.ScreenID, inverseCLUT{2});   % Apply gamma table
            end
%             for Eye = 1:2
%                 load(Display.CLUT{Eye});                                                                          % Load gamma lookup tables
%                 Screen('LoadNormalizedGammaTable', Display.ScreenID, inverseCLUT{Eye},[],Display.ScreenID{Eye});	% Apply gamma table
%             end
        catch
            fprintf('\nERROR: Attempt to upload CLUT %s to graphics adapter failed!\n', Display.CLUT);
            GammaCorrect.On = 0;
        end
    end
    
  	%=================== SET PHOTODIODE TARGET
    Photodiode.On = 1;                                      % Photodiode target defaults to 'on'
 	Photodiode.Diameter = 0.018*Display.Pixels_per_m(1); 	% size of target (pixels)
    Photodiode.Position = 1;                                % bottom left corner (2 = bottom right corner)
    if Photodiode.On == 0
       	Photodiode.OnColour =  Display.Background;
        Photodiode.OffColour = Display.Background;             
    elseif Photodiode.On == 1
        Photodiode.OffColour = [128 128 128];
        Photodiode.OnColour = [0 0 0];   
    end
    Photodiode.Rect{1} = [Display.Rect(1),Display.Rect(4)-Photodiode.Diameter, Display.Rect(1)+Photodiode.Diameter, Display.Rect(4)];   % Bottom LEFT for RIGHT eye
    Photodiode.Rect{2} = [Display.Rect(3)-Photodiode.Diameter,Display.Rect(4)-Photodiode.Diameter, Display.Rect(3), Display.Rect(4)];   % Bottom RIGHT for LEFT eye
    
	Display.Rect = Screen('rect',Display.ScreenID);         % Get dual screen (horizontal span) resolution 
    if Debug.On == 1, Display.Rect = Display.Rect/2; end  	% For debugging, decrease onscreen window size
    
    %================== DEFINE MATLAB KEYBOARD INPUTS
    KbName('UnifyKeyNames');
    Key.Exit        = KbName('Escape');
    Key.Penalize    = KbName('p');
    Key.LightsOut   = KbName('b');
    Key.LightsOn    = KbName('g');
    Key.StartTone   = KbName('s');
    Key.Audio       = KbName('a');
    Key.Diode       = KbName('d');
    Key.Training    = KbName('t');
    Key.EPI         = KbName('e');
    Key.Menu        = KbName('m');
    Key.MinInterval = 0.2;
    Key.LastPress   = GetSecs;
    
    %================= SET PATH FOR IMAGES
    EPI.On = 0;
    rootDir = fileparts(mfilename('fullpath'));
%     ImageDir = fullfile(rootDir,'ALL IMAGES');
%     ImageDir = 'P:\murphya\Stimuli\Images';
%     ImageDir = 'P:\dengc\stimuli\Images\Processed';
%     ImageDir = 'P:\murphya\Stimuli\Processed';
    ImageDir = 'P:\koyanok\Stimuli\Categories_Stim10k_small';
    MotionDir = 'P:\murphya\Stimuli\MotionClips';

    %=========================== OPEN PTB WINDOW ==============================
    HideCursor;                                         % Hide mouse pointer
    if usejava('jvm')==1
        ListenChar(2);                                 	% Ignore keyboard input to command line 
    end
    rand('seed',GetSecs);                               % Seed random number generator from current time
    warning off all;                                    % Turn off Matlab warning output
    Screen('Preference', 'VisualDebugLevel', 1);        % Make initial screen black instead of white
    if UsePsychImaging == 1
        PsychImaging('PrepareConfiguration');                                       % Setup psychimaging pipeline
%         PsychImaging('AddTask', 'AllViews', 'RestrictProcessing', CenterRect([0 0 512 512], Screen('Rect', Display.ScreenID)));
        [Display.win, Display.Rect] = PsychImaging('OpenWindow', Display.ScreenID, Display.Background(1), Display.Rect, [], [], Display.Stereomode);
%         Screen('OpenWindow', Display.SlaveScreenID, Display.Background(1), [], [], [], Display.Stereomode);
    else
        [Display.win, Display.Rect] = Screen('OpenWindow', Display.ScreenID, Display.Background(1),Display.Rect,[],[], Display.Stereomode);
    end
    Screen('Preference', 'SkipSyncTests', 0);                                       % Do not skip sync tests
    priorityLevel = MaxPriority(Display.win);                                       % Set high priority
    Priority(priorityLevel);                                                        
    Display.FlipInt = Screen('GetFlipInterval', Display.win);                       % Get flip interval
    Screen('BlendFunction', Display.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    	% Enable alpha channel   
    Screen('Preference', 'DefaultVideocaptureEngine', 3);                           % 0 = Quicktime; 3 = GStreamer
    Screen('Flip',Display.win);       

	%======================== INITIALIZE AUDIO SETTINGS =====================
    Audio.On = 0;                                       % Play audio feedback tones/ movie sound?
    [Audio.Beep, Audio.Noise] = AuditoryFeedback;       % Generate tones
    Audio.Beep = Audio.Beep/3;                          % Reduce volume
    Audio.Noise = Audio.Noise/3;
    if ~IsWin
        Audio.PsychSound = 1;
        Speak('Audio initialized.');                    % Inform user that audio is ready
    elseif IsWin
        Audio.PsychSound = 0;
        Audio.a = actxserver('SAPI.SpVoice.1');
%         Audio.a.Speak('Audio initialized.');
        Audio.a.Speak('Initializing experiments on stim S 4.');
    end
    StartTone = [Audio.Beep(1,:), zeros(1,numel(Audio.Beep(1,:))),Audio.Beep(1,:), zeros(1,numel(Audio.Beep(1,:))), Audio.Beep(1,:)];
    Snd('Play',StartTone);  
    if ~ismac
        Audio.Error = Audio.Beep(2,:);
        Audio.Penalty = Audio.Noise;
    end
    
    %========================== LOAD ICONS AND LOGO =========================
	Button.IconNames = {'ButtonOn.png','ButtonOff.png','ButtonHighlight.png','SpeakerOn.png','SpeakerOff.png','Penalty.png','Photodiode.png','GammaCorrect.png','Sleep.png','Debug.png','EPI.png','Metronome.png','Exit.png','NIMHlogo.png','BadMonkey.png'};%,'APMsubfunctions_logo.png'};
    Button.IconFunctions = {'Audio','Penalty','Photodiode','GammaCorrect','Sleep','Debug','EPI','Metronome','Exit'};
    Button.On = zeros(1,numel(Button.IconFunctions));
    Sleep.On = 1;
    Exit.On = 0;
    Metronome.On = 0;
    
  	NIMHlogoIcon = [];
    c = []; Alpha = [];
    SpeakerOnIcon = [];
    SpeakerOffIcon = [];
    PenaltyIcon = [];
    PhotodiodeIcon = [];
    APMsubfunctions_logoIcon = [];
    BadMonkeyIcon = [];
    GammaCorrectIcon = [];
    SleepIcon = [];
    ButtonOnIcon = [];
    ButtonOffIcon = [];
    ButtonHighlightIcon = [];
    DebugIcon = [];
    EPIIcon = [];
    MetronomeIcon = [];
    ExitIcon = [];
    
    for i = 1:numel(Button.IconNames)
        Button.IconName{i} = strcat(Button.IconNames{i}(1:end-4),'Icon');
        eval(sprintf('[%s,c,Alpha] = imread(''%s'');', Button.IconName{i},Button.IconNames{i}));              % Load icon
        if ~isempty(Alpha) 
            eval(sprintf('%s(:,:,4) = Alpha;',Button.IconName{i}));
        else
            eval(sprintf('%s(:,:,4) = repmat(255, size(%s,1),size(%s,2))-double(%s(:,:,3));',Button.IconName{i},Button.IconName{i},Button.IconName{i},Button.IconName{i}));
        end 
    end
    Button.Icon(1) = Screen('MakeTexture',Display.win, ButtonOffIcon);
    Button.Icon(2) = Screen('MakeTexture',Display.win, ButtonOnIcon);
    Button.Highlight(1) = Screen('MakeTexture',Display.win, zeros(size(ButtonOnIcon)));
    Button.Highlight(2) = Screen('MakeTexture',Display.win, ButtonHighlightIcon);
    Audio.SpeakerIcon(1) = Screen('MakeTexture',Display.win, SpeakerOffIcon);
    Audio.SpeakerIcon(2) = Screen('MakeTexture',Display.win, SpeakerOnIcon);
    Penalty.BadMonkey = Screen('MakeTexture',Display.win, BadMonkeyIcon);
 	APMTex = Screen('MakeTexture',Display.win, NIMHlogoIcon);%APMsubfunctions_logoIcon);
    LogoRect = [0 0 size(NIMHlogoIcon,2),size(NIMHlogoIcon,1)]/2;
%     LogoRect = [0 0 size(APMsubfunctions_logoIcon,2),size(APMsubfunctions_logoIcon,1)];
 	LogoDestRect = CenterRect(LogoRect, Display.Rect);
    LogoOffset = -0.3*Display.Pixels_per_deg(1);                                                % Set binocular disparity of logo
    for e = 2:-1:1
        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, e-1);
        Offset = (round(e-1.5)*LogoOffset);
        Screen('DrawTexture',Display.win,APMTex,[],LogoDestRect+[Offset,0,Offset,0]);           % Display main logo
    end
    
    Button.DestRect(1,:) = [40 20 100 80];
    for b = 2:numel(Button.IconFunctions)
        eval(sprintf('%s.Icon = Screen(''MakeTexture'',Display.win, %s);', Button.IconFunctions{b},Button.IconName{b+4}));  % Make textures
        Button.DestRect(b,:) = Button.DestRect(b-1,[3 2 3 4])+[10 0 Button.DestRect(1,3)-Button.DestRect(1,1)+10, 0];       % Calculate rect
    end
    for b = 1:numel(Button.IconFunctions)
        eval(sprintf('Button.On(b) = %s.On;',Button.IconFunctions{b}));                                                     % Check function status
        Screen('DrawTexture',Display.win,Button.Icon(Button.On(b)+1),[],Button.DestRect(b,:));                              % Draw appropriate button background
        Screen('DrawTexture',Display.win,Button.Highlight(1),[],Button.DestRect(b,:));                                      % Draw appropriate button highlight
    end
    Screen('DrawTexture',Display.win,Audio.SpeakerIcon(Audio.On+1),[],Button.DestRect(1,:));                                % Display default audio settings icon
    for b = 2:numel(Button.IconFunctions)
        eval(sprintf('Screen(''DrawTexture'',Display.win,%s.Icon,[],Button.DestRect(b,:));', Button.IconFunctions{b}));     % Display default settings for all functions
    end
    
    
    %================ Set up socket and wait for connection and input =========
    if ~isfield(Display.Settings, 'IPaddress')
        Display.Settings.IPaddress = '(IP unknown)';
    end
    Instruction=sprintf(['Execute ESS interface control...\n\n\nSystem executable: essapm\n\n',...
        'Remote system IP: %s\n\nRemote system name: %s'], Display.Settings.IPaddress, Display.Settings.CompName);
    Screen('TextFont', Display.win, 'Arial');
    for e = 1:2
        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, e-1);                     % Draw to screen
        DrawFormattedText(Display.win, Instruction, 40, Display.Rect(4)-200, [0 0 0], []);      % Draw text
    end
    Screen('Flip',Display.win);

    % readTimeOut = 0.001;
    readTimeOut = 2;                                        % how long to wait for the socket info
    sockcon=pnet('tcpsocket',4610);
    pnet(sockcon,'setreadtimeout',readTimeOut);
    con=pnet(sockcon,'tcplisten');
    if con~=-1
        pnet(con,'setreadtimeout',readTimeOut);
    end


%% ==================== WAIT FOR COMMAND FROM QNX =========================
% run blocking readline and response commands
    while(1)
        if(pnet(con,'status')==0||con==-1)
            con=pnet(sockcon,'tcplisten'); 
            pnet(con,'setreadtimeout',2);
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
            if keyIsDown && keyCode(Key.Exit)               % if the esc key is pressed abort
                break;
            end
            continue;
        end
        commandIn = pnet(con,'readline');
        if strcmp(commandIn,'')
        	%================= Check for commands from STIM PC keyboard
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
            if keyIsDown && keyCode(Key.Exit)
                break;
            end
            if keyIsDown && secs > Key.LastPress+Key.MinInterval
                Key.LastPress = secs;
                if keyCode(Key.Penalize)
                    R_punishStim;
                elseif keyCode(Key.StartTone)                               % Play tone to alert monkey to experiment start
                    PlayAudio(Audio, StartTone);
                elseif keyCode(Key.Menu)  
                    ShowMenu;
            	elseif keyCode(Key.EPI)
                    EPI.On = ~EPI.On;
                    PlayEPI(EPI.On);
                end
            end
            continue;
        end
        if(commandIn=='q')
            break
        end
        if(pnet(con,'status')==0)
            con=pnet(sockcon,'tcplisten');
            continue;         
        end
        try
            reply = '\n';
            parseCommand(commandIn); % determines what function/subfuction to run and then executes it. 
        catch
    %             psychrethrow(psychlasterror);
            rethrow(lasterror)
        end
        pnet(con,'printf', reply);                              % indicate completion by sending the reply string back to the QNX computer.
    end

    %======================= FINISH AND CLOSE CONNECTION ======================
    if isfield(Display,'OriginalGamma')
        Screen('LoadNormalizedGammaTable', Display.ScreenID, Display.OriginalGamma);
    end
    Screen('Close');                                            % Close all open PTB textures
    Screen('CloseAll');                                         % Close PTB window
    ListenChar(0);                                              % Restore command line keyboard output
    ShowCursor;                                                 % Show cursor
    pnet('closeall');                                           % Close socket connection
catch
    if isfield(Display,'OriginalGamma')
        Screen('LoadNormalizedGammaTable', Display.ScreenID, Display.OriginalGamma);
    end
  	Screen('Close');                                            % Close all open PTB textures
    Screen('CloseAll');                                         % Close PTB window
    if usejava('jvm')==1
        ListenChar(0);                                              % Restore command line keyboard output
    end
    ShowCursor;                                                 % Show cursor
    pnet('closeall');                                           % Close socket connection
    Error = lasterror                                          % Get last error message
%     opentoline(mfilename,  Error.stack(1).line);
    fprintf('\n\nERROR: in  <a href="matlab: opentoline(which(''%s.m''),%d)">%s.m line %d</a> : %s\n\n', mfilename, Error.stack(1).line, mfilename, Error.stack(1).line, Error.message);
    rethrow(lasterror);                                         
end


%% ========================= parseCommand(command)=========================
% Reads in commands sent through the socket and creates function calls out of them.
% parseCommand is called by by qpcsAPM at the end of the while loop.
% command should be a single text string
%   - the first "word" in the string should be a function to call
%	- each "word" after that will be read in as an arguement to the function
%   - each arguement should be separated by a space
%==========================================================================
    function parseCommand(command) 
        Time = datestr(now);
        fprintf('%s :Parsing command string: %s\n', Time(end-7:end), command);
        commandArray = textscan(command,'%s');                              % Convert command string to a cell array broken up by spaces
        switch(sprintf(cell2mat(commandArray{1}(1))))                       % If the first word is a command...
            case 'ping'
                reply = sprintf('pong %s\n',sprintf(cell2mat(commandArray{1}(2))));
                return
            case 'setsystem'                                                % Initialize the appropriate functions for the experiment to run
                currentSystem = sprintf(cell2mat(commandArray{1}(2)));      % Convert experiment name from cell to string
                try
                    currentSystemHandle = eval(strcat(currentSystem,'()'));	% Call experiment function and return a cell with the list of subfunctions
                catch
                    psychrethrow(psychlasterror); 
                end
                return;
            case 'q'                                                        % Close the session
                pnet('closeall');
                return;
        end
        %******************************************************************
        % The rest of the commands reference the system (persistent currentSystem)
        %******************************************************************
        getComd = @(x) sprintf(cell2mat(commandArray{1}(x)));       % returns command at index x
        comdEval = '';                                              % initializes comdEval to be filled with arguments for the function to be called
        for i =2:size(commandArray{1},1);
            if i==2 
                comdEval=str2num(getComd(2));                       % resets comdEval with first sent argument
                continue;
            end
           % comdEval = strcat([comdEval,' ',getComd(i)]); % actually fills the variable with the sent arguments 
            comdEval = [comdEval str2num(getComd(i))];
        end
        for i=1:size(currentSystemHandle,1)                         % loops through all the subfunctions of the current experiment
            if regexp(func2str(currentSystemHandle{i}),getComd(1),'once') > size(currentSystem,2)
                if size(comdEval)==0
                    currentSystemHandle{i}();                       % calls current subfunction of experiment with no arguments
                else
                    currentSystemHandle{i}(comdEval);               % calls current subfunction of experiment with arguments
                end
                return;
            end
        end
        clear getComd comdEval 
    end


%% ========================== StereoCalib =================================
% A 9-point calibration routine for manual calibration of eye tracking
% devices.
%
%==========================================================================

    function stereocalibHandle = apm_calib()

        stereocalibHandle = ...
            {...
                @R_punishStim;...
                @R_fixOn;...
                @R_setFixColor;...
                @R_remoteTargetOn;...
                @R_getStmParamByName;...
                @R_getStmParamName;...
                @R_fixOff;...
                @R_getStmNParams;...
                @R_getStmParam;...
                @R_getBlock;...
                @R_nextTrial;...
                @R_clearscreen;...
                @R_reset;...
                @R_clearstim;...
                @R_getFixPosX;...
                @R_getFixPosY;...
                @R_getFixEye;...
                @newStmParam;...
                @stereocalibReset;...
            };

        persistent stereocalib;
        persistent Fix;
        
        Screen('Flip', Display.win);    % Calibration is default experiment, so clear screen once socket connection is made
        stereocalib.eye = 1;            
        Fix.Type = 0;                   
        
        %======================== PRESENT FIXATION MARKER =================
        function R_fixOn (varargin)
         	params=varargin{1,1};
            Fix.Offset = [params(1) -params(2)]*Display.Pixels_per_deg(1);
            Fix.Eye = params(3);
            stereocalib.eye = params(3);
            Fix.Size = params(4)*Display.Pixels_per_deg(1);
            stereocalib.EccDegX = params(5);
            stereocalib.EccDegY = params(6);
            Fix = fixOn(Fix);
%             if isfield(Display,'DestRect')
%                 Screen('SelectStereoDrawBuffer', Display.win, abs(floor(Fix.Eye-0.5)));
%                 array2Flip = Screen('GetImage', Display.win, Display.Rect, 'backBuffer'); 
%                 TextTexture = Screen('MakeTexture', Display.win, array2Flip);
%                 Screen('SelectStereoDrawBuffer', Display.win, abs(floor(Fix.Eye-0.5)));
%                 Screen('DrawTexture',Display.win,TextTexture, [], Display.DestRect);
%             end
        	[VBL FixOn] = Screen('Flip', Display.win);
        end

        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)];
        end

        function R_remoteTargetOn ()

        end

        function R_getStmParamByName (varargin)
            params=varargin{1,1};
            if strcmp(params,'xpos')
                R_getFixPosX();
            elseif strcmp(params,'ypos')
                R_getFixPosY();
            elseif strcmp(params,'eye')
                R_getFixEye();
            else
                reply = 'ERROR\n';
            end
        end

        function R_getStmParamName ()
            reply = 'xpos\n';
        end

        function R_fixOff()
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background);
            end 
            [VBL FixOff] = Screen('Flip', Display.win);   
        end

        function R_getStmNParams ()
        	reply = '3\n';
        end

        function R_getStmParam ()
        	reply = '0.0\n';
        end

        function R_getBlock(n)
            stereocalib.trial_num = 0;
            A = perms([-1,0,1]);
            A = [A(:,1:2); -1,-1;0,0;1,1];
            A = [A, zeros(9,1); A, ones(9,1)];
            A = [A,randperm(18)'];
            A = sortrows(A,4);
            stereocalib.block = A(:,1:3);
            reply = sprintf('%d\n',size(stereocalib.block,1));
        end

        function R_nextTrial()
             if(stereocalib.trial_num < size(stereocalib.block,1)-1)
                stereocalib.trial_num = stereocalib.trial_num + 1;
             else
                 stereocalib.trial_num = -1;
                 reply = sprintf('-1\n');
                 return;
             end
             stereocalib.xpos = stereocalib.block(stereocalib.trial_num,1);
             stereocalib.ypos = stereocalib.block(stereocalib.trial_num,2);
             reply = sprintf('%d %d %d %d\n',stereocalib.trial_num,stereocalib.eye,stereocalib.xpos,stereocalib.ypos);
        end

        function R_clearscreen()         

        end

        function R_reset (varargin)
            R_fixOff;
        end

        function R_clearstim ()

        end

        function R_getFixPosX(varargin)
            stereocalib.EccDegX = varargin{1,1};
            reply = sprintf('%d\n',stereocalib.xpos*stereocalib.EccDegX); 
        end

        function R_getFixPosY(varargin)
            stereocalib.EccDegY = varargin{1,1};
            reply = sprintf('%d\n',stereocalib.ypos*stereocalib.EccDegY);           
        end

        function R_getFixEye()
            reply = sprintf('%d\n',stereocalib.eye);
        end

        function newStmParam (varargin)

        end

        function stereocalibReset()

        end
    end

%% =========================== apm_Tonotopy ===============================

     function ToneHandle = apm_Tonotopy()

            ToneHandle = ...
                {...
                    @R_nextTrial;...
                    @R_getBlock;...
                    @R_loadTones;...
                    @R_toneOn;...
                    @R_toneOff;...
                    @R_fixOn;...
                    @R_fixOff;...
                };

            persistent Event;
            persistent Trial_num;
            persistent Block;
            persistent Tone;


         %========================= LOAD TONES  
         function R_loadTones(varargin)
         	params = varargin{1,1};
            Tone.NoFreq = params(1);                % How many different frequencies
            Tone.MinFreq = params(2);               % Minimum frequency (Hz)
            Tone.MaxFreq = params(3);               % Maximum frequency (Hz)
            Tone.Duration = params(4);              % Tone duration (ms)
            Tone.ScaleType = params(5);             % 1 = use linear scale, 2 = use log scale
            Tone.Type = 1;                          % Tone type (1 = sine wave; 2 = band filtered white noise)
            Snd('Open');                            
            Tone.SampleRate = Snd('DefaultRate');	% Sample rate (Hz)
            if Tone.ScaleType == 1
                Tone.AllFreq = linspace(Tone.MinFreq,Tone.MaxFreq, Tone.NoFreq);
            elseif Tone.ScaleType == 2
                 Tone.AllFreq = logspace(2,4, Tone.NoFreq);
            end
            for f = 1:Tone.NoFreq
                LoadingText = sprintf('Loading %d ms, %d Hz tone (tone %d of %d)...', Tone.Duration, Tone.AllFreq(f), f, Tone.NoFreq);
                DrawFormattedText(Display.win, LoadingText, 20, Display.Rect(4)-40, [0 0 0], []);   % Draw text
            	Screen('Flip',Display.win);
                Tone.SineWav{f} = repmat(MakeBeep(Tone.AllFreq(f), Tone.Duration/1000), [2,1]);
            end  
            Screen('Flip',Display.win);
         end
            
         %========================= PLAY TONE   
         function R_toneOn(varargin)
         	params = varargin{1,1};
            FreqIndx = params(1)+1;                 % Which frequency tone to play?
            EarIndx = randi(3); %params(2);      	% Which side? (1 = left, 2 = right, 3 = stereo)
            Volume = params(3);                     % Volume (0-1)
            
            if EarIndx == 3
                EarIndx = [1,2];
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
            end
            [VBL ToneOn] = Screen('Flip', Display.win);
            Snd('Play', Tone.SineWav{FreqIndx}(EarIndx,:)*Volume);
         end

          %========================= STOP TONE   
         function R_toneOff(varargin)
%             Snd('Quiet');
        	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ToneOff] = Screen('Flip', Display.win);
         end

         function R_fixOn(varargin)
             
         end
         function R_fixOff
             
         end
         
     end

%% =========================== apm_picTuning ==============================
% Presents static stimuli from a specified directory containing .jpg images
%
    function ImageTuningHandle = apm_picTuning()

        ImageTuningHandle = ...
            {...
                @R_setFixColor;...
                @R_nextTrial;...
                @R_getBlock;...
                @R_getFixPosX;...
                @R_getFixPosY;...
                @R_fixOn;...
                @R_loadPic;...
                @R_preloadPics;...
                @R_picOn;...
                @R_picOff;...
                @R_fixOff;...
                @R_punishStim;... 
            };
        persistent ImageTexture;
        persistent MaskTex;
        persistent Event;
        persistent Trial_num;
        persistent Block;
        persistent Pic;
        persistent Mask;
        persistent Fix;

        Screen('BlendFunction', Display.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        
        %====================== SET STIMULUS DEFAULT PARAMETERS ==============
        Fix.Type = 1;                                                   % 0 = dot; 1 = cross; 2 = square; 3 = crosshairs
        Fix.Eye = 2;                                                    % Present fixation binocularly
     	Fix.Drawn = 0;
        Event = [0 0 0];
        if isfield(Fix,'Texture')
           Fix = rmfield(Fix,'Texture');
        end
        
     	%====================== LOAD IMAGE FILES AS TEXTURES ==================  
        function R_preloadPics(varargin)
            params = varargin{1,1};
            Pic.ImageIndx   = params(1):params(2);
            NoImages        = numel(Pic.ImageIndx);                                           
            ImageType       = params(3);                                    % Image category (e.g. snakes, faces, etc)
            ImageRadius     = 5*Display.Pixels_per_deg(1);              	% Default radius of image size BEFORE drawing based on size input                 
            TextureWindow   = [2*ImageRadius, 2*ImageRadius];            	% Set size of the stimulus before mask is applied (degrees)
        
            if exist('ImageTexture','var')                                  
                Screen('Close', ImageTexture);                           	% close any previously loaded image textures
            end
            ImageFormat = '.png';
            reverseStr = '';
            for n = 1:NoImages
                if Debug.On == 1
                    LoadingText = sprintf('Loading image %d of %d... (#%d - %d)',n,NoImages,Pic.ImageIndx(1),Pic.ImageIndx(end));
                    fprintf([reverseStr, LoadingText]);
                    reverseStr = repmat(sprintf('\b'), 1, length(LoadingText));
                    DrawFormattedText(Display.win, LoadingText, 20, Display.Rect(4)-40, [0 0 0], []);% Draw text
                    Screen('Flip',Display.win);
                end
                if any(strfind(ImageDir, 'koyano')) == 1
                    ImageNames{n} = sprintf('cat%03d%s', Pic.ImageIndx(n), ImageFormat);
                else
                    ImageNames{n} = [num2str(Pic.ImageIndx(n)),ImageFormat];                  	% Construct filename of next image
                end
                [Img ImgCmap] = imread(fullfile(ImageDir, ImageNames{n}));                  % Load image file
                ImgDim = size(Img);                                                         % Get image dimensions
                if min(ImgDim([1,2])) ~= 2*ImageRadius                                    	% if smallest image dimesnion is not requested size...
                    scale = 2*ImageRadius/min(ImgDim([1 2]));                               % Resize image so smallest dimension fits
                    Img = imresize(Img, scale);
                    ImgDim = size(Img);                                                     % Get new Img image dimensions
                end         
                if max(ImgDim([1,2])) > 2*ImageRadius                                   	% If largest image dimension is too large...
                    Crop = (max(ImgDim([1,2]))-(2*ImageRadius))/2;
                    if find(ImgDim==max(ImgDim([1,2])))==1
                        Img(end-Crop:end,:,:) = [];
                        Img(1:Crop,:,:) = [];
                    elseif find(ImgDim==max(ImgDim([1,2])))==2
                        Img(:,end-Crop:end,:) = [];
                        Img(:,1:Crop,:) = [];
                    end
                end
                Pic.ImageRect = [0,0,size(Img,2),size(Img,1)];
                ImageTexture(n) = Screen('MakeTexture', Display.win, double(Img));          % Convert image to PTB texture
            end
            fprintf('\nImages loaded sucessfully!\n\n');
            Screen('Fillrect', Display.win, Display.Background, Display.Rect);
            Screen('Flip', Display.win, [], 1);
            
            %========================= CREATE ALPHA CHANNEL MASK
            Mask.Dim = TextureWindow+4;                                                 % Mask will be 2 pixels larger than the stimulus window
            Mask.ApRadius = (min(Mask.Dim)/2)-8;
            Mask.Colour = Display.Background(1);
            Mask.Edge = 2;                                                              % Cosine mask edge
            Mask.Taper = 0.2;                                                           % Cosine edge tapers off over 20% of aperture radius
            MaskTex = GenerateAlphaMask(Mask, Display);
            Mask.Rect = [0 0 Mask.Dim];                                                 % Mask size
            Mask.DestRect = CenterRect(Mask.Rect, Display.Rect);                     	% Destination for mask is centred in screen window
        end

        %======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            if Fix.Drawn == 0
                Fix = fixOn(Fix);
                Fix.Drawn = 1;
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
              	Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
            end
        	[VBL FixOn] = Screen('Flip', Display.win, [], 1);
        	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            When = FixOn +(1/Display.RefreshRate)+0.001;
         	[VBL FixOff] = Screen('Flip', Display.win, When,[],Display.DontSync,Display.MultiFlip);
            if Debug.On == 1
                fprintf('Photodiode target was on for %.3f ms\n', (FixOff-FixOn)*1000);
            end
        end

        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win, [],[],Display.DontSync,Display.MultiFlip); 
        end

        %=========================== DRAW TEXTURE =========================
        function R_picOn(varargin)
            params = varargin{1,1};
            Eye = params(1);
            Pic.Rotation = params(2);
            if Eye < 2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);                     % Select blank eye
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
               	Eye = abs(floor(Eye-0.5));                                                              % For PTB, 0 = Left, 1 = Right
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);                     % Select which eye to draw to
                Screen('DrawTexture', Display.win, ImageTexture(Pic.Number==Pic.ImageIndx), [], Pic.Rect, Pic.Rotation, [], Pic.Contrast);     
                Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Pic.Rotation); 	% Apply Gaussian aperture mask
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye+1});
            elseif Eye == 2
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                    Screen('DrawTexture', Display.win, ImageTexture(Pic.Number==Pic.ImageIndx), [], Pic.Rect, Pic.Rotation, [], Pic.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Pic.Rotation); 	% Apply Gaussian aperture mask
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                 	Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
                end
            end
            [VBL ImageOn] = Screen('Flip', Display.win, [], 0, Display.DontSync,Display.MultiFlip);%, ImageOff+(ISI/1000)-FlipInt);         % Flip image frame
            Event(end+1,:) = [1, ImageOn, Event(end,2)-ImageOn];
        end

        %=========================== CLEAR TEXTURE ========================
        function R_picOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Mask.DestRect);
                Screen('DrawTexture', Display.win, Fix.Texture);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);%, [], 0, Display.DontSync,Display.MultiFlip);                                       % Clear stimulus  
            Event(end+1,:) = [0, ImageOff, Event(end,2)-ImageOff];
        end

        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
            disp(Event);
        end

        %========================= LOAD IMAGE VARIABLES ===================
        function R_loadPic(varargin)
            params=varargin{1,1};
            Pic.Number = params(1);
            Pic.ScaleX = params(2);
            Pic.ScaleY = params(3);
            Pic.PosX = params(4);
            Pic.PosY = params(5);
            Pic.Contrast = params(6);
            Mask.DestRect = round([0, 0, Pic.ScaleX, Pic.ScaleY]*Display.Pixels_per_deg(1))+[0 0 2 2];
            Pic.Rect = round([0, 0, Pic.ScaleX, Pic.ScaleY]*Display.Pixels_per_deg(1));
            Pic.Rect = CenterRect(Pic.Rect,Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            Mask.DestRect = CenterRect(Mask.DestRect, Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];           
        end
        
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
    end


%% ========================== apm_StereoFaces =============================
% Presents static 3D stimuli from a specified directory containing stereo 
% images in side-by-side format
%
    function ImageTuningHandle = apm_StereoFaces()

        ImageTuningHandle = ...
            {...
                @R_setFixColor;...
                @R_nextTrial;...
                @R_getBlock;...
                @R_getNextStim;...
                @R_getFixPosX;...
                @R_getFixPosY;...
                @R_fixOn;...
                @R_loadPic;...
                @R_preloadPics;...
                @R_picOn;...
                @R_picOff;...
                @R_fixOff;...
                @R_punishStim;... 
            };
        persistent ImageTexture;
        persistent Event;
        persistent Trial_num;
        persistent Block;
        persistent Pic;
        persistent Fix;

        Screen('BlendFunction', Display.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        %====================== SET STIMULUS DEFAULT PARAMETERS ==============
        Fix.Type    = 1;                                                    % 0 = dot; 1 = cross; 2 = square; 3 = crosshairs
        Fix.Eye     = 2;                                                 	% Present fixation binocularly
     	Fix.Drawn   = 0;
        Event = [0 0 0];
        if isfield(Fix,'Texture')
           Fix = rmfield(Fix,'Texture');
        end
        
     	%====================== LOAD IMAGE FILES AS TEXTURES ==================  
        function R_preloadPics(varargin)
            ExpType = 4;    % <<<<<<<<<<<<  Days 1-2 = ExpType 1; Day 3 = ExpType 2; Days 4-5 = ExpType 3
            
            if ~isfield(Pic, 'ConditionMatrix')
                FaceImageDir    = 'P:\murphya\MacaqueFace3D\BlenderFiles\Renders\';
                if ExpType < 4
                    MonkeyIDs   = 1;
                elseif ExpType >= 4
                    MonkeyIDs 	= [1,2,3,4,5];
                end

             	Elevations  = [-30, 0, 30];
             	Distances   = [-20, 0, 20];
                if ExpType == 1
                    Expressions = {'neutral'};
                    Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
                    Scales      = [10, 12, 20];
                elseif ExpType == 2
                  	Expressions = {'fear'};
                    Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
                    Scales      = [10, 12, 20];
                elseif ExpType == 3
                    Expressions = {'neutral','fear','lipsmack','threat'};
                    Azimuths    = [-30, 0, 30];
                    Scales      = [12];
                elseif ExpType == 4
                    Expressions = {'neutral'};
                    Azimuths    = [-30, 0, 30];
                    Scales      = [12];
                end
                Indx        = 1;
                for exp = 1:numel(Expressions)
                    for az = 1:numel(Azimuths)
                        for el = 1:numel(Elevations)
                            for d = 1:numel(Distances)
                                for s = 1:numel(Scales)
                                    for m = 1:numel(MonkeyIDs)
                                        Pic.ImgFilenames{Indx} = fullfile(FaceImageDir, sprintf('Monkey_%d', MonkeyIDs(m)), sprintf('Macaque_Id%d_%s_az%d_el%d_dist%d_sc%d.png',  MonkeyIDs(m), Expressions{exp}, Azimuths(az), Elevations(el), Distances(d), Scales(s)));
                                        if numel(MonkeyIDs) == 1
                                            Pic.ConditionMatrix(Indx,:) = [exp, az, el, d, s];
                                        elseif numel(MonkeyIDs)>1
                                            Pic.ConditionMatrix(Indx,:) = [m, az, el, d, s];
                                        end
                                        Indx = Indx+1;
                                    end
                                end
                            end
                        end
                    end
                end
                %<<<<<<<<< To include stereoscopic depth profile as a
                %variable (1) stereo congruent; 0) mono; -1) inconruent:
                if ExpType == 3
                    Pic.ConditionMatrix = [Pic.ConditionMatrix, ones(size(Pic.ConditionMatrix,1),1); Pic.ConditionMatrix, zeros(size(Pic.ConditionMatrix,1),1); Pic.ConditionMatrix, -1*ones(size(Pic.ConditionMatrix,1),1)];
                end
            end
            Pic.Order = GenerateDesign(size(Pic.ConditionMatrix,1), 100);
           
            params = varargin{1,1};
            if ~isfield(Pic, 'ImageIndx') || (params(1) ~= Pic.ImageIndx(1) || params(2) ~= Pic.ImageIndx(end))
                if ExpType <= 2
                    Pic.ImageIndx   = params(1):params(2);
                elseif ExpType == 3
                    Pic.ImageIndx   = 1:numel(Pic.ImgFilenames);      % <<<< TEMPORARY HACK FOR PRESENTING ALL EXPRESSIONS AT ORIGINAL SCALE ONLY
                elseif ExpType == 4
                    Pic.ImageIndx   = 1:numel(Pic.ImgFilenames); 
                end
                NoImages        = numel(Pic.ImageIndx);
                if exist('ImageTexture','var')                                  
                    Screen('Close', ImageTexture);                           	% close any previously loaded image textures
                end
                reverseStr = '';
                for n = 1:NoImages
                        LoadingText = sprintf('Loading image %d of %d... (#%d - %d)',n,NoImages,Pic.ImageIndx(1),Pic.ImageIndx(end));
    %                     fprintf([reverseStr, LoadingText]);
                        reverseStr = repmat(sprintf('\b'), 1, length(LoadingText));
                        DrawFormattedText(Display.win, LoadingText, 20, Display.Rect(4)-40, [0 0 0], []);% Draw text
                        Screen('Flip',Display.win);

                    [im, map, alpha] = imread(Pic.ImgFilenames{Pic.ImageIndx(n)});
                    im(:,:,4) = alpha;
                    ImageTexture(n) = Screen('MakeTexture', Display.win, double(im)/max(double(im(:)))*255);
                end
              	if size(Pic.ConditionMatrix,2) == 6                     % If stereo depth profile is a variable...
                    Pic.ImageIndx   = 1:(numel(Pic.ImgFilenames)*3);
                end
                Pic.SourceRect{1} = [0, 0, size(im,2)/2, size(im,1)];
                Pic.SourceRect{2} = Pic.SourceRect{1}+[size(im,2)/2, 0, size(im,2)/2, 0];
                fprintf('\nImages loaded sucessfully!\n\n');
                Screen('Fillrect', Display.win, Display.Background, Display.Rect);
                Screen('Flip', Display.win, [], 1);
            end

        end

        %======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset  = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size    = params(3)*Display.Pixels_per_deg(1);
            if Fix.Drawn == 0
                Fix = fixOn(Fix);
                Fix.Drawn = 1;
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
              	Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
            end
        	[VBL FixOn] = Screen('Flip', Display.win, [], 1);
        	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            When = FixOn +(1/Display.RefreshRate)+0.001;
         	[VBL FixOff] = Screen('Flip', Display.win, When,[],Display.DontSync,Display.MultiFlip);
            if Debug.On == 1
                fprintf('Photodiode target was on for %.3f ms\n', (FixOff-FixOn)*1000);
            end
        end

        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win, [],[],Display.DontSync,Display.MultiFlip); 
        end

        %========================== CLEAR FIXATION ========================
        function R_getNextStim(varargin)
            params      = varargin{1,1};
            TrialNumber = params(1)+1;
            if TrialNumber > numel(Pic.Order)
                TrialNumber = TrialNumber-numel(Pic.Order);
            end
            reply = sprintf('%d\n',Pic.Order(TrialNumber));
            
        end
            
        %=========================== DRAW TEXTURE =========================
        function R_picOn(varargin)
            params = varargin{1,1};
            Pic.Rotation = params(2);
            
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                if Pic.Stereo == 1
                    Screen('DrawTexture', Display.win, ImageTexture(Pic.Number), Pic.SourceRect{Eye}, [], Pic.Rotation, [], Pic.Contrast); 
                elseif Pic.Stereo == 0
                    Screen('DrawTexture', Display.win, ImageTexture(Pic.Number), Pic.SourceRect{1}, [], Pic.Rotation, [], Pic.Contrast); 
%                     Screen('FillOval', Display.win, [255 0 0], Photodiode.Rect{ceil(abs(Eye-2.5))});    % <<<<<< TEMP DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                elseif Pic.Stereo == -1
                    Screen('DrawTexture', Display.win, ImageTexture(Pic.Number), Pic.SourceRect{ceil(abs(Eye-2.5))}, [], Pic.Rotation, [], Pic.Contrast); 
                end
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOn] = Screen('Flip', Display.win, [], 0, Display.DontSync,Display.MultiFlip);%, ImageOff+(ISI/1000)-FlipInt);         % Flip image frame
            Event(end+1,:) = [1, ImageOn, Event(end,2)-ImageOn];
        end

        %=========================== CLEAR TEXTURE ========================
        function R_picOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Display.Rect);
                Screen('DrawTexture', Display.win, Fix.Texture);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);%, [], 0, Display.DontSync,Display.MultiFlip);                                       % Clear stimulus  
            Event(end+1,:) = [0, ImageOff, Event(end,2)-ImageOff];
        end

        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
        end

        %========================= LOAD IMAGE VARIABLES ===================
        function R_loadPic(varargin)
            params=varargin{1,1};
            Pic.Number  = params(1);
            Pic.ScaleX  = params(2);
            Pic.ScaleY  = params(3);
            Pic.PosX    = params(4);
            Pic.PosY    = params(5);
            Pic.Contrast = params(6);         

            if size(Pic.ConditionMatrix,2) == 6
                Pic.Stereo  = Pic.ConditionMatrix(Pic.Number,6);
                if Pic.Number > numel(Pic.ImgFilenames)
                    if rem(Pic.Number, numel(Pic.ImgFilenames)) ~= 0
                        Pic.Number = rem(Pic.Number, numel(Pic.ImgFilenames));
                    else
                        Pic.Number = numel(Pic.ImgFilenames);
                    end
                end
            else
             	Pic.Stereo  = params(7);
            end
        end
        
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
    end

%% ========================== apm_picTuningStereo =========================
% Presents static stimuli from a specified directory containing .png images
%
    function ImageTuningHandle = apm_picTuningStereo()

        ImageTuningHandle = ...
            {...
                @R_setFixColor;...
                @R_nextTrial;...
                @R_getBlock;...
                @R_getFixPosX;...
                @R_getFixPosY;...
                @R_fixOn;...
                @R_loadPic;...
                @R_preloadPics;...
                @R_picOn;...
                @R_picOff;...
                @R_fixOff;...
                @R_punishStim;... 
            };
        persistent ImageTexture;
        persistent MaskTex;
        persistent Event;
        persistent Trial_num;
        persistent Block;
        persistent Pic;
        persistent Mask;
        persistent Fix;
        persistent Params;
        persistent Cube

        Screen('BlendFunction', Display.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        
        %====================== SET STIMULUS DEFAULT PARAMETERS ==============
        Fix.Type    = 1;                                                   % 0 = dot; 1 = cross; 2 = square; 3 = crosshairs
        Fix.Eye     = 2;                                                    % Present fixation binocularly
     	Fix.Drawn   = 0;
        Event       = [0 0 0];
        if isfield(Fix,'Texture')
           Fix = rmfield(Fix,'Texture');
        end
        
        %==================== SET STEREO PARAMETERS
        Params.NoTrials             = 9;
        Params.IPD                  = 0.035;
        Params.Background           = [127 127 127];
        Params.PhysicalSizes       	= [0.06, 0.1, 0.16];                                        % Virtual physical sizes of images (m)
        Params.PositionsInDepth     = [0.10, 0, -0.10];                                         % Positions in depth relative to screen (m)
        Params.Distances            = repmat(Display.D, [1, numel(Params.PositionsInDepth)]) - Params.PositionsInDepth;
        Params.VergenceAngles       = atand(repmat(Params.IPD/2, [1, numel(Params.PositionsInDepth)])./Params.Distances);
        Params.DisparityOffsetM     = tand(Params.VergenceAngles).*Params.PositionsInDepth;     % Disparity offsets (m)
        Params.DisparityOffsetP     = Params.DisparityOffsetM*Display.Pixels_per_m(1);          % Disparity offsets (pixels)
        Params.DisparityOffsetD     = Params.DisparityOffsetP/Display.Pixels_per_deg(1);        % Disparity offsets (degrees)
        Params.PerspectiveScaling   = atand(repmat(0.1/2, [1, numel(Params.Distances)])./Params.Distances)*2;
        Params.Design               = nan(Params.NoTrials, 3);
        Params.Design(:,1)          = [1 1 1 2 2 2 3 3 3];
        Params.Design(:,2)          = repmat(1:numel(Params.DisparityOffsetP), [1, 3])';
        for p = 1:numel(Params.PhysicalSizes)
            Params.RetinalSizesD(p, :) = atand(repmat(Params.PhysicalSizes(p)/2, [1,numel(Params.PositionsInDepth)])./(Display.D-Params.PositionsInDepth));
            Params.RetinalSizesP(p, :) = Params.RetinalSizesD(p, :)*Display.Pixels_per_deg(1);
            for d = 1:numel(Params.DisparityOffsetP)
                Params.SourceRect{p, d}  = [0, 0, Params.RetinalSizesP(p, d), Params.RetinalSizesP(p, d)];
                Params.DestRect{1, p, d} = CenterRect(Params.SourceRect{p}, Display.Rect)+[Params.DisparityOffsetP(d), 0, Params.DisparityOffsetP(d), 0];
                Params.DestRect{2, p, d} = CenterRect(Params.SourceRect{p}, Display.Rect)-[Params.DisparityOffsetP(d), 0, Params.DisparityOffsetP(d), 0];
            end
        end
        
        
     	%====================== LOAD IMAGE FILES AS TEXTURES ==================  
        function R_preloadPics(varargin)
%             params          = varargin{1,1};
%             Pic.ImageIndx   = params(1):params(2);
%             ImageDir            = 'P:\murphya\PositionInDepth\2D_images';
            ImageDir            = 'P:\murphya\Stimuli\Processed';
         	ImageFormat         = '.png';
            ImageCategoryDirs   = dir(ImageDir);
            ImageCategoryDirs   = {ImageCategoryDirs([ImageCategoryDirs.isdir]).name};
            AllImageFiles       = wildcardsearch(ImageDir, ['*',ImageFormat]);
            NoImages            = numel(AllImageFiles);
            
            ImageRadius     = 5*Display.Pixels_per_deg(1);                      % Default radius of image size BEFORE drawing based on size input                 
            TextureWindow   = [2*ImageRadius, 2*ImageRadius];               	% Set size of the stimulus before mask is applied (degrees)
        
            if exist('ImageTexture','var')                                  
                Screen('Close', ImageTexture);                                  % close any previously loaded image textures
            end

            reverseStr = '';
            for n = 1:NoImages
%                 if Debug.On == 1
                    LoadingText = sprintf('Loading image %d of %d...',n,NoImages);
%                     fprintf([reverseStr, LoadingText]);
%                     reverseStr = repmat(sprintf('\b'), 1, length(LoadingText));
                   	DrawFormattedText(Display.win, LoadingText, 20, Display.Rect(4)-40, [0 0 0], []);% Draw text
                    Screen('Flip',Display.win);
%                 end
                [Img ImgCmap] = imread(AllImageFiles{n});                                   % Load image file
                ImgDim = size(Img);                                                         % Get image dimensions
                if min(ImgDim([1,2])) ~= 2*ImageRadius                                    	% if smallest image dimesnion is not requested size...
                    scale = 2*ImageRadius/min(ImgDim([1 2]));                               % Resize image so smallest dimension fits
                    Img = imresize(Img, scale);
                    ImgDim = size(Img);                                                     % Get new Img image dimensions
                end         
                if max(ImgDim([1,2])) > 2*ImageRadius                                   	% If largest image dimension is too large...
                    Crop = (max(ImgDim([1,2]))-(2*ImageRadius))/2;
                    if find(ImgDim==max(ImgDim([1,2])))==1
                        Img(end-Crop:end,:,:) = [];
                        Img(1:Crop,:,:) = [];
                    elseif find(ImgDim==max(ImgDim([1,2])))==2
                        Img(:,end-Crop:end,:) = [];
                        Img(:,1:Crop,:) = [];
                    end
                end
                Pic.ImageRect = [0,0,size(Img,2),size(Img,1)];
                ImageTexture(n) = Screen('MakeTexture', Display.win, double(Img));          % Convert image to PTB texture
            end
            fprintf('\nImages loaded sucessfully!\n\n');
            Screen('Fillrect', Display.win, Display.Background, Display.Rect);
            Screen('Flip', Display.win, [], 1);
            
            %========================= CREATE ALPHA CHANNEL MASK
            Mask.Dim            = TextureWindow+4;                                	% Mask will be 2 pixels larger than the stimulus window
            Mask.ApRadius       = (min(Mask.Dim)/2)-8;
            Mask.Colour         = Display.Background(1);
            Mask.Edge           = 2;                                              	% Cosine mask edge
            Mask.Taper          = 0.2;                                            	% Cosine edge tapers off over 20% of aperture radius
            MaskTex             = GenerateAlphaMask(Mask, Display);
            Mask.Rect           = [0 0 Mask.Dim];                                 	% Mask size
            Mask.DestRect       = CenterRect(Mask.Rect, Display.Rect);            	% Destination for mask is centred in screen window
            
            %========================= CREATE BACKGROUND BORDER
%             Square.Density      = 0.7;                                            	% proportion of the square area to fill with squares
%             Square.InnerBorder  = 4*Display.Pixels_per_deg;                         % size of border between surrounding squares and stimulus (pixels)
%             Square.Size         = 0.5*Display.Pixels_per_deg(1);                    % dimensions of each square (pixels)
%             Square.Filled       = 0;                                                % 0 = outline; 1 = solid
%             Stim.Window         = Mask.DestRect+Params.DisparityOffsetP(d);       
%             Stim.Background     = Display.Background;
%             Pic.BackgroundTexture   = BackgroundSquares(Display, Stim, Square);
            
            
            %==================== GENERATE SURROUND TEXTURE
            Cube.IPD            = Params.IPD;
            Cube.Density        = 1;                                            % proportion of the texture area to fill with squares
            Cube.Number        = 300;
            Cube.BlankRect      = Mask.DestRect;                                % central rectangle in which no cubes will be drawn
            Cube.InnerBorder    = 2*Display.Pixels_per_deg;                     % size of border between surrounding cubes and stimulus (pixels)
            Cube.Size           = 0.015*Display.Pixels_per_m;                   % dimensions of each cube (pixels)
            Cube.Texture        = 2;                                            % 0 = wire frame; 1 = solid; 2 = window texture; 3 = Rubiks texture; 4 = custom texture
            Cube.DepthRange     = [-0.2, 0.2];                                  % [near, far] depth limits (metres)
            Cube.Background     = [Params.Background, 255];                   	% background color (RGBA)
            [Cube.BackgroundTextures] = BackgroundCubes(Display, Cube);      % Generate background texture
%             BackgroundCubes(Display, Cube);
        end

        %======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            if Fix.Drawn == 0
                Fix = fixOn(Fix);
                Fix.Drawn = 1;
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
%                 Screen('DrawTexture', Display.win, Pic.BackgroundTexture);
                Screen('DrawTexture', Display.win, Cube.BackgroundTextures(Eye));
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
              	Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
            end
        	[VBL FixOn] = Screen('Flip', Display.win, [], 1);
        	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            When = FixOn +(1/Display.RefreshRate)+0.001;
         	[VBL FixOff] = Screen('Flip', Display.win, When,[],Display.DontSync,Display.MultiFlip);
            if Debug.On == 1
                fprintf('Photodiode target was on for %.3f ms\n', (FixOff-FixOn)*1000);
            end
        end

        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
%                 Screen('DrawTexture', Display.win, Pic.BackgroundTexture);
                Screen('DrawTexture', Display.win, Cube.BackgroundTextures(Eye));
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
            end 
            [VBL FixOff] = Screen('Flip', Display.win, [],[],Display.DontSync,Display.MultiFlip); 
        end

        %=========================== DRAW TEXTURE =========================
        function R_picOn(varargin)
            params          = varargin{1,1};
            Eye             = params(1);
            Pic.Rotation    = params(2);
            trial           = params(3)+1;
            if Eye < 2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);                     % Select blank eye
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
               	Eye = abs(floor(Eye-0.5));                                                              % For PTB, 0 = Left, 1 = Right
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);                     % Select which eye to draw to
                
                Screen('DrawTexture', Display.win, ImageTexture(Pic.Number), [], Pic.Rect, Pic.Rotation, [], Pic.Contrast);     
                Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Pic.Rotation); 	% Apply Gaussian aperture mask
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye+1});
            elseif Eye == 2
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
%                     Screen('DrawTexture', Display.win, Pic.BackgroundTexture);
                 	Screen('DrawTexture', Display.win, Cube.BackgroundTextures(Eye));
                    Screen('DrawTexture', Display.win, ImageTexture(Pic.Number), [], Params.DestRect{Eye, Params.Design(trial,1), Params.Design(trial,2)}, Pic.Rotation, [], Pic.Contrast);
%                     Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Params.DestRect{Eye, Params.Design(trial,1), Params.Design(trial,2)}, Pic.Rotation); 	% Apply Gaussian aperture mask
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                 	Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
                end
            end
            [VBL ImageOn] = Screen('Flip', Display.win, [], 0, Display.DontSync,Display.MultiFlip);%, ImageOff+(ISI/1000)-FlipInt);         % Flip image frame
            Event(end+1,:) = [1, ImageOn, Event(end,2)-ImageOn];
        end

        %=========================== CLEAR TEXTURE ========================
        function R_picOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, Cube.BackgroundTextures(Eye));
%                 Screen('DrawTexture', Display.win, Pic.BackgroundTexture);
%                 Screen('FillRect', Display.win, Display.Background, Mask.DestRect);
                Screen('DrawTexture', Display.win, Fix.Texture);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);%, [], 0, Display.DontSync,Display.MultiFlip);                                       % Clear stimulus  
            Event(end+1,:) = [0, ImageOff, Event(end,2)-ImageOff];
        end

        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
            disp(Event);
        end

        %========================= LOAD IMAGE VARIABLES ===================
        function R_loadPic(varargin)
            params=varargin{1,1};
            Pic.Number      = params(1);
            Pic.ScaleX      = params(2);
            Pic.ScaleY      = params(3);
            Pic.PosX        = params(4);
            Pic.PosY        = params(5);
            Pic.Contrast    = params(6);
            Mask.DestRect   = round([0, 0, Pic.ScaleX, Pic.ScaleY]*Display.Pixels_per_deg(1))+[0 0 2 2];
            Pic.Rect        = round([0, 0, Pic.ScaleX, Pic.ScaleY]*Display.Pixels_per_deg(1));
            Pic.Rect        = CenterRect(Pic.Rect,Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            Mask.DestRect   = CenterRect(Mask.DestRect, Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];           
        end
        
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
    end

%% =========================== apm_motion =================================
% Presents motion stimuli from a specified directory containing .avi movies
%==========================================================================
    function ImageTuningHandle = apm_motion()

        ImageTuningHandle = ...
            {...
                @R_setFixColor;...
                @R_nextTrial;...
                @R_getBlock;...
                @R_getFixPosX;...
                @R_getFixPosY;...
                @R_fixOn;...
                @R_loadMov;...
                @R_preloadMovs;...
                @R_movOn;...
                @R_movOff;...
                @R_fixOff;...
                @R_punishStim;... 
            };
        persistent MovieTextures;
%         persistent MaskTex;
        persistent Trial_num;
        persistent Block;
        persistent Mov;
        persistent Mask;
        persistent Fix;
        
        
        %====================== SET STIMULUS DEFAULT PARAMETERS ==============
        if isfield(Fix,'Texture')
           Fix = rmfield(Fix,'Texture');
        end
        Fix.Drawn = 0;
        Fix.Type = 1;                                                   % 0 = dot; 1 = cross; 2 = square; 3 = crosshairs
        Fix.Eye = 2;                                                    % Present fixation binocularly
        Fix.FixAlwaysOn = 1;                                            % Keep fixation marker on during stimulus poresentation?
        Mov.Loop = 1;                                                   
      	Mov.FullScreen = 0;                                             % use full screen?
        Mov.Rate = 1;                                                   % Playback speed (proportion: -1 = reverse, <1 = slow, 1 = normal, >1 = fast)
        Mov.Volume = 0;                                                 % Audio track volume (1 = full, 0 = mute)
        Mov.MaintainAR = 0;                                             % Adjust aspect ratio of movie to fit Mov.Dims?
        Mov.Background = Display.Background;
        Mov.StartTime = 0;                                              % Always begin playback from frame 1
        Mov.TimeInFrames = 1;                                           % Time is specified in frames
        Mov.PreloadSecs = 2;                                            % Number of seconds to preload
        Mov.Blocking = 1;                                               % Wait for next frame to become available
        Mov.Async = 4;                                                  % 4 = Allow asynchronous (i.e. background) loading of movies 
        Mov.PreloadToTextures = 1;                                      % Preload movie frames into texture?
        
     	%====================== LOAD IMAGE FILES AS TEXTURES ==================  
        function R_preloadMovs(varargin)
            params = varargin{1,1};             
            Mov.Indx = params(1):params(2);                                 % Specify which movie files to load
            Mov.NoClips = numel(Mov.Indx);                                           
            Mov.PlaybackTime = 2;                                           % Set playback duration (seconds)
            Mov.WindowSize = [10 10]*Display.Pixels_per_deg(1);          	% Default movie screen size BEFORE drawing based on size input                 
            Mov.Dims = Mov.WindowSize;                                       % Movie dimensions (if Mov.Fullscreen = 0)

            if exist('MovieTextures','var')    
                for n = 1:numel(MovieTextures)
                    Screen('Close', MovieTextures{n});                           	% close any previously loaded image textures
                    if isfield(Mov,'Handle')
                        Screen('CloseMovie',Mov.Handle{n});
                    end
                end
            end
            Mov.Format = '.avi';
            reverseStr = '';
            for n = 1:Mov.NoClips
%                 if Debug.On == 1
                    LoadingText = sprintf('Loading movie clip %d of %d... (Movie #%d)',n,Mov.NoClips,Mov.Indx(n));
                    fprintf([reverseStr, LoadingText]);
                    reverseStr = repmat(sprintf('\b'), 1, length(LoadingText));
                    DrawFormattedText(Display.win, LoadingText, 20, Display.Rect(4)-40, [0 0 0], []);           % Draw text
                    Screen('Flip',Display.win);
%                 end
                Mov.FileName{n} = fullfile(MotionDir,[num2str(Mov.Indx(n)),Mov.Format]);                        % Construct filename of next image
                
                [Mov.Handle{n}, Mov.duration{n}, Mov.fps{n}, Mov.width{n}, Mov.height{n}, Mov.Count{n}, Mov.AR]= Screen('OpenMovie', Display.win, Mov.FileName{n},Mov.Async,Mov.PreloadSecs); 
                
                %============== PRELOAD FRAMES TO TEXTURES
                if Mov.PreloadToTextures == 1
                    MovieTextures{n} = zeros(1, Mov.Count{n});                                              	% preallocate cells for texture handles
                    Screen('SetMovieTimeIndex', Mov.Handle{n}, Mov.StartTime, Mov.TimeInFrames);
                    Screen('PlayMovie',Mov.Handle{n}, 100, 0, 0);                                               % For async decoding, play movie at x100 speed with no looping or sound
                    for f = 1:Mov.Count{n}
                        MovieTextures{n}(f) = Screen('GetMovieImage', Display.win, Mov.Handle{n},Mov.Blocking);
                    end
                end
                
                %============== CALCULATE MOVIE AND SCREEN RECTANGLES
                Mov.SourceRect{2} = [0 0 Mov.width{n}, Mov.height{n}];
                if Mov.MaintainAR == 0
                    if Mov.FullScreen == 1
                        Mov.DestRect{n} = Display.Rect;
                    elseif Mov.FullScreen == 0
                        Mov.DestRect{n} = [0 0 Mov.Dims]*Display.Pixels_per_deg(1);
                    end
                elseif Mov.MaintainAR == 1
                    if Mov.FullScreen == 1
                        Mov.WidthDeg = Display.Rect(3);
                    else
                        Mov.WidthDeg = Mov.Dims(1)*Display.Pixels_per_deg(1);
                    end
                    Mov.DestRect{n} = (Mov.SourceRect{2}/Mov.width{n})*Mov.WidthDeg;
                end
                if ~isempty(find(Mov.DestRect{n} > Display.Rect))
                    Mov.DestRect{n} = Mov.DestRect{n}*min(Display.Rect([3, 4])./Mov.DestRect{n}([3, 4]));
                    fprintf('Requested movie size > screen size! Defaulting to maximum size.\n');
                end
                Mov.DestRect{n} = CenterRect(Mov.DestRect{n}, Display.Rect);
            end
            fprintf('\nMotion stimuli loaded sucessfully!\n\n');
            Screen('Fillrect', Display.win, Display.Background, Display.Rect);
            Screen('Flip', Display.win, [], 1);
            
%             %========================= CREATE ALPHA CHANNEL MASK
%             Mask.Dim = TextureWindow+4;                                                 % Mask will be 2 pixels larger than the stimulus window
%             Mask.ApRadius = (min(Mask.Dim)/2)-8;
%             Mask.Colour = Display.Background(1);
%             Mask.Edge = 2;                                                              % Cosine mask edge
%             Mask.Taper = 0.2;                                                           % Cosine edge tapers off over 20% of aperture radius
%             MaskTex = GenerateAlphaMask(Mask, Display);
%             Mask.Rect = [0 0 Mask.Dim];                                                 % Mask size
%             Mask.DestRect = CenterRect(Mask.Rect, Display.Rect);                     	% Destination for mask is centred in screen window
        end

        
        %======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            if Fix.Drawn == 0
                Fix = fixOn(Fix);
                Fix.Drawn = 1;
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
            end
        	[VBL FixOn] = Screen('Flip', Display.win, [], 1);
        	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            When = FixOn +(1/Display.RefreshRate)+0.001;
         	[VBL FixOff] = Screen('Flip', Display.win, When,[],Display.DontSync,Display.MultiFlip);
            if Debug.On == 1
                fprintf('Photodiode target was on for ~%.3f ms\n', (FixOff-FixOn)*1000);
            end
        end

        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win, [],[],Display.DontSync,Display.MultiFlip); 
        end

        %========================== PRESENT STIMULUS ======================
        function R_movOn(varargin)
            params = varargin{1,1};
            Eye = params(1);
            n = find(Mov.Indx==Mov.Number);
            Screen('SetmovieTimeIndex',Mov.Handle{n},Mov.StartTime,Mov.TimeInFrames);
            Screen('PlayMovie',Mov.Handle{n}, Mov.Rate, Mov.Loop, Mov.Volume);
            
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            EndMovie = 0;
            f = 1;
            while EndMovie == 0
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    parseCommand(commandIn);                % parse command
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn,'R_movOff')         	% Exits loop if movie is being turned off.     
                    EndMovie =1;
                    break
                end
                if Mov.PreloadToTextures == 1               % If movies were preloaded to PTB textures...
                    MovieTex = MovieTextures{n}(floor(f));
                    f = f + (1/(Display.RefreshRate/Mov.fps{n}));
                    if floor(f) > numel(MovieTextures{n})
                        f = 1;
                    end
               	elseif Mov.PreloadToTextures == 0           
                    MovieTex = Screen('GetMovieImage', Display.win, Mov.Handle{n},Mov.Blocking);
                end
                if MovieTex < 0
                    break;
                end
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                    Screen('FillRect',Display.win, Mov.Background);
                    if MovieTex > 0
                        Screen('DrawTexture', Display.win, MovieTex);%, Mov.SourceRect{Eye});%, Mov.DestRect); % <<<<<<<<<<<<<< TEMP FUDGE
                        Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
                    end
                    if Fix.FixAlwaysOn == 1
                        Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect, Fix.Rotation(Eye)); 
                    end
                end
                [VBL FrameOn] = Screen('Flip', Display.win);
%             	Screen('Close', MovieTex);
                [keyIsDown,secs,keyCode] = KbCheck;                     % Check keyboard for 'escape' press        
                if keyIsDown && keyCode(Key.Exit) == 1                   % Press Esc for abort
                    EndMovie = 1;
                end
            end
            R_movOff;
        end
        
        %=========================== CLEAR TEXTURE ========================
        function R_movOff
            if nargin == 1
                Screen('PlayMovie', Mov.Handle{Mov.Number}, 0, Mov.Loop, Mov.Volume);       % Stop playback
%                 Screen('CloseMovie', Mov.Handle{Mov.Number});
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Mov.Background); 
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect, Fix.Rotation(Eye)); 
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus  
        end

        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
        end

        %========================= LOAD MOVIE VARIABLES ===================
        function R_loadMov(varargin)
            params=varargin{1,1};
            Mov.Number = params(1);
            Mov.ScaleX = params(2);
            Mov.ScaleY = params(3);
            Mov.PosX = params(4);
            Mov.PosY = params(5);
            Mov.Contrast = params(6);
            Mask.DestRect = round([0, 0, Mov.ScaleX, Mov.ScaleY]*Display.Pixels_per_deg(1))+[0 0 2 2];
            Mov.Rect = round([0, 0, Mov.ScaleX, Mov.ScaleY]*Display.Pixels_per_deg(1));
            Mov.Rect = CenterRect(Mov.Rect,Display.Rect) + [Mov.PosX,Mov.PosY,Mov.PosX,Mov.PosY];
            Mask.DestRect = CenterRect(Mask.DestRect, Display.Rect) + [Mov.PosX,Mov.PosY,Mov.PosX,Mov.PosY];           
        end
        
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
    end



%% =========================== apm_binorivOKN =============================
% Presents a static image to one eye and a moving image to the other eye,
% allowing perception to be deduced from the presence (or absence) of
% optokinetic nystagmus eve movements.
%
% OPTIONAL SETTINGS:
%       Intermittent:   Presents stimuli intermittently for perceptual stabilization
%       Flicker:        Flickers stimuli at different rates in each eye for frequency tagging
%       OccularSwitch: 	Switches stimuli rapidly between the two eyes for stimulus rivalry
%
%==========================================================================

    function BinorivHandle = apm_binoriv()

        BinorivHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_dotsAndPicOn;...
                @R_setPic;...
                @R_setDots;...
                @R_setDotsMove;...
                @R_setDotsDir;...
                @R_dotsAndPicOff;...
                @R_getBlock;...
            };
        
        persistent Fix;
        persistent Dots;
        persistent Mask;
        persistent Frame;
        persistent Loop;
        persistent DotsTexture;
        persistent ImageTexture;
        persistent BackgroundTexture;
        persistent MaskTex;
        persistent Flicker;
        persistent Intermittent;
        persistent BackgroundOn;

        %======================= SET STIMULUS PARAMETERS ======================
        BackgroundOn = 0;
        Fix.On = 0;                                     % Keep fixation on during stimulus presentation?
        Fix.Type = 3;                                   % square fixation
        Fix.Eye = 2;                                  	% Present fixation binocularly
        Fix.LineWidth = 3;

        Pic.Rotation = 0;
        Frame = 1;
        FrameTotal = 1;
        Loop = 0;
        
        %============= FLICKER FREQUENCY TAGGING?
        Flicker.On = 0;                                 % Add frequency tagging to stimuli?
        Flicker.RateHzL = 6;                            % Set flicker rate (Hz)
        Flicker.RateHzR = 15;     
        Flicker.DurL = (1/Flicker.RateHzL)*Display.RefreshRate; % Convert to frames
        Flicker.DurR = (1/Flicker.RateHzR)*Display.RefreshRate; % Convert to frames
        Flicker.Eye1 = 1;
        Flicker.Eye2 = 1;
        
        %============= INTERMITTENT PRESENTATION?
        Intermittent.On = 0;                
        Intermittent.OnDur = 1;             % Set duration of stimulus on period (s)
        Intermittent.OffDur = 1;            % Set duration of blank period (s)
        Intermittent.StimOn = 1;            % Start trial with stimulus 'on'
         
     	Stim.Window = [0 0 500 500];
        Stim.Background = Display.Background;

        
     	%======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL FixOn] = Screen('Flip', Display.win, [], 1);
            Fix = fixOn(Fix);
         	[VBL FixOn] = Screen('Flip', Display.win);
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
         %========================== CLEAR FIXATION ========================
        function R_fixOff
            if isfield(Fix,'DestRect')
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                    Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
                end 
                [VBL FixOff] = Screen('Flip', Display.win, [], 1);                                         % Clear stimulus 
            end
        end
        
        %======================== PRESENT STIMULI =========================
        function R_dotsAndPicOn(varargin)
            params=varargin{1,1};
            Eye = params(1);
            Intermittent.On = params(2);
            Intermittent.OnDur = params(3)/1000;
            Intermittent.OffDur = params(4)/1000;
            
            Frame = 1;
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            
            while 1
             	commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    parseCommand(commandIn);                % parse command
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn,'R_dotsandPicOff')     % Exits loop if trial is quit or ended
                    break
                end
                 
                %===================== DRAW STIMULI =======================
                E = Eye;
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                if Flicker.Eye1 == 1 && Intermittent.StimOn == 1
                    Screen('DrawTexture', Display.win, DotsTexture(Frame), [], Pic.Rect, Dots.DirMean, [], Dots.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Dots.DirMean);                  % Apply aperture mask
                    Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{E+1});
               	else
                    Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{E+1});
                end
                E = abs(floor(Eye-0.5));
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E);                                       % Select which eye to draw to
                Screen('DrawTexture', Display.win, BackgroundTexture);
                if Flicker.Eye2 == 1 && Intermittent.StimOn == 1
                    Screen('DrawTexture', Display.win, ImageTexture, [], Pic.Rect, Pic.Rotation, [], Pic.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Pic.Rotation);                  % Apply aperture mask
                    Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{E+1});
                else
                    Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{E+1});
                end
                if Fix.On == 1 || Intermittent.StimOn == 0
                    for E = 1:2
                        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E-1);
                        Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect, Fix.Rotation(E)); 
                    end
                end
                [VBL FrameOn] = Screen('Flip', Display.win);
                if Frame == 1 && Loop == 0
                    Intermittent.SwitchOn = FrameOn;
                end

                %=================== FREQUENCY TAGGING
                if Flicker.On == 1
                    if mod(ceil(FrameTotal/Flicker.DurL),2)
                        Flicker.Eye1 = 1;                   % Flicker on
                    else
                        Flicker.Eye1 = 0;                   % Flicker off
                    end
                    if mod(ceil(FrameTotal/Flicker.DurR),2)
                        Flicker.Eye2 = 1;                   % Flicker on
                    else
                        Flicker.Eye2 = 0;                   % Flicker off
                    end
                end
                
                %================== INTERMITTENT PRESENTATION
                if Intermittent.On == 1
                    if Intermittent.StimOn == 1 && GetSecs >= Intermittent.SwitchOn+Intermittent.OnDur
                        Intermittent.StimOn = 0;
                        Intermittent.SwitchOff = GetSecs;
                    elseif Intermittent.StimOn == 0 && GetSecs >= Intermittent.SwitchOff+Intermittent.OffDur
                        Intermittent.StimOn = 1;
                        Intermittent.SwitchOn = GetSecs;
                    end
                end
                
                %================== ANIMATION LOOP
                Frame = Frame+1;
                FrameTotal = (Loop*Dots.NrFrames)+Frame;
                if Frame == Dots.NrFrames
                    Frame = 1;
                    Loop = Loop+1;
                end
                [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                if keyIsDown && keyCode(Key.Exit)
                    break;
                end
            end
            R_dotsAndPicOff;
            if exist('DotsTexture','var')
                Screen('Close', DotsTexture);
            end
            Screen('Close', ImageTexture);
        end
        
        %====================== SET IMAGE PARAMETERS ======================
        function R_setPic(varargin)
            params=varargin{1,1};
            if isfield(Pic, 'ScaleX')                                   % If Scale fields already exist
                if params(2) ~= Pic.ScaleX || params(3) ~= Pic.ScaleY   % If requested image size has changed...
                    NewBackground = 1;
                else
                    if BackgroundOn == params(7)
                        NewBackground = 0;
                    else
                        NewBackground = 1;
                    end
                end
            else
                NewBackground = 1;
            end
            Pic.Number = params(1);
            Pic.ScaleX = params(2);
            Pic.ScaleY = params(3);
            Pic.PosX = params(4);
            Pic.PosY = params(5);
            Pic.Contrast = params(6);
            BackgroundOn = params(7);
            Pic.Category = params(8);
            ImageRadius = round((Pic.ScaleX/2)*Display.Pixels_per_deg(1));
            
            AllImageTypes = regexp(genpath(ImageDir),['[^;]*'],'match');                % Get list of image subcategory folders
            CategoryDir = AllImageTypes{Pic.Category+1};                                % Get path of requested image category
            [x ImageCategory] = fileparts(CategoryDir);                                 % Get text string for category
            cd(CategoryDir);                                                            % Change to that folder
            Images = dir('*.jpg');                                                      % Get all jpg images
            [Img ImgCmap] = imread(fullfile(CategoryDir, Images(Pic.Number).name));
            ImgDim = size(Img);                                                         % Get image dimensions
            if min(ImgDim([1,2])) ~= 2*ImageRadius                                    	% if smallest image dimesnion is not requested size...
                scale = 2*ImageRadius/min(ImgDim([1 2]));                               % Resize image so smallest dimension fits
                Img = imresize(Img, scale);
                ImgDim = size(Img);                                                     % Get new Img image dimensions
            end         
            if max(ImgDim([1,2])) > 2*ImageRadius                                   	% If largest image dimension is too large...
                Crop = (max(ImgDim([1,2]))-(2*ImageRadius))/2;
                if find(ImgDim==max(ImgDim([1,2])))==1
                    Img(end-Crop:end,:,:) = [];
                    Img(1:Crop,:,:) = [];
                elseif find(ImgDim==max(ImgDim([1,2])))==2
                    Img(:,end-Crop:end,:) = [];
                    Img(:,1:Crop,:) = [];
                end
            end
            Pic.ImageRect = [0,0,size(Img,2),size(Img,1)];
            Pic.Rect = CenterRect(Pic.ImageRect,Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            ImageTexture = Screen('MakeTexture', Display.win, double(Img));             % Convert image to PTB texture
            
            %========================= CREATE ALPHA CHANNEL MASK ======================
            TextureWindow = [2*ImageRadius, 2*ImageRadius];                     % Set size of the stimulus before mask is applied (degrees)
            Dots.Window = TextureWindow;
            Stim.Background = Display.Background(1);
            Mask.Dim = TextureWindow+4;                                         % Mask will be 2 pixels larger than the stimulus window
            Mask.ApRadius = (min(Mask.Dim)/2)-8;
            Mask.Colour = Display.Background(1);
            Mask.Edge = 0;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
            MaskTex = GenerateAlphaMask(Mask, Display);
            Mask.Rect = [0 0 Mask.Dim];                                      	% Mask size
            Mask.DestRect = CenterRect(Mask.Rect, Display.Rect);            	% Destination for mask is centred in screen window
            Mask.DestRect = Mask.DestRect + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            Stim.Background = Display.Background;
          	Stim.Window = [0 0 Mask.Dim];
            if BackgroundOn == 1 
                if NewBackground == 1
                    BackgroundTexture = BackgroundSquares(Display, Stim);        	% Generate squares background
             	end
            elseif BackgroundOn == 0
                BackgroundTexture = Screen('MakeTexture', Display.win, ones(Display.Rect([4,3]))*Display.Background(1));
            end
        end
        
        %========================== SET DOT PARAMETERS ====================
        function R_setDots(varargin)
            params=varargin{1,1};
            Dots.Num = params(1);
            Dots.Colour = params([2,3,4])*255;
            Dots.Contrast = params(4);
            if Penalty.On == 1
                Dots.RandColour = 1;        % Set to 1 for random coloured dots for training
            else
                Dots.RandColour = 0;        % Set to 1 for random coloured dots for training
            end
        end
        
        function R_setDotsMove(varargin)
            params=varargin{1,1};
            Dots.Signal = params(1);
            Dots.SpeedMean = params(2);
            Dots.SpeedStd = params(3);
            Dots.LifeMean = params(4);
            Dots.LifeStd = params(5);
            Dots.Lifetime = Dots.LifeMean;
        end
        
        function R_setDotsDir(varargin)
            params=varargin{1,1};
            Dots.DirMean = params(1)-90;
            Dots.DirStd = params(2);
            Dots.SizeMean = params(3)*Display.Pixels_per_deg(1);
            Dots.SizeStd = params(4)*Display.Pixels_per_deg(1);
            
         	%================ GENERATE TRANSPARENT MOTION TEXTURES ============
            Dots.Velocity = Dots.SpeedMean*Display.Pixels_per_deg(1);
            Dots.Background = Display.Background;
            Dots.DrawAngle = Dots.DirMean;
            Dots.Size = ([Dots.SizeMean, Dots.SizeStd, 1]);
            DotsTexture = GenerateRDK(Dots, Display);
            Dots.NrFrames = numel(DotsTexture);
        end
        
        %=========================== CLEAR TEXTURE ========================
        function R_dotsAndPicOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win, [], 1);                                       % Clear stimulus 
            
            fprintf('There are %d dot textures open...\n', numel(DotsTexture));
            for i = 1:numel(DotsTexture)
                Screen('Close', DotsTexture(i));
            end
            Screen('Close', DotsTexture);
        end
        
        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
        end
    end


%% =========================== apm_binorivPA ==============================
% Physically alternates presentation of different images to the two eyes,
% simulating perception during binocular rivalry. Timing of alternations
% can either be randomized within a specified range or replicate the timing
% of alternations recorded during a previous binocular rivalry block.
%
%
%==========================================================================

    function BinorivHandle = apm_binorivPA()

        BinorivHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_stimOn;...
                @R_setPic;...
                @R_setDots;...
                @R_setDotsMove;...
                @R_setDotsDir;...
                @R_dotsAndPicOff;...
                @R_getBlock;...
                @R_Transition;...
            };
        
        persistent Fix;
        persistent Dots;
        persistent Pic;
        persistent Mask;
        persistent Frame;
        persistent Loop;
        persistent DotsTexture;
        persistent ImageTexture;
        persistent BackgroundTexture;
        persistent MaskTex;
        persistent Flicker;
        persistent Intermittent;
        persistent Transition;

        %======================= SET STIMULUS PARAMETERS ====================== 
        Frame = 1;
        BackgroundOn = 0;
        Fix.Type = 3;
        Fix.Eye = 2;
        Fix.On = 0;
        Pic.Rotation = 0;
        Pic.MinContrast = 0;
        Dots.MinContrast = 0;
        Flicker.On = 0;                 % Add frequency tagging to stimuli?
        Flicker.RateHzL = 6;            % Set flicker rate (Hz)
        Flicker.RateHzR = 15;     
        Flicker.DurL = (1/Flicker.RateHzL)*Display.RefreshRate; % Convert to frames
        Flicker.DurR = (1/Flicker.RateHzR)*Display.RefreshRate; % Convert to frames
        Flicker.Eye1 = 1;
        Flicker.Eye2 = 1;
        Frame = 1;
        FrameTotal = 1;
        Loop = 0;
        Intermittent.On = 0;                
        Intermittent.OnDur = 1;             % Set duration of stimulus on period (s)
        Intermittent.OffDur = 1;            % Set duration of blank period (s)
        Intermittent.StimOn = 1;            % Start trial with stimulus 'on'
        BackgroundTexture = Screen('MakeTexture', Display.win, ones(Display.Rect([4,3]))*Display.Background(1));

        %======================== PRESENT FIXATION ========================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL FixOn] = Screen('Flip', Display.win, [], 1);
            Fix = fixOn(Fix);
        	[VBL FixOn] = Screen('Flip', Display.win);
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
         %========================== CLEAR FIXATION ========================
        function R_fixOff
%             for Eye = 1:2
%                 currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
%                 Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
%             end 
%             [VBL FixOff] = Screen('Flip', Display.win, [], 1);                                         % Clear stimulus  
        end
        
        %======================== PREPARE TRANSITION ======================
        function R_Transition(TransitionDuration)
            disp(TransitionDuration)
            R_stimOn(TransitionDuration);
%             Transition.NoFrames = TransitionDuration/1000*Display.RefreshRate;
%             Transition.Frames = 1:Transition.NoFrames;
%             switch Transition.Type 
%                 case 1 %============== LINEAR TRANSITION
%                     Transition.Increments = Pic.Contrast/Transition.Frames;
%                     
%                 case 2 %============== CUMULATIVE GAUSSIAN
%                     Transition.Mu = Transition.NoFrames/2;
%                     Transition.Sigma = Transition.NoFrames/5;
%                     Transition.Increments = normcdf(Transition.Frames,Transition.Mu,Transition.Sigma);
%                     
%                 case 3 %============== SINUSOID
%                     
%             end
        end
        
        %======================== PRESENT STIMULI =========================
        function R_stimOn(varargin)
            params = varargin{1,1};
            Eye = params(1);
            Stim = params(2);
            Transition.Type = params(3);        % 
            Transition.Duration = params(4);    % Duration of transition period (ms)
            Transition.Frames = round(Transition.Duration/1000*Display.RefreshRate);
            switch Transition.Type 
                case 1 %============== IMMEDIATE STEP TRANSITION
                    Transition.Increments = repmat(1,[1,Transition.Frames]);
                    
                case 2 %============== LINEAR TRANSITION
                    Transition.Increments = repmat(1/Transition.Frames,[1,Transition.Frames]);
                    
                case 3 %============== GAUSSIAN
                    Transition.Mu = Transition.Frames/2;
                    Transition.Sigma = Transition.Frames/5;
                    Transition.Increments = normpdf(1:Transition.Frames,Transition.Mu,Transition.Sigma);
                    
                case 4 %============== SINUSOID
                    
            end


            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            while 1
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    parseCommand(commandIn);                % parse command
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn, 'R_dotsAndPicOff')
                    break
                end
 
                %==================== INITIALIZE VARIABLES ================
                if Frame == 1 && Loop == 0
                    Intermittent.SwitchOn = GetSecs;            % Set switch time to start time
                    Pic.CurrentContrast = Pic.MinContrast;      % Start with pic off
                    Dots.CurrentContrast = Dots.MaxContrast;    % Start with dots on
                    Pic.ContrastRange = Pic.MaxContrast-Pic.MinContrast;
                    Dots.ContrastRange = Dots.MaxContrast-Dots.MinContrast;
                else
                    %================ YOKED TRANSITION ====================
                    if Frame <= Transition.Frames
                        if Stim == 1                                    % DOTS ON
                            if Pic.CurrentContrast > Pic.MinContrast
                                Pic.CurrentContrast = Pic.CurrentContrast- (Transition.Increments(Frame)*Pic.ContrastRange);
                            end
                            if Dots.CurrentContrast < Dots.MaxContrast
                                Dots.CurrentContrast = Dots.CurrentContrast+ (Transition.Increments(Frame)*Dots.ContrastRange);
                            end
                        elseif Stim == 2                                % IMAGE ON
                            if Pic.CurrentContrast < Pic.MaxContrast
                                Pic.CurrentContrast = Pic.CurrentContrast+ (Transition.Increments(Frame)*Pic.ContrastRange);
                            end
                            if Dots.CurrentContrast > Dots.MinContrast
                                Dots.CurrentContrast = Dots.CurrentContrast- (Transition.Increments(Frame)*Dots.ContrastRange);
                            end
                        end
                    end
                end
                
                %===================== DRAW STIMULI =======================
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                if Flicker.Eye1 == 1 && Intermittent.StimOn == 1
                    Screen('DrawTexture', Display.win, DotsTexture(Frame), [], Pic.Rect, Dots.DirMean, [], Dots.CurrentContrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);                  % Apply Gaussian aperture mask
                    Screen('FillOval', Display.win, [Photodiode.OnColour, Dots.CurrentContrast*(1/Dots.MaxContrast)*255], Photodiode.Rect{Eye+1});
                end
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, abs(floor(Eye-0.5)));         % Select which eye to draw to
                Screen('DrawTexture', Display.win, BackgroundTexture);
                if Flicker.Eye2 == 1 && Intermittent.StimOn == 1
                    Screen('DrawTexture', Display.win, ImageTexture, [], Pic.Rect, Pic.Rotation, [], Pic.CurrentContrast);
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);                  % Apply Gaussian aperture mask
                    Screen('FillOval', Display.win, [Photodiode.OnColour Pic.CurrentContrast*(1/Pic.MaxContrast)*255], Photodiode.Rect{abs(floor(Eye-0.5))+1});
                end
                if Fix.On == 1 || Intermittent.StimOn == 0
                    for E = 1:2
                        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E-1);
                        Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect, Fix.Rotation(E)); 
                    end
                end
                [VBL FrameOn] = Screen('Flip', Display.win);


                %=================== FREQUENCY TAGGING
                if Flicker.On == 1
                    if mod(ceil(FrameTotal/Flicker.DurL),2)
                        Flicker.Eye1 = 1;                   % Flicker on
                    else
                        Flicker.Eye1 = 0;                   % Flicker off
                    end
                end
                
                %================== INTERMITTENT PRESENTATION
                if Intermittent.On == 1
                    if Intermittent.StimOn == 1 && GetSecs >= Intermittent.SwitchOn+Intermittent.OnDur
                        Intermittent.StimOn = 0;
                        Intermittent.SwitchOff = GetSecs;
                    elseif Intermittent.StimOn == 0 && GetSecs >= Intermittent.SwitchOff+Intermittent.OffDur
                        Intermittent.StimOn = 1;
                        Intermittent.SwitchOn = GetSecs;
                    end
                end
                
                %================== ANIMATION LOOP
                Frame = Frame+1;
                FrameTotal = (Loop*Dots.NrFrames)+Frame;
                if Frame == Dots.NrFrames
                    Frame = 1;
                    Loop = Loop+1;
                end
                [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                if keyIsDown && keyCode(Key.Exit)
                    break;
                end
            end
%             Screen('Close', DotsTexture);
%             Screen('Close', ImageTexture);
         	R_dotsAndPicOff;
        end
        

        %====================== SET IMAGE PARAMETERS ======================
        function R_setPic(varargin)
            params=varargin{1,1};
            if isfield(Pic, 'ScaleX')                                   % If Scale fields already exist
                if params(2) ~= Pic.ScaleX || params(3) ~= Pic.ScaleY   % If requested image size has changed...
                    NewBackground = 1;
                else
                    NewBackground = 0;
                end
            else
                NewBackground = 1;
            end
            Pic.Number = params(1);
            Pic.ScaleX = params(2);
            Pic.ScaleY = params(3);
            Pic.PosX = params(4);
            Pic.PosY = params(5);
            Pic.MaxContrast = params(6);
%             BackgroundOn = params(7);
            Pic.Category = params(8);
            
            ImageRadius = round((Pic.ScaleX/2)*Display.Pixels_per_deg(1));
            AllImageTypes = regexp(genpath(ImageDir),['[^;]*'],'match');                % Get list of image subcategory folders
            CategoryDir = AllImageTypes{Pic.Category+1};                                % Get path of requested image category
            [x ImageCategory] = fileparts(CategoryDir);                                 % Get text string for category
            cd(CategoryDir);                                                            % Change to that folder
            Images = dir('*.jpg');                                                      % Get all jpg images
            [Img ImgCmap] = imread(fullfile(CategoryDir, Images(Pic.Number).name));
            ImgDim = size(Img);                                                         % Get image dimensions
            if min(ImgDim([1,2])) ~= 2*ImageRadius                                    	% if smallest image dimesnion is not requested size...
                scale = 2*ImageRadius/min(ImgDim([1 2]));                               % Resize image so smallest dimension fits
                Img = imresize(Img, scale);
                ImgDim = size(Img);                                                     % Get new Img image dimensions
            end         
            if max(ImgDim([1,2])) > 2*ImageRadius                                   	% If largest image dimension is too large...
                Crop = (max(ImgDim([1,2]))-(2*ImageRadius))/2;
                if find(ImgDim==max(ImgDim([1,2])))==1
                    Img(end-Crop:end,:,:) = [];
                    Img(1:Crop,:,:) = [];
                elseif find(ImgDim==max(ImgDim([1,2])))==2
                    Img(:,end-Crop:end,:) = [];
                    Img(:,1:Crop,:) = [];
                end
            end
            Pic.ImageRect = [0,0,size(Img,2),size(Img,1)];
            Pic.Rect = CenterRect(Pic.ImageRect,Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            ImageTexture = Screen('MakeTexture', Display.win, double(Img));             % Convert image to PTB texture
            
            %========================= CREATE ALPHA CHANNEL MASK ======================
            TextureWindow = [2*ImageRadius, 2*ImageRadius];                     % Set size of the stimulus before mask is applied (degrees)
            Dots.Window = TextureWindow;
            Stim.Background = Display.Background(1);
            Mask.Dim = TextureWindow+4;                                         % Mask will be 2 pixels larger than the stimulus window
            Mask.ApRadius = (min(Mask.Dim)/2)-8;
            Mask.Colour = Display.Background(1);
            Mask.Edge = 0;
            MaskTex = GenerateAlphaMask(Mask, Display);
            Mask.Rect = [0 0 Mask.Dim];                                      	% Mask size
            Mask.DestRect = CenterRect(Mask.Rect, Display.Rect);            	% Destination for mask is centred in screen window
            Mask.DestRect = Mask.DestRect + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            Stim.Background = Display.Background;
          	Stim.Window = [0 0 Mask.Dim];
%             if BackgroundOn == 1 
%                 if NewBackground == 1
%                     BackgroundTexture = BackgroundSquares(Display, Stim);        	% Generate squares background
%              	end
%             elseif BackgroundOn == 0
                BackgroundTexture = Screen('MakeTexture', Display.win, ones(Display.Rect([4,3]))*Display.Background(1));
%             end
        end
        
        %========================== SET DOT PARAMETERS ====================
        function R_setDots(varargin)
            params=varargin{1,1};
            Dots.Num = params(1);
            Dots.Colour = params([2,3,4]);
            Dots.MaxContrast = params(5);
          	if Penalty.On == 1
                Dots.RandColour = 1;        % Set to 1 for random coloured dots for training
            else
                Dots.RandColour = 0;        % Set to 1 for random coloured dots for training
            end
        end
        
        function R_setDotsMove(varargin)
            params=varargin{1,1};
            Dots.Signal = params(1);
            Dots.SpeedMean = params(2);
            Dots.SpeedStd = params(3);
            Dots.LifeMean = params(4);
            Dots.LifeStd = params(5);
            Dots.Lifetime = Dots.LifeMean;
        end
        
        function R_setDotsDir(varargin)
            params=varargin{1,1};
            Dots.DirMean = params(1)-90;
            Dots.DirStd = params(2);
            Dots.SizeMean = params(3)*Display.Pixels_per_deg(1);
            Dots.SizeStd = params(4)*Display.Pixels_per_deg(1);
            
         	%================ GENERATE TRANSPARENT MOTION TEXTURES ============
            TextureWindow = [300, 300]; %<<<<<< FUDGE FOR MISSING VARIABLE!!!!
            Dots.Window = TextureWindow;
            Dots.Velocity = Dots.SpeedMean*Display.Pixels_per_deg(1);
            Dots.Background = Display.Background;
            Dots.DrawAngle = Dots.DirMean;
            Dots.Size = ([Dots.SizeMean, Dots.SizeStd, 1]);
            DotsTexture = GenerateRDK(Dots, Display);
            Dots.NrFrames = numel(DotsTexture);
        end
        
        
        %=========================== CLEAR TEXTURE ========================
        function R_dotsAndPicOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, BackgroundTexture); 
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win, [], 1);                                       % Clear stimulus  
        end
        
        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
        end

    end


%% ========================== apm_binorivSLOF =============================
% Presents different static images to the two eyes, which intermittently
% move in opposite directions, generating short-latency occular following
% eye movements in the direction of the perceptually dominant stimulus's
% motion.
%
% OPTIONAL SETTINGS:
%       Intermittent:   Presents stimuli intermittently for perceptual stabilization
%       Flicker:        Flickers stimuli at different rates in each eye for frequency tagging
%       OccularSwitch: 	Switches stimuli rapidly between the two eyes for stimulus rivalry
%       SLOF:           Brief sudden movements of stimuli to induce short-latency occular following
%
%==========================================================================

    function BinorivHandle = apm_binorivSLOF()

        BinorivHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_dotsAndPicOn;...
                @R_setPic;...
                @R_setDots;...
                @R_setDotsMove;...
                @R_setDotsDir;...
                @R_dotsAndPicOff;...
                @R_getBlock;...
            };
        
        persistent Fix;
        persistent Dots;
        persistent Mask;
        persistent Frame;
        persistent Loop;
        persistent DotsTexture;
        persistent ImageTexture;
        persistent BackgroundTexture;
        persistent MaskTex;
        persistent Flicker;
        persistent Intermittent;
        persistent BackgroundOn;
        persistent Probe;

        %======================= SET STIMULUS PARAMETERS ======================
        BackgroundOn = 0;
        Fix.On = 0;                                     % Keep fixation on during stimulus presentation?
        Fix.Type = 0;
        Fix.Eye = 2;
        
        Pic.Rotation = 0;
        Dots.Contrast = 1;
        Frame = 1;
        FrameTotal = 1;
        Loop = 0;
        StimType = 2;                   % 1 = Gratings; 2 = Images
        
        %============= FLICKER FREQUENCY TAGGING?
        Flicker.On = 0;                 % Add frequency tagging to stimuli?
        Flicker.RateHzL = 6;            % Set flicker rate (Hz)
        Flicker.RateHzR = 15;     
        Flicker.DurL = (1/Flicker.RateHzL)*Display.RefreshRate; % Convert to frames
        Flicker.DurR = (1/Flicker.RateHzR)*Display.RefreshRate; % Convert to frames
        Flicker.Eye1 = 1;
        Flicker.Eye2 = 1;
        
        %============= INTERMITTENT PRESENTATION?
        Intermittent.On = 0;                
        Intermittent.OnDur = 1;             % Set duration of stimulus on period (s)
        Intermittent.OffDur = 1;            % Set duration of blank period (s)
        Intermittent.StimOn = 1;            % Start trial with stimulus 'on'
        
        %============= SHORT LATENCY OCCULAR FOLLOWING?
        TrialDuration = 300;
        Probe.On = 1;                                    % 
        Probe.Range = [1 2];                             % Set time range for frequency of motion probe (seconds)
        Probe.Duration = 0.2;                            % Set duration of motion probe (seconds)
        Probe.Speed = 50*Display.Pixels_per_deg(1);      % Set speed of probe motion (degrees/second)
        Probe.Onsets = cumsum((rand(1,TrialDuration/min(Probe.Range))*diff(Probe.Range))+Probe.Range(1));
        Probe.Onsets(Probe.Onsets>TrialDuration) = [];
        Probe.Onsets(end+1) = TrialDuration+1;
        Probe.ShiftPerFrame = Probe.Speed/Display.RefreshRate;    % Pixels per frame
        NrFrames = Probe.Duration*Display.RefreshRate;

     	Stim.Window = [0 0 500 500];
        Stim.Background = Display.Background;
        
%         BackgroundTexture = BackgroundSquares(Display, Stim); 
%         if Stereobackground == 0
%             Screen('FillRect', Display.win, Display.Background, Display.Rect);
%             Screen('Flip', Display.win);
%         end
        
     	%======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            Fix = fixOn(Fix);
          	[VBL FixOn] = Screen('Flip', Display.win);
        end
        
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            %params=sscanf(varargin{1,1},'%f');
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
         %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
%                 Screen('DrawTexture', Display.win, BackgroundTexture);
                Screen('FillRect', Display.win, Display.Background);
            end 
            [VBL FixOff] = Screen('Flip', Display.win);                                         % Clear stimulus  
        end
        
        %======================== PRESENT STIMULI =========================
        function R_dotsAndPicOn(Eye)
            Frame = 1;
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            
            for E = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E-1);     	% Select which eye to draw to
                Screen('DrawTexture', Display.win, BackgroundTexture);
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect); 
            end
            [VBL FixOn] = Screen('Flip', Display.win); 
            
            xoffsetL = 0;
            xoffsetR = 0;
            while 1
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    parseCommand(commandIn);                % parse command
                    pnet(con,'printf', reply);              % return control to the socket
                end
%                 if strcmp(commandIn,'R_dotsandPicOff')     % Exits loop if trial is quit or ended
%                     break
%                 end
                 
                %===================== DRAW STIMULI =======================
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                if Flicker.Eye1 == 1 && Intermittent.StimOn == 1
                    
                    if StimType == 1
                        xoffsetL = mod((Frame-1)*shiftperframe,pixelsPerCycle);
                        xoffsetR = mod((NrFrames-(DotDirection*Frame))*shiftperframe,pixelsPerCycle);
                    elseif StimType == 2
                        xoffsetL = xoffsetL + ((Frame-1)*Probe.ShiftPerFrame);
                        xoffsetR = xoffsetR - ((Frame-1)*Probe.ShiftPerFrame);
                    end
                    Probe(Eye).srcRect =[xoffsetL 0 xoffsetL + ImageDim(1) ImageDim(1)];
                    Probe(abs(floor(Eye-0.5))).srcRect =[xoffsetR 0 xoffsetR + ImageDim(1) ImageDim(1)];
                    
                    Screen('DrawTexture', Display.win, TextureR, Probe(Eye).srcRect, [], Probe(Eye).RotationAngle+Probe(Eye).OffsetAngle, [], CatchAlphaL);
                    
                    
%                     Screen('DrawTexture', Display.win, DotsTexture(Frame), [], Pic.Rect, Dots.DirMean, [], Dots.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);                  % Apply Gaussian aperture mask
                end
                if Fix.On == 1
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect); 
                end
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, abs(floor(Eye-0.5)));     % Select which eye to draw to
                Screen('DrawTexture', Display.win, BackgroundTexture);
                if Flicker.Eye2 == 1 && Intermittent.StimOn == 1
                    Screen('DrawTexture', Display.win, ImageTexture, [], Pic.Rect, Pic.Rotation, [], Pic.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);                  % Apply Gaussian aperture mask
                end
                if Fix.On == 1
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect); 
                end
                [VBL FrameOn] = Screen('Flip', Display.win);
                if Frame == 1 && Loop == 0
                    Intermittent.SwitchOn = FrameOn;
                end

                %=================== FREQUENCY TAGGING
                if Flicker.On == 1
                    if mod(ceil(FrameTotal/Flicker.DurL),2)
                        Flicker.Eye1 = 1;                   % Flicker on
                    else
                        Flicker.Eye1 = 0;                   % Flicker off
                    end
                    if mod(ceil(FrameTotal/Flicker.DurR),2)
                        Flicker.Eye2 = 1;                   % Flicker on
                    else
                        Flicker.Eye2 = 0;                   % Flicker off
                    end
                end
                
                %================== INTERMITTENT PRESENTATION
                if Intermittent.On == 1
                    if Intermittent.StimOn == 1 && GetSecs >= Intermittent.SwitchOn+Intermittent.OnDur
                        Intermittent.StimOn = 0;
                        Intermittent.SwitchOff = GetSecs;
                    elseif Intermittent.StimOn == 0 && GetSecs >= Intermittent.SwitchOff+Intermittent.OffDur
                        Intermittent.StimOn = 1;
                        Intermittent.SwitchOn = GetSecs;
                    end
                end
                
                %================== MOTION PROBE
                if Probe.On == 1
                    if GetSecs-StartTime >= Probe.Onsets(NextProbe) && MotionOn == 0   % If the next probe is due...
                        MotionOn = 1;
                        MotionOnTime = GetSecs;
                        NextProbe = NextProbe+1;
                    elseif GetSecs >= MotionOnTime+Probe.Duration && MotionOn == 1     % If the probe duration has elapsed...
                        MotionOn = 0;
                        MotionOnTime = inf;
                        Frame = 1;
                        xoffsetL = 0; xoffsetR = 0;
                    end
                end
                
                %================== ANIMATION LOOP
                Frame = Frame+1;
                FrameTotal = (Loop*Dots.NrFrames)+Frame;
                if Frame == Dots.NrFrames
                    Frame = 1;
                    Loop = Loop+1;
                end
                [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                if keyIsDown && keyCode(Key.Exit)
                    break;
                end
            end
            
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, BackgroundTexture);
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus
            Screen('Close', DotsTexture);
            Screen('Close', ImageTexture);
        end
        
        %====================== SET IMAGE PARAMETERS ======================
        function R_setPic(varargin)
            params=varargin{1,1};
            Pic.Number = params(1);
            Pic.ScaleX = params(2);
            Pic.ScaleY = params(3);
            Pic.PosX = params(4);
            Pic.PosY = params(5);
            Pic.Contrast = params(6);
            
            ImageRadius = round((Pic.ScaleX/2)*Display.Pixels_per_deg(1));
            cd(ImageDir);
            Images = dir('*.jpg');
            [Img ImgCmap] = imread(fullfile(ImageDir, Images(Pic.Number).name));
            ImgDim = size(Img);                                                         % Get image dimensions
            if min(ImgDim([1,2])) ~= 2*ImageRadius                                    	% if smallest image dimesnion is not requested size...
                scale = 2*ImageRadius/min(ImgDim([1 2]));                               % Resize image so smallest dimension fits
                Img = imresize(Img, scale);
                ImgDim = size(Img);                                                     % Get new Img image dimensions
            end         
            if max(ImgDim([1,2])) > 2*ImageRadius                                   	% If largest image dimension is too large...
                Crop = (max(ImgDim([1,2]))-(2*ImageRadius))/2;
                if find(ImgDim==max(ImgDim([1,2])))==1
                    Img(end-Crop:end,:,:) = [];
                    Img(1:Crop,:,:) = [];
                elseif find(ImgDim==max(ImgDim([1,2])))==2
                    Img(:,end-Crop:end,:) = [];
                    Img(:,1:Crop,:) = [];
                end
            end
            Pic.ImageRect = [0,0,size(Img,2),size(Img,1)];
            Pic.Rect = CenterRect(Pic.ImageRect,Display.Rect) + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            ImageTexture = Screen('MakeTexture', Display.win, double(Img));             % Convert image to PTB texture
            
            %========================= CREATE ALPHA CHANNEL MASK ======================
            TextureWindow = [2*ImageRadius, 2*ImageRadius];                     % Set size of the stimulus before mask is applied (degrees)
            Dots.Window = TextureWindow;
            Stim.Background = Display.Background(1);
            Mask.Dim = TextureWindow+4;                                         % Mask will be 2 pixels larger than the stimulus window
            Mask.ApRadius = (min(Mask.Dim)/2)-8;
            Mask.Colour = Display.Background(1);
            Mask.Edge = 0;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
            MaskTex = GenerateAlphaMask(Mask, Display);
            Mask.Rect = [0 0 Mask.Dim];                                      	% Mask size
            Mask.DestRect = CenterRect(Mask.Rect, Display.Rect);            	% Destination for mask is centred in screen window
            Mask.DestRect = Mask.DestRect + [Pic.PosX,Pic.PosY,Pic.PosX,Pic.PosY];
            Stim.Background = Display.Background;
          	Stim.Window = [0 0 Mask.Dim];
            if BackgroundOn == 1
                BackgroundTexture = BackgroundSquares(Display, Stim);        	% Generate squares background
            elseif BackgroundOn == 0
                BackgroundTexture = Screen('MakeTexture', Display.win, ones(Display.Rect([4,3]))*Display.Background(1));
            end
        end
        
        %========================== SET DOT PARAMETERS ====================
        function R_setDots(varargin)
            params=varargin{1,1};
            Dots.Num = params(1);
            Dots.Colour = params([2,3,4])*255;
        end
        
        function R_setDotsMove(varargin)
            params=varargin{1,1};
            Dots.Signal = params(1);
            Dots.SpeedMean = params(2);
            Dots.SpeedStd = params(3);
            Dots.LifeMean = params(4);
            Dots.LifeStd = params(5);
            Dots.Lifetime = Dots.LifeMean;
        end
        
        function R_setDotsDir(varargin)
            params=varargin{1,1};
            Dots.DirMean = params(1)-90;
            Dots.DirStd = params(2);
            Dots.SizeMean = params(3)*Display.Pixels_per_deg(1);
            Dots.SizeStd = params(4)*Display.Pixels_per_deg(1);
            
         	%================ GENERATE TRANSPARENT MOTION TEXTURES ============
            Dots.Velocity = Dots.SpeedMean*Display.Pixels_per_deg(1);
            Dots.Background = Display.Background;
            Dots.DrawAngle = Dots.DirMean;
            Dots.Size = ([-Dots.SizeStd, Dots.SizeStd]+Dots.SizeMean);
            DotsTexture = GenerateRDK(Dots, Display);
            Dots.NrFrames = numel(DotsTexture);
        end
        
        %=========================== CLEAR TEXTURE ========================
        function R_dotsAndPicOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('DrawTexture', Display.win, BackgroundTexture);
                Screen('DrawTexture', Display.win, Fix.Texture); 
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus  
        end
        
        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
        end

    end


%% ============================== apm_movie ===============================
% Loads and presents regular or 3D movies from the Movie.Dir directory.
%
%==========================================================================

    function MovieHandle = apm_movie()

        MovieHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_CalibrateMovie;...
                @R_videoOn;...
                @R_videoOff;...
                @R_getBlock;...
                @R_LoadMovie;...
            };
        
        persistent Fix;
        persistent mov;
        persistent Movie
        
%         Screen('BlendFunction', Display.win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);     % Change alpha blend mode (<<< FIX for new PTB contrast issue added by APM 18/10/2015)
        
        %====================== SET DEFAULT PARAMETERS ====================
        Fix.Type            = 0;  
        Fix.FixAlwaysOn     = 0;
        Fix.Eye             = 2;                                
        Movie.Contrast      = 1;                    % Play at full contrast
        Movie.Mirror        = 1;                   	% Mirror invert video to display correctly? (Warning: increases processing!)
        Movie.StartTime     = 0;                  	% Specify start time (seconds)
        Movie.Eye           = 2;                 	% Movie presentation defaults to binocular
        Movie.Loop          = 1;                 	% Loop movie after last frame?
%         Movie.Background = 128;                	% Set movie background color (RGB)
        Movie.TimeInFrames = 0;                     % 0 = Specify time in seconds; 1 = specify time in frames
        
        %============================= OPEN MOVIE =========================
        function R_LoadMovie(varargin) 
            params=varargin{1,1};
            Movie.Number            = params(1);         	% Specify which movie file to play
            Movie.PlaybackTime      = params(2);            % Set playback duration (seconds)
            Movie.Dims              = params([3,4]);      	% Movie dimensions (if Movie.Fullscreen = 0)
            Movie.FullScreen        = params(5);            % use full screen?
            Movie.Rate              = params(6)/100;       	% Playback speed (proportion: -1 = reverse, <1 = slow, 1 = normal, >1 = fast)
            if Audio.On == 1
                Movie.Volume        = params(7)/100;    	% Audio track volume (1 = full, 0 = mute)
            elseif Audio.On == 0
                Movie.Volume        = 0;
            end
            Movie.StaticFrames      = params(8);            % Present movie frames as sequential static images
            Movie.StaticFPS         = params(9);            % Set frame rate for static presentation
            Movie.Type              = params(10);          	% Set movie type: 1 = Monkey Thieves; 2 = Russ movies; 4 = 3D movies
            Movie.MaintainAR        = params(11);           % Adjust aspect ratio of movie to fit movie.Dims?
            Movie.Background        = params(12);
            Movie.Show3D            = 0;
            switch Movie.Type 
                case 1
                    Movie.Dir = 'C:\Users\lab\Desktop\Movies\MonkeyThieves';
%                     Movie.Dir = 'C:\Documents and Settings\lab\Desktop\Movies\MonkeyThieves';
                case 2
                    Movie.Dir = 'C:\Users\lab\Desktop\Movies\RussMoviesWMV';
%                     Movie.Dir = 'C:\Documents and Settings\lab\Desktop\Movies\BriansMovies';
                case 3
                    Movie.Dir = 'C:\Users\lab\Desktop\Movies\LubekWMV';
                case 4
                    Movie.Dir = 'C:\Users\lab\Desktop\Movies\3DMovies';
                    Movie.Show3D = 1;
            end
            if Movie.Type == 2 || Movie.Type == 3
                if Movie.Mirror == 0
                    MacaqueMovie = dir(fullfile((Movie.Dir),['Movie',num2str(Movie.Number),'.wmv']));
                    Movie.Filename = fullfile(Movie.Dir, MacaqueMovie(1).name);
                elseif Movie.Mirror == 1
                    MacaqueMovie = dir(fullfile((Movie.Dir),['Movie',num2str(Movie.Number),'_mirrored.mp4']));
                    Movie.Filename = fullfile(Movie.Dir, MacaqueMovie(1).name);
                end
            elseif Movie.Type ~= 2
                MacaqueMovies = dir(Movie.Dir);
            	Movie.Number = Randi(numel(MacaqueMovies)-2)+2;
            	Movie.Filename = fullfile(Movie.Dir, MacaqueMovies(Movie.Number).name);
            end
            
            %=============== LOAD SCENE TRASNITION DATA
            TransitionsFile = dir([Movie.Filename(1:strfind(Movie.Filename,'.')-1), '_*.mat']);
%             TransitionsFile = wildcardsearch(Movie.Dir, [Movie.Filename(1:strfind(Movie.Filename,'.')-1), '_*.mat']);
            if isempty(TransitionsFile)
                Movie.TransitionFrames = [];
            else
                Trans = load(fullfile(Movie.Dir, TransitionsFile(1).name));
                Movie.TransitionFrames = Trans.ThreshFrames;
                Movie.TransitionFrames(end+1) = inf;
            end
            
            Movie.Format = Movie.Filename([end-5, end-4]);          % Get movie file format
            Movie.PixelFormat = 4;                                  % 3 = RGB; 4 = RGBA (default)
%             [mov, Movie.duration, Movie.fps, Movie.width, Movie.height, Movie.count, Movie.AR]= Screen('OpenMovie', Display.win, Movie.Filename, [], [], [], Movie.PixelFormat); 
            [mov, Movie.duration, Movie.fps, Movie.width, Movie.height, Movie.count, Movie.AR]= Screen('OpenMovie', Display.win, Movie.Filename);
            Movie.SourceRect{2} = [0 0 Movie.width, Movie.height];
            if Movie.Type == 2
                Movie.SourceRect{2} = [5 0 Movie.width-5, Movie.height];
            end
            Screen('SetmovieTimeIndex',mov,Movie.StartTime,Movie.TimeInFrames);
            Screen('PlayMovie',mov,Movie.Rate, Movie.Loop, Movie.Volume);

            if Movie.MaintainAR == 0
                if Movie.FullScreen == 1
                    Movie.DestRect = Display.Rect;
                elseif Movie.FullScreen == 0
                    Movie.DestRect = [0 0 Movie.Dims]*Display.Pixels_per_deg(1);
                end
            elseif Movie.MaintainAR == 1
                if Movie.FullScreen == 1
                    Movie.WidthDeg = Display.Rect(3);
                else
                    Movie.WidthDeg = Movie.Dims(1)*Display.Pixels_per_deg(1);
                end
                Movie.DestRect = (Movie.SourceRect{2}/Movie.width)*Movie.WidthDeg;
            end
            if ~isempty(find(Movie.DestRect > Display.Rect))
                Movie.DestRect = Movie.DestRect*min(Display.Rect([3, 4])./Movie.DestRect([3, 4]));
                fprintf('Requested movie size > screen size! Defaulting to maximum size.\n');
            end
            Movie.DestRect = CenterRect(Movie.DestRect, Display.Rect);
            if Movie.Show3D == 1
                if strcmpi(Movie.Format, 'LR')          % Horizontal split screen
                    Movie.SourceRect{1} = Movie.SourceRect{2}./[1 1 2 1];
                    Movie.SourceRect{2} = Movie.SourceRect{1}+[Movie.SourceRect{1}(3),0, Movie.SourceRect{1}(3),0];     
                elseif strcmpi(Movie.Format, 'RL')          
                    Movie.SourceRect{2} = Movie.SourceRect{2}./[1 1 2 1];
                    Movie.SourceRect{1} = Movie.SourceRect{2}+[Movie.SourceRect{2}(3),0, Movie.SourceRect{2}(3),0];  
                elseif strcmpi(Movie.Format, 'TB')      % Vertical split screen
                    Movie.SourceRect{1} = Movie.SourceRect{2}./[1 1 1 2];
                    Movie.SourceRect{2} = Movie.SourceRect{1}+[0,Movie.SourceRect{1}(4),0, Movie.SourceRect{1}(4)];
                else
                    fprintf('\nERROR: 3D movie format must be specified in filename!\n');
                end
            else
                Movie.SourceRect{1} = Movie.SourceRect{2};
            end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Movie.Background);           % Set background to black
            end
           
        end
        
        %======================= RECALIBRATE EYE POS ======================
        function R_CalibrateMovie(varargin)
            fprintf('\nRunning 9-point calibration...\n');
            params=varargin{1,1};
            MaxXpos = params(1)/2;
            MaxYpos = params(2)/2;
            PointDur = 0.5;                     % Duration of each target position (seconds)
            AllPos = [0, -1, -1, -1, 0, 1, 1, 1, 0;...
                      1, 1, 0, -1, -1, -1, 0, 1, 0]';
        	AllPos = AllPos.*repmat([MaxXpos, MaxYpos],[size(AllPos, 1),1]);
            Fix.Size = 0.5*Display.Pixels_per_deg(1);
            Fix.Colour = [0 0 0];
            Fix.Eye = 2;
            Fix.Type = 0;
            FixOn = GetSecs;
            for p = 1:size(AllPos, 1)
                Fix.Offset = AllPos(p,:)*Display.Pixels_per_deg(1);     % Set fixation poisition
                Fix = fixOn(Fix);                                       % Darw fixation
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                    Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
                end
            	while GetSecs < FixOn+PointDur
                    
                end
                [VBL FixOn] = Screen('Flip', Display.win);
            end
        end
        
        %========================= PRESENT MOVIE ==========================
        function R_videoOn(varargin)
            params = varargin{1,1};
            Eye = params(1);
            
            Screen('SetmovieTimeIndex',mov,Movie.StartTime,Movie.TimeInFrames);
            Screen('PlayMovie',mov, Movie.Rate, Movie.Loop, Movie.Volume);
            
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            EndMovie = 0;
            t = 1;
            AllFrameOnsets = [];
            while EndMovie == 0
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    parseCommand(commandIn);                % parse command
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn,'R_videoOff')         	% Exits loop if movie is being turned off.     
                    EndMovie =1;
                    break
                end
                [MovieTex, Timeindex] = Screen('GetMovieImage', Display.win, mov, 1);
                FrameIndex = Timeindex*Movie.fps;
                if MovieTex < 0
                    break;
                end
%                 if Movie.Mirror == 1
%                     SourceRectScale = [1 1 2 1];
%                     Eye = 1;
%                     Screen('DrawTexture', Display.win, MovieTex, Movie.SourceRect{Eye}, Movie.SourceRect{1});
%                     array2Flip = Screen('GetImage', Display.win, Movie.SourceRect{Eye}.*SourceRectScale, 'backBuffer');
%                     FlippedArray = array2Flip(:,end:-1:1,:);
%                     MovieTex = Screen('MakeTexture', Display.win, FlippedArray); 
%                     Screen('FillRect',Display.win, Movie.Background);
%                 else
                    SourceRectScale = [1 1 1 1];
%                 end
                if ~isempty(Movie.TransitionFrames) 
                    if FrameIndex >= Movie.TransitionFrames(t)      % Mark scene transitions with slight change in photodiode luminance
                        PhotodiodeColour = [64 64 64];
                        t = t +1;
                    else
                        PhotodiodeColour = Photodiode.OnColour;
                    end
                else
                    PhotodiodeColour = Photodiode.OnColour;
                end
                
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                    Screen('FillRect',Display.win, Movie.Background);
                    if MovieTex > 0
                        Screen('DrawTexture', Display.win, MovieTex, (Movie.SourceRect{Eye}.*SourceRectScale), Movie.DestRect);
                        Screen('FillOval', Display.win, PhotodiodeColour, Photodiode.Rect{Eye});
                        array2Flip = Screen('GetImage', Display.win, Movie.DestRect, 'backBuffer');
                    end
                    if Fix.FixAlwaysOn == 1
                        Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect, Fix.Rotation(Eye)); 
                    end
                end
                Screen('DrawingFinished', Display.win);
                [VBL FrameOn] = Screen('Flip', Display.win);
            	Screen('Close', MovieTex);
                if Movie.StaticFrames == 1
                    while GetSecs < FrameOn+(1/Movie.StaticFPS)
                        
                    end
                end
                [keyIsDown,secs,keyCode] = KbCheck;                     % Check keyboard for 'escape' press        
                if keyIsDown && keyCode(Key.Exit) == 1                   % Press Esc for abort
                    EndMovie = 1;
                end
            end
            R_videoOff(mov);
        end
        
        %=========================== CLEAR TEXTURE ========================
        function R_videoOff(mov)
            if nargin == 1
                Screen('PlayMovie',mov, 0, Movie.Loop, Movie.Volume);   % Stop playback
                Screen('CloseMovie', mov);
            end
            
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Movie.Background); 
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus  
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
        %======================== PRESENT FIXATION ============================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Size = params(4)*Display.Pixels_per_deg(1);
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Eye = params(3);
            Fix.FixAlwaysOn = params(5);
            Fix = fixOn(Fix);
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
         	[VBL FixOn] = Screen('Flip', Display.win);
        end

        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Movie.Background, Display.Rect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win, [], 1);                                         % Clear stimulus  
        end
        
        function R_getBlock(n)
            Trial_num = 0;
            Block = 0;
        end
        
    end


%% ============================ apm_RFMapping =============================
% Present moving dot or grating patches across the visual field for
% receptive field mapping of LIP neurons.
%
% Optimal stimulus diameter is calculated based on the requested eccentricity,
% using the linear regression equations from:
%   Hamed et al (2001) for LIP neurons [y = 0.9x + 5.842]
%   Petersen, Robinson & Keys (1985) for Pulvinar
%
% REFERENCES
%   Hamed SB, Duhamel JR, Bremmer F & Graf W (2001). Representation of the
%       visual field in the lateral intraparietal area of macaque monkeys: a
%       quantitative receptive field analysis. Exp Brain Res, 140: 127-144;
%   Petersen SE, Robinson DL, Keys W (1985). Pulvinar nuclei of the behaving 
%       rhesus monkey: visual responses and their modulation. J
%       Neurophysiol, 54: 867-886.
%
%==========================================================================

    function RFHandle = apm_RFMappingNew()

        RFHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_stimOn;...
                @R_stimOff;...
                @R_makeStim;...
                @R_getBlock;...
            };
        
        persistent Dots;
        persistent Grating;
        persistent Stim;
        persistent StimTexture;
        persistent Fix;
        persistent Mask;
        persistent MaskTex;
        
        %======================= SET STIMULUS PARAMETERS ======================
        Fix.Type    = 1;                                        % 1 = cross; 
        Fix.Eye     = 2;
        Stim.Type   = 4;                                        % 1 = moving dots, 2 = Static Gabor, 3 = Static blob, 4 = Dynamic noise
        Stim.ROI    = 'PI';                                     % Specify ROI being tested
        Stim.Contrast = 1;                                  
        if Stim.Type < 3
            Stim.NoDirections   = 4;
            Stim.Directions     = 0:(360/Stim.NoDirections):(360-(360/Stim.NoDirections));
            Stim.Eccentricities = [0 4 8 10];               % Set eccentiricities of stimulus center (deg)
            Stim.NoPolarAngles  = [1 8 8 8];
        elseif Stim.Type >= 3
            Stim.NoDirections   = 1;
            Stim.Directions     = 0;
            Stim.Eccentricities = [0.8 2.2 4 6.5 10];       % Set eccentiricities of stimulus centers (deg)
            Stim.NoPolarAngles  = [4 8 12 14 16];          	% Set the number of polar angles for each eccentricity
        end
        Stim.NoConditions = sum(Stim.NoPolarAngles)*Stim.NoDirections;
        Stim.Params = zeros(1,8);
        
        %============= Calculate optimal stimulus diameter for specified eccentricity
        switch Stim.ROI
            case 'LIP'          % based on Hamed et al (2001) 
                b = 0.9;        % regression fit slope
                a = 5.842;      % Y intercept (degrees)
            case 'PI'           % based on Petersen, Robinson & Keys (1985)
                b = 0.3235;
                a = 1;
            case 'PL'
                b = 0.5357;
                a = 0.5;
            case 'Pdm'
                b = 0.4412;
                a = 3.5;
            case 'LGN'
                b = 0.2;
                a = 0.2;
        end
        Stim.Diameters = (b*Stim.Eccentricities + a)*Display.Pixels_per_deg(1);
            
        Stim.Variables = zeros(Stim.NoConditions,3);        % N x 3 Matrix, 1 = eccentricity, 2 = polar, 3 = direction
        for e = 1:numel(Stim.Eccentricities)
         	Stim.PolarAngles{e} = 0:(360/Stim.NoPolarAngles(e)):(360-(360/Stim.NoPolarAngles(e)));
            x = 1+sum(Stim.NoPolarAngles(1:e-1));
            Indx = x:(x+Stim.NoPolarAngles(e)*Stim.NoDirections-1);
            Stim.Variables(Indx,1) = Stim.Eccentricities(e);
            for p = 1:Stim.NoPolarAngles(e)
                x = 1+sum(Stim.NoPolarAngles(1:e-1))+((p-1)*Stim.NoDirections);
                Indx = x:(x+Stim.NoDirections-1);
                Stim.Variables(Indx,2) = Stim.PolarAngles{e}(p);
            end
        end
        Stim.Variables(:,3) = repmat(Stim.Directions,[1,Stim.NoConditions/Stim.NoDirections]);

        %================ Dot motion parameters
        Dots.Num        = 300;
        Dots.Signal     = 1;
        Dots.SpeedMean  = 5;             % Degrees per second
        Dots.SpeedStd   = 0;
        Dots.LifeMean   = 10;
        Dots.LifeStd    = 0;
        Dots.DirStd     = 0;
        Dots.SizeMean   = 0.07*Display.Pixels_per_deg(1);
        Dots.SizeStd    = 0.05*Display.Pixels_per_deg(1);

      	Grating.CyclesPerDeg = 1;       % Number of grating cycles per degree
        Grating.Phase = 0;              % Phase offset
        
        %======================== PRESENT FIXATION ========================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
%             if ~isfield(Fix,'Texture')
                Fix = fixOn(Fix);
%                 FixDrawn = 1;
%             else
%                 FixDrawn = 0;
%             end
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
%                 if FixDrawn == 0
%                     Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
%                 end
            end
        	[VBL FixOn] = Screen('Flip', Display.win, [], 1);
        	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);  
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            When = FixOn +(1/Display.RefreshRate)+0.001;
         	[VBL FixOff] = Screen('Flip', Display.win, When,[],Display.DontSync,Display.MultiFlip);
            if Debug.On == 1
                fprintf('Photodiode target was on for %.3f ms\n', (FixOff-FixOn)*1000);
            end
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});
            end 
            [VBL FixOff] = Screen('Flip', Display.win);                             % Clear stimulus  
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end 
          	Screen('Flip', Display.win);                                           % Clear photodiode target
        end
        
        %========================= GENERATE STIMULI =======================
        function R_makeStim(varargin)
            params = varargin{1,1};
            Stim.Window = params([4,5])*Display.Pixels_per_deg(1);                      % Update stimulus size
            Stim.Type = params(6);
          	Stim.Rect = [0 0 Stim.Window];                                          
            switch Stim.Type
                case 1                                                                  % For moving dot stimuli...
                    Mask.Edge = 0;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
                    if Stim.Params(1:5) ~= params(1:5)                                  % If parameters have changed...
                        Dots.Colour = params([1,2,3])*255;
                        Dots.Window = Stim.Window;
                        Dots.Velocity = Dots.SpeedMean*Display.Pixels_per_deg(1);       % Velocity in pixels per second
                        Dots.Background = Display.Background;   
                        Dots.Size = [Dots.SizeMean, Dots.SizeStd, 1];
                        StimTexture = GenerateRDK(Dots, Display);                       % Generate motion dot texture franes
                        Stim.NrFrames = numel(StimTexture);     
                    end

                case 2                                                                  % For static Gabor patch stimuli...
                    Mask.Edge = 1;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
                    if Stim.Params(1:5) ~= params(1:5)
                        Grating.Colour = params([1,2,3])*255;
                        Grating.Dim = Stim.Window;
                        Grating.CyclesPerDeg = 1;
                        StimTexture(1) = GenerateSineGrating(Grating, Display);
                        Stim.NrFrames = 1;
                    end
                
                case 3                                                                  % For a static colored blob...
                    Mask.Edge = 0;
                    if Stim.Params(1:3) ~= params(1:3) 
                        Stim.Colour = params([1,2,3])*255;
                        TextureMat = ones(Stim.Window,3);
                        for RGB = 1:3
                            TextureMat(:,:,RGB) = Stim.Colour(RGB);
                        end
                        StimTexture(1) = Screen('MakeTexture',Display.win, TextureMat);
                        Stim.NrFrames = 1;
                    end
                    
                case 4
                    Mask.Edge = 0;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
                    if Stim.Params(1:5) ~= params(1:5)
                        Stim.NoFrames = 10;
                        Stim.NoDots = 400;
                        Stim.Background = Display.Background;
                        Stim.DotColor = [0 255];
                        Stim.DotSize = [0.05 0.15]*Display.Pixels_per_deg(1);
                        StimTexture = GenerateDynamicNoise(Stim,Display);
                        Stim.NrFrames = numel(StimTexture);     
                    end
            end
            %================= CREATE ALPHA CHANNEL MASK ==================
            if Stim.Params([4,5]) ~= params([4,5])                                  % If requested stimulus size has changed...
                Mask.Dim = Stim.Window+4;                                           % Mask will be 2 pixels larger than the stimulus window
                Mask.ApRadius = (min(Mask.Dim)/2)-8;
                Mask.Colour = Display.Background(1);
                Mask.S = Mask.ApRadius/3;                                           % Standard deviation of Gaussian envelope
                MaskTex = GenerateAlphaMask(Mask, Display);
                Mask.Rect = [0 0 Mask.Dim];                                      	% Mask size
            end
%             if Stim.Params([6,7]) ~= params([6,7])                                  % If requested stimulus position has changed...
                Stim.WindowPos = [0 0];
                Stim.DestRect = CenterRect(Stim.Rect, Display.Rect) + [Stim.WindowPos, Stim.WindowPos];
                Mask.DestRect = CenterRect(Mask.Rect, Display.Rect) + [Stim.WindowPos, Stim.WindowPos]; 
%             end
            Stim.Params = params;                                                   % Save newest stimulus parameters to structure
            
        end
        
        
        %======================== PRESENT STIMULI =========================
        function R_stimOn(varargin)
            params      = varargin{1,1};
            Eye         = params(1);                        
            Condition   = params(2);
            if Eye==2 
                while Condition > size(Stim.Variables,1)
                    Condition = Condition - Stim.NoConditions;
                end
            end
            Eccentricity    = Stim.Variables(Condition,1)*Display.Pixels_per_deg(1);
            PolarAngle      = Stim.Variables(Condition,2);
            Orientation     = Stim.Variables(Condition,3);
            Diameter        = Stim.Diameters(find(Stim.Variables(Condition,1)==Stim.Eccentricities));
          
            X = Eccentricity*sind(PolarAngle);
            Y = Eccentricity*cosd(PolarAngle);
            Stim.DestRect = CenterRect([0 0 Diameter Diameter],Display.Rect)+[X Y X Y];
            Mask.DestRect = Stim.DestRect;
            
            Frame = 1;
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            while 1
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    fprintf(strcat(commandIn,'\n'));
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn,'R_stimOff')            % Exits loop if trial is quit or ended
                    break
                end

                %===================== DRAW STIMULI =======================
                if Eye < 2
                  	currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);                       % Select which eye to draw to
                    Screen('DrawTexture', Display.win, StimTexture(Frame), Stim.Rect, Stim.DestRect, Orientation, [], Stim.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Orientation); 	% Apply Gaussian aperture mask
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                    Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye+1});
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, abs(floor(Eye-0.5)));   	% Select which eye to draw to
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                elseif Eye == 2
                    for E = 1:2
                        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E-1);                       % Select which eye to draw to
                        Screen('DrawTexture', Display.win, StimTexture(Frame), Stim.Rect, Stim.DestRect, Orientation, [], Stim.Contrast); 
                        Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect, Orientation); 	% Apply Gaussian aperture mask
                        Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                        Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{E});
                    end
                end
                [VBL FrameOn] = Screen('Flip', Display.win);

                %================== ANIMATION LOOP
                Frame = Frame+1;
                if Frame >= Stim.NrFrames
                    Frame = 1;
                end
                [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                if keyIsDown && keyCode(Key.Exit)
                    break;
                end
            end
            R_stimOff;
        end
        
        
        %=========================== CLEAR TEXTURE ========================
        function R_stimOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background);
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus  
        end
        
        function R_getBlock(n)

        end
    end


%% ============================ apm_attention =============================
% Variant of the Posner cued attention paradigm.
%
%
%==========================================================================

    function AttendHandle = apm_attention()

        AttendHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_setDisksColor;...
                @R_lineLocation;...
                @R_lineWidth;...
                @R_disksOnFC;...
                @R_disksOnFT;...
                @R_disksOnAll;...
                @R_lineOn;...
                @R_diskChange;...
                @R_fixOff;...
                @R_disksOff;...
                @R_lineOff;...
                @R_getBlock;...
            };
        
        persistent Disk;
        persistent Fix;
        persistent Line;

     	%======================= SET STIMULUS PARAMETERS ======================
        Fix.Type = 0;                                   % 0 = circle; 1 = Cross; 2 = Crosshair
        Fix.LineWidth = 3;                              % Line width in pixels
        Line.Eye = 1;
        Disk.Number = 8;                                % Total number of target disks
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
        %======================== PRESENT FIXATION ========================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(4)*Display.Pixels_per_deg(1);
            Fix.Eye = params(3);
%             Fix.Eye = abs(floor(params(3)-0.5));
            Fix = fixOn(Fix);
            [VBL FixOn] = Screen('Flip', Display.win, [], 1);
        end
        
      	%========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win, [], 1);                                         % Clear stimulus  
        end
        
        %====================== GET DISK PARAMETERS =======================
        function R_setDisksColor (varargin)
            params=varargin{1,1};
            Disk.Colour = [params(1) params(2) params(3)]*255;
            Disk.Ecc = params(4)*Display.Pixels_per_deg(1);
            Disk.Size = params(5)*Display.Pixels_per_deg(1);
            Disk.Rect = CenterRect([0 0 Disk.Size Disk.Size], Display.Rect);
            for d = 1:Disk.Number
            	Angle = (d-1)*360/Disk.Number;
                Disk.Positions(d,:) = Disk.Ecc*[cosd(Angle) sind(Angle)];
                Disk.DestRects(d,:) = Disk.Rect+[Disk.Positions(d,:),Disk.Positions(d,:)];
            end
        end
        
        %======================== GET CUE PARAMETERS ======================
        function R_lineLocation(varargin)
            params = varargin{1,1};
            Line.Locations = params;
            Line.XCoord = params([1,4])*Display.Pixels_per_deg(1)+Display.Centre(1);
            Line.YCoord = params([5,8])*Display.Pixels_per_deg(1)+Display.Centre(2);
        end

        function R_lineWidth(varargin)
            params = varargin{1,1};
            Line.Width = params(1);
        end
        
        function R_lineOn(varargin)
            params = varargin{1,1};
            Line.Eye = abs(floor(params(1)-0.5));
            Line.Colour = params([2,3,4]);
            currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Line.Eye);
            Screen('DrawLine', Display.win, Line.Colour, Line.XCoord(1), Line.YCoord(1), Line.XCoord(2), Line.YCoord(2), Line.Width); 
            Screen('Flip', Display.win, [], 1);
        end
        
        function R_lineOff(varargin)
            if isfield(Line, 'XCoord')
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Line.Eye);
                Screen('DrawLine', Display.win, Display.Background, Line.XCoord(1), Line.YCoord(1), Line.XCoord(2), Line.YCoord(2), Line.Width+2); 
                Screen('Flip', Display.win, [], 1);
            end
        end
        
        %==================== DRAW ALL DISKS ==============================
        function R_disksOnFC(Eye)
            Disk.Eye = abs(floor(Eye-0.5));
            currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Disk.Eye);
        	Screen('FillOval', Display.win, Disk.Colour, Disk.DestRects(1:2:end,:)', Disk.Size+1);
            Screen('Flip', Display.win, [], 1);
        end
        
      	function R_disksOnFT(Eye)
            Disk.Eye = abs(floor(Eye-0.5));
            currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Disk.Eye);
        	Screen('FillOval', Display.win, Disk.Colour, Disk.DestRects(2:2:end,:)', Disk.Size+1);
            Screen('Flip', Display.win, [], 1);
        end
        
     	function R_disksOnAll(Eye)
            Disk.Eye = abs(floor(Eye-0.5));
            currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Disk.Eye);
        	Screen('FillOval', Display.win, Disk.Colour, Disk.DestRects', Disk.Size+1);
            Screen('Flip', Display.win, [], 1);
        end
        
        %==================== DRAW TARGET DISK ============================
        function R_diskChange(varargin)
            params = varargin{1,1};
            Disk.TargetDisk = params(1);
            Disk.TargetColour = params([2,3,4]);
            currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Disk.Eye);
        	Screen('FillOval', Display.win, Disk.TargetColour, Disk.DestRects(Disk.TargetDisk,:)', Disk.Size+1);
            Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Disk.Eye+1});
            Screen('Flip', Display.win, [], 1);
        end

        %========================= CLEAR DISKs ============================
        function R_disksOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Display.Rect);
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end 
            [VBL FixOff] = Screen('Flip', Display.win);                                         % Clear stimulus  
        end
    end

%% ============================== apm_SFM =================================
% Generates animations of ambiguous and disparity-defined unambiguous 
% structure-from-motion (SFM) stimuli, based on input structure 'Stim'.  
% Returns animation frames as PsychToolbox textures, with each animation
% contained within a separate cell of the arrays for left and right eye
% images.
%
%       MaxDisparity            Disparity in arcminutes of dots at the
%                               closest point of the object to the observer.
%       IPD                     Observer's interpupillary distance (metres)
%
% STIMULUS PARAMETERS:
%       Stim.Type:              0 = sphere, 1 = cylinder
%       Stim.Hemi:              0 = full object, 1 = half object
%       Stim.CVCC               1 = convex, 2 = concave
%       Stim.AngularVelocity:   angular velocity of rotation (°/s)
%       Stim.MaxDisparity       maximum disparity for catch period stimuli (arcmin)
%      	Stim.Veridical          0 = use max disparity, 1 = use geometrically correct disparity    
%     	Stim.Radius             radius of sphere or cylinder (°)
%     	Stim.Height             stimulus height (°)
%    	Stim.NoDots             total number of dots
%      	Stim.DotSize            dot size (°)
%      	Stim.DotColour          0 = black & white, 1 = single colour
%       Stim.DotColourRGB       RGB value to use for single colour dots
%     	Stim.NoFrames           = round((360/Stim.AngularVelocity)*Display.RefreshRate);   
%      	Stim.RadPerFrame        rotation between consecutive frames in radians
%      	Stim.Background         background colour  
%       Stim.DotType            0 = circular, 1 = square, 2 = antialiased
%       Stim.DotLifetime        Default: inf = do not replace. (Frames)
%       Stim.Tilt               axis tilt (CW/ACW) in degrees from vertical
%       Stim.Slant              axis slant (forward/ backward) in degrees from vertical
%       Stim.Shading            Matrix containing shading model
%       Stim.ShadingNoise       Amount of noise to add to shading signal (%)        
%       Stim.Filter             1 = apply Gaussian filter to dots. NB. This
%                               can drastically increase time taken to
%                               generate stimuli.  Consider alternative
%
%==========================================================================

    function SFMHandle = apm_SFM()

            SFMHandle = ...
                {...
                    @R_punishStim;... 
                    @R_fixOn;...
                    @R_setFixColor;...
                    @R_fixOff;...
                    @R_setDots;...
                    @R_StimOn;...
                    @R_StimOff;...
                    @R_getBlock;...
                };

            persistent Stim;
            persistent Fix;
            persistent IPD;
            persistent Intermittent;
            persistent BackgroundTexture;
            persistent BackgroundOn;
            persistent SFMTextureL;
            persistent SFMTextureR;
            
            BackgroundOn = 1;
            BackgroundTexture = Screen('MakeTexture', Display.win, ones(Display.Rect([4,3]))*Display.Background(1));
            IPD = 0.035;                                        % Average rhesus macaque interpupillary distance (metres)
            Fix.Type = 0;
            Fix.Eye = 2;
            
            %================= SET STIMULUS PARAMETERS ====================
            Stim.Type = 0;                                      %0 = sphere, 1 = cylinder
            Stim.Hemi = 0;                                      %0 = full object, 1 = half object
            Stim.CVCC  = 1;                                     %1 = convex, 2 = concave

            Stim.MaxDisparity = 10;                              % maximum disparity for unambiguous stimuli (arcmin)
            Stim.Veridical = 1;                                 %0 = use max disparity, 1 = use geometrically correct disparity    
            Stim.Direction = 2;
            Stim.Contrast = 1*255;
            
            Stim.Background = Display.Background;               %background colour  
            Stim.DotType = 2;                                   %0 = circular, 1 = square, 2 = antialiased
            Stim.DotLifetime = inf;                             %Default: inf = do not replace. (Frames)

        	Intermittent.On = 0;                
            Intermittent.OnDur = 1;             % Set duration of stimulus on period (s)
            Intermittent.OffDur = 1;            % Set duration of blank period (s)
            Intermittent.StimOn = 1;            % Start trial with stimulus 'on'
 
            %====================== GET FIXATION COLOUR ===================
            function R_setFixColor (varargin)
                params=varargin{1,1};
                Fix.Colour = [params(1) params(2) params(3)]*255;
            end
          	%======================== PRESENT FIXATION ========================
            function R_fixOn (varargin)
                params=varargin{1,1};
                Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
                Fix.Size = params(3)*Display.Pixels_per_deg(1);
                Fix = fixOn(Fix);
            	[VBL FixOn] = Screen('Flip', Display.win);
            end

            %================== GENERATE ANIMATION FRAMES =================
            function R_setDots(varargin)
                params = varargin{1,1};
                
                Stim.NoDots = params(1);
                Stim.DotColourRGB = params([2,3,4]);
                Stim.DotSize = params(5)*Display.Pixels_per_deg(1);
                Stim.AngularVelocity = params(6);
                Stim.Tilt = params(7);
                Stim.Hemi = params(8);
                Stim.Radius = params(9)*Display.Pixels_per_deg(1)/2;    %radius of sphere or cylinder (°)
                Stim.Height = params(10)*Display.Pixels_per_deg(1);  	%stimulus height (°)
                BackgroundOn = params(11);
                
                Stim.DotColour = 1;     % 1 = use RGB value 
                Stim.NoFrames = round((360/Stim.AngularVelocity)*Display.RefreshRate); 
                Stim.RadPerFrame = 2*pi/Stim.NoFrames;
                Stim.Rect = [0 0 Stim.Height, Stim.Height]+[0 0 Stim.DotSize, Stim.DotSize]*2;
                Stim.DestRect = CenterRect(Stim.Rect,Display.Rect);
                
                if Debug.On == 1
                    DisplayText('Loading...',Display);
                end
                MaxDisparity = 20;
                Display.Pixels_per_m = Display.Pixels_per_m(1);
                [SFMTextureL, SFMTextureR] = GenerateSFM(Display, Stim, MaxDisparity, IPD);
                
                if BackgroundOn == 1
                    BkgParam.Window = [0 0 500 500];
                    BkgParam.Background = Display.Background(1);
                    BackgroundTexture = BackgroundSquares(Display, BkgParam);        	% Generate squares background
                elseif BackgroundOn == 0
                    BackgroundTexture = Screen('MakeTexture', Display.win, ones(Display.Rect([4,3]))*Display.Background(1));
                end
            end
        
            %================= PRESENT STIMULUS ANIMATION =================
            function R_StimOn(Eye)
                Frame = 1;
                pnet(con,'setreadtimeout',.001);
                pnet(con,'printf', reply);

                for E = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, E-1);
                    Screen('DrawTexture', Display.win, BackgroundTexture);
                end
                Screen('Flip', Display.win, [], 1);
                
                while 1
                    commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                    if ~isempty(commandIn)
                        fprintf(strcat(commandIn,'\n'));
                        pnet(con,'printf', reply);              % return control to the socket
                    end
                    if strcmp(commandIn,'R_StimOff')            % Exits loop if trial is quit or ended
                        break
                    end
                
                    %===================== DRAW STIMULI =======================
                    Eye = 0;
                    if Intermittent.StimOn == 1
                        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);
                        Screen('DrawTexture', Display.win, BackgroundTexture);
                        if Stim.MaxDisparity == 0
                            Screen('DrawTexture', Display.win, SFMTextureL{1}(Frame), Stim.Rect, Stim.DestRect, [], [], Stim.Contrast);
                        else
                            Screen('DrawTexture', Display.win, SFMTextureL{Stim.Direction}(Frame), Stim.Rect, Stim.DestRect, [], [], Stim.Contrast);
                        end
                        Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye+1});
                    elseif Intermittent.StimOn == 0
                        Screen('FillRect', Display.win, Display.Background, Stim.DestRect);
                        Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye+1});
                    end
                    if Stim.MaxDisparity ~= 0
                        Eye = 1;
                        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);
                        Screen('DrawTexture', Display.win, BackgroundTexture);
                        Screen('DrawTexture', Display.win, SFMTextureR{Stim.Direction}(Frame), Stim.Rect, Stim.DestRect, [], [], Stim.Contrast);
                        Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye+1});
                    end
                    [VBL FrameOn] = Screen('Flip', Display.win);

                    %================== ANIMATION LOOP
                    Frame = Frame+1;
                    if Frame == Stim.NoFrames
                        Frame = 1;
                    end
                    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                    if keyIsDown && keyCode(Key.Exit)
                        break;
                    end  
                    
                    %================= INTERMITTENT PRESENTATION
                    if Intermittent.On == 1
                        if Intermittent.StimOn == 1 && GetSecs >= Intermittent.SwitchOn+Intermittent.OnDur
                            Intermittent.StimOn = 0;
                            Intermittent.SwitchOff = GetSecs;
                        elseif Intermittent.StimOn == 0 && GetSecs >= Intermittent.SwitchOff+Intermittent.OffDur
                            Intermittent.StimOn = 1;
                            Intermittent.SwitchOn = GetSecs;
                        end
                    end
                end
                R_StimOff;
            end
        
            %================= PRESENT STIMULUS ANIMATION =================
            function R_StimOff()
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                    Screen('DrawTexture', Display.win, BackgroundTexture);
                    Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
                end
                [VBL StimOff] = Screen('Flip', Display.win);
                for n = 1:3
                    Screen('Close', SFMTextureL{n});
                    Screen('Close', SFMTextureR{n});
                end
            end
    end


%% ============================ apm_SFMPA ==============================


%% =========================== StereoTest =================================
% A simple test for stereoscopic vision in monkeys/ infants, based on a
% 9-point calibration routine. On each trial, a random dot stereogram is 
% presented, in which a target (defined exclusively by binocular disparity)
% appears in one of 9 possible spatial locations. The observer's task is to
% saccade to the target.
%
%==========================================================================

    function StereoTestHandle = apm_StereoTest()

        StereoTestHandle = ...
            {...
                @R_punishStim;...
                @R_fixOn;...
                @R_setParams;...
                @R_setFixColor;...
                @R_remoteTargetOn;...
                @R_getStmParamByName;...
                @R_getStmParamName;...
                @R_fixOff;...
                @R_getStmNParams;...
                @R_getStmParam;...
                @R_getBlock;...
                @R_nextTrial;...
                @R_clearscreen;...
                @R_reset;...
                @R_clearstim;...
                @R_getFixPosX;...
                @R_getFixPosY;...
                @R_getFixEye;...
                @newStmParam;...
                @stereocalibReset;...
            };

        persistent stereocalib_xpos;
        persistent stereocalib_ypos;
        persistent stereocalib_trial_num;
        persistent stereocalib_eye;
        persistent stereocalib_block;
        persistent stereocalib_fixcolor;
        persistent Dots;
        persistent Target;
        persistent BorderSquares;
        persistent Cue;
        
        
        %===================== GENERATE RANDOM DOT POSITIONS ==============
%         Target.Diparities   = [10/60, 20/60 0.5, 0.75, 1];                  % 
%         Target.Disp       	= 0.25*Display.Pixels_per_deg(1);               % +ve = near, -ve = far
                    Cue.On      = 1;                                         	% Cue to target?
        Cue.Width   = 4;                                            % Cue width
        Cue.Colour  = [255 0 0 20];                              	% Cue color RGBA

        Fix.Colour  = [0 0 0];
        Fix.Size    = 1*Display.Pixels_per_deg(1);
        Fix.Rect    = CenterRect([0 0 Fix.Size Fix.Size], Display.Rect);
        
        
        function R_setParams(varargin)
            params = varargin{1,1};
            Dots.Size           = 3;
            Dots.Density        = params(2);
            Target.Disp       	= params(1)/60*Display.Pixels_per_deg(1);    
            Target.Ecc          = params(3);
            Dots.Rect           = round([0 0 18 18]*Display.Pixels_per_deg(1));
            Dots.DestRect       = CenterRect(Dots.Rect, Display.Rect);
            Dots.TotalArea      = Dots.Rect(3)*Dots.Rect(4);
            Dots.AreaPerDot     = 2*pi*Dots.Size/2;
            Dots.PerFrame       = round(Dots.TotalArea/Dots.AreaPerDot*Dots.Density);
            Dots.Centre         = Dots.Rect([3,4])/2;
            Dots.Colour         = round(rand(1,Dots.PerFrame))*255;             % Dots are black and white
            Dots.Colour         = repmat(Dots.Colour,[3,1]);                                                      
            for P = 1:9
               Dots.x{P} = (rand([1,Dots.PerFrame])*(Dots.DestRect(3)-Dots.DestRect(1)))-Dots.Centre(1); 
               Dots.y{P} = (rand([1,Dots.PerFrame])*(Dots.DestRect(4)-Dots.DestRect(2)))-Dots.Centre(2);
            end

            Stim.Window         = Dots.Rect;
            Stim.Background     = Display.Background;
            BorderSquares       = BackgroundSquares(Display, Stim);
        end
        

        %======================== PRESENT FIXATION MARKER =================
        function R_fixOn (varargin)
            
            %================== Present central fixation
         	for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);               % Select which eye to draw to
                Screen('DrawTexture', Display.win, BorderSquares);                                  % Draw background texture
                Screen('DrawDots', Display.win, [Dots.x{1}; Dots.y{1}], Dots.Size, Dots.Colour, Display.Centre, 2);      % Draw random dot field
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});      	% Draw Photodiode target
                Screen('FillOval', Display.win, [0 0 0], CenterRect([0 0 20 20], Display.Rect));   	% Draw explicit cue to target location
            end
            [CentralFixOnset, VBL] = Screen('Flip', Display.win);
            

            if (size(varargin,1)~=0)
                params=varargin{1,1};
                Target.Pos      = [params(1); params(2)]*Display.Pixels_per_deg(1);     % Coordinates of target center (pixels)
                Target.Size     = params(4)*Display.Pixels_per_deg(1);                  % Target diameter (pixels) 
                Target.Disp     = params(5)/60*Display.Pixels_per_deg(1);           	% Target disparity (convert from arcmin to pixels) (+ve = near, -ve = far)
                stereocalib_eye = params(3);
                Cue.Colour(4)   = params(6)*255;
            end
            P = randi(numel(Dots.x));
            for d = 1:Dots.PerFrame         % Add horizontal disparity to dots located inside the target
                if (Dots.x{P}(d)-Target.Pos(1))^2 + (Dots.y{P}(d)-Target.Pos(2))^2 <= (Target.Size/2)^2
                    X{P,1}(d) = Dots.x{P}(d)+Target.Disp/2;
                    X{P,2}(d) = Dots.x{P}(d)-Target.Disp/2;
                else
                    X{P,1}(d) = Dots.x{P}(d);
                    X{P,2}(d) = Dots.x{P}(d);
                end
            end
            
            Cue.Rect = [Display.Centre,Display.Centre]+[-Target.Size/2+Target.Pos(1),-Target.Size/2+Target.Pos(2),Target.Size/2+Target.Pos(1),Target.Size/2+Target.Pos(2)];
            Cue.DestRect{1} = Cue.Rect+[Target.Disp/2,0,Target.Disp/2,0];
            Cue.DestRect{2} = Cue.Rect-[Target.Disp/2,0,Target.Disp/2,0];
            
            for Eye = 1:2
                XY = [X{P,Eye}; Dots.y{P}];
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);               % Select which eye to draw toca
                
                Screen('DrawTexture', Display.win, BorderSquares);                                  % Draw background texture
                Screen('DrawDots', Display.win, XY, Dots.Size, Dots.Colour, Display.Centre, 2);      % Draw random dot field
%                 Screen('FillOval', Display.win, Fix.Colour, Fix.Rect);                               % Draw central fixation marker
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye});            % Draw Photodiode target
                if Cue.On == 1
                    Screen('FillOval', Display.win, Cue.Colour, Cue.DestRect{Eye});                  	% Draw explicit cue to target location
                end
            end
            CentralFixDur = 0.15;
            while GetSecs < CentralFixOnset + CentralFixDur
                
            end
            Screen('Flip', Display.win);
        end

        function R_setFixColor (varargin)
            params=varargin{1,1};
            stereocalib_fixcolor = [params(1) params(2) params(3)];
        end

        function R_remoteTargetOn ()

        end

        function R_getStmParamByName (varargin)
            params=varargin{1,1};
            if strcmp(params,'xpos')
                R_getFixPosX();
            elseif strcmp(params,'ypos')
                R_getFixPosY();
            elseif strcmp(params,'eye')
                R_getFixEye();
            else
                reply = 'ERROR\n';
            end
        end

        function R_getStmParamName ()
            reply = 'xpos\n';
        end

%         function R_fixCenter() % Central fix on
%             for Eye = 1:2
%                 currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);               % Select which eye to draw to
%                 Screen('DrawTexture', Display.win, BorderSquares);                                  % Draw background texture
%                 Screen('DrawDots', Display.win, [Dots.x{1}; Dots.y{1}], Dots.Size, Dots.Colour, Display.Centre, 2);      % Draw random dot field
%                 Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});      	% Draw Photodiode target
%                 Screen('FillOval', Display.win, [0 0 0], CenterRect([0 0 20 20], Display.Rect));   	% Draw explicit cue to target location
%             end
%             Screen('Flip', Display.win);
%         end

        function R_fixOff() % Fixation off
             for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);               % Select which eye to draw to
                    Screen('DrawTexture', Display.win, BorderSquares);                                  % Draw background texture
                    Screen('DrawDots', Display.win, [Dots.x{1}; Dots.y{1}], Dots.Size, Dots.Colour, Display.Centre, 2);      % Draw random dot field
                    Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});      	% Draw Photodiode target
             end
             Screen('Flip', Display.win);
        end
        
        function R_getStmNParams ()
                reply = '3\n';
        end

        function R_getStmParam ()
                reply = '0.0\n';
        end

        function R_getBlock(n)
            stereocalib_trial_num = 0;
            A = perms([-1,0,1]);
            A = [A(:,1:2); -1,-1;0,0;1,1];
            A = [A, zeros(9,1); A, ones(9,1)];
            A = [A,randperm(18)'];
            A = sortrows(A,4);
            stereocalib_block = A(:,1:3);
            reply = sprintf('%d\n',size(stereocalib_block,1));
        end

        function R_nextTrial()
             if(stereocalib_trial_num < size(stereocalib_block,1)-1)
                stereocalib_trial_num = stereocalib_trial_num + 1;
             else
                 stereocalib_trial_num = -1;
                 reply = sprintf('-1\n');
                 return;
             end
             stereocalib_xpos = stereocalib_block(stereocalib_trial_num,1);
             stereocalib_ypos = stereocalib_block(stereocalib_trial_num,2);
             reply = sprintf('%d %d %d %d\n',stereocalib_trial_num,stereocalib_eye,stereocalib_xpos+1,stereocalib_ypos+1);
        end

        function R_clearscreen()         

        end

        function R_reset ()

        end

        function R_clearstim ()

        end

        function R_getFixPosX(Ecc) 
            reply = sprintf('%d\n',stereocalib_xpos*Ecc); 
        end

        function R_getFixPosY(Ecc)
            reply = sprintf('%d\n',stereocalib_ypos*Ecc);           
        end

        function R_getFixEye()
            reply = sprintf('%d\n',stereocalib_eye);
        end
    end


%% ============================ apm_DelaySac =============================
% Present moving dot or grating patches across the visual field for
% delayed saccade mapping.
%
%==========================================================================

    function DelaySacHandle = apm_DelaySac()

        DelaySacHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_stimOn;...
                @R_stimOff;...
                @R_makeStim;...
                @R_getBlock;...
            };
        
        persistent Dots;
        persistent DotsTexture;
        persistent Fix;
        persistent Mask;
        persistent MaskTex;
        persistent Stim;
        
        %======================= SET STIMULUS PARAMETERS ======================
        Fix.Type = 0;
        Stim.Type = 2;          % 1 = moving dots; 2 = static disk
        
        Dots.Num = 300;
        Dots.Signal = 1;
        Dots.Contrast = 1;
        Dots.SpeedMean = 5;     % Degrees per second
        Dots.SpeedStd = 0;
        Dots.LifeMean = 10;
        Dots.LifeStd = 0;
        Dots.DirStd = 0;
        Dots.SizeMean = 0.07*Display.Pixels_per_deg(1);
        Dots.SizeStd = 0.05*Display.Pixels_per_deg(1);
        Dots.Params = zeros(1,8);
        
        %======================== PRESENT FIXATION ========================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            Fix.Eye = 2;
            Fix.On = 1;
            Fix = fixOn(Fix);
          	[VBL FixOn] = Screen('Flip', Display.win);
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win);                                         % Clear stimulus  
        end
        
        %========================== SET DOT PARAMETERS ====================
        function R_makeStim(varargin)
            params = varargin{1,1};
            
            if Dots.Params(1:6) ~= params(1:6)                                           % If parameters have changed...
                Dots.Colour = params([1,2,3])*255;
                Dots.DirMean = params(4)-90;
                Dots.Window = [params(5) params(5)]*Display.Pixels_per_deg(1);
                Dots.Rect = [0 0 Dots.Window];
                
                %============== GENERATE TRANSPARENT MOTION TEXTURES ==========
                Dots.Velocity = Dots.SpeedMean*Display.Pixels_per_deg(1);           % Velocity in pixels per second
                Dots.Background = Display.Background;   
                Dots.DrawAngle = Dots.DirMean;
                Dots.Size = [Dots.SizeMean, Dots.SizeStd, 1];
                DotsTexture = GenerateRDK(Dots, Display);
                Dots.NrFrames = numel(DotsTexture);
            end

            %================= CREATE ALPHA CHANNEL MASK ==================
            if Dots.Params(5) ~= params(5)                                          % If requested stimulus size has changed...
                Mask.Dim = Dots.Window+4;                                           % Mask will be 2 pixels larger than the stimulus window
                Mask.ApRadius = (min(Mask.Dim)/2)-8;
                Mask.Colour = Display.Background(1);
                Mask.Edge = 0;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
                MaskTex = GenerateAlphaMask(Mask, Display);
                Mask.Rect = [0 0 Mask.Dim];                                      	% Mask size
            end
            if Dots.Params([6,7]) ~= params([6,7])                                  % If requested stimulus position has changed...
                Dots.WindowPos = params([6,7])*Display.Pixels_per_deg(1);
                Dots.DestRect = CenterRect(Dots.Rect, Display.Rect) + [Dots.WindowPos, Dots.WindowPos];
                Mask.DestRect = CenterRect(Mask.Rect, Display.Rect) + [Dots.WindowPos, Dots.WindowPos]; 
            end
            Dots.Params = params;
        end
        
        %======================== PRESENT STIMULI =========================
        function R_stimOn(Eye)
            Frame = 1;
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            
            while 1
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    fprintf(strcat(commandIn,'\n'));
                    pnet(con,'printf', reply);              % return control to the socket

                    if strcmp(commandIn,'R_stimOff')            % Exits loop if trial is quit or ended
                        break
                    elseif strcmp(commandIn,'R_fixOff')
                        R_fixOff
                        Fix.On = 0;
                    end
                end
                
                %===================== DRAW STIMULI =======================
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);
                if Stim.Type == 1
                    Screen('DrawTexture', Display.win, DotsTexture(Frame), Dots.Rect, Dots.DestRect, Dots.DirMean, [], Dots.Contrast); 
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);                  % Apply Gaussian aperture mask
                elseif Stim.Type == 2
                    Screen('FillOval', Display.win, Dots.Colour, Dots.DestRect);
                end
                if Fix.On==1
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                end
                Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{Eye+1});
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, abs(floor(Eye-0.5)));     % Select which eye to draw to
                if Fix.On==1
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);
                end
                [VBL FrameOn] = Screen('Flip', Display.win);

                %================== ANIMATION LOOP
                Frame = Frame+1;
                if Frame == Dots.NrFrames
                    Frame = 1;
                end
                [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                if keyIsDown && keyCode(Key.Exit)
                    break;
                end
            end
            R_stimOff
        end
                
        %=========================== CLEAR TEXTURE ========================
        function R_stimOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background);
                Screen('DrawTexture', Display.win, Fix.Texture); 
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus  
        end
        
        function R_getBlock(n)

        end
    end


%% ============================== apm_CFS =================================
% Present continuous flash suppresion mask to one eye and a static image to
% the other eye. 
%
%==========================================================================

    function CFSHandle = apm_CFS()

        CFSHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_stimOn;...
                @R_stimOff;...
                @R_makeCFSStim;...
                @R_getBlock;...
            };
        
        persistent CFS;
        persistent CFSTextures;
        persistent Mask;
        persistent MaskTex;
        persistent Fix;
        persistent Target;
        persistent ImageTextures;
        persistent BorderTexture;
        

        %======================== PRESENT FIXATION ========================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            Fix.Eye = 2;
            Fix.On = 1;
            Fix = fixOn(Fix);
          	[VBL FixOn] = Screen('Flip', Display.win);
            CFS.On = 0;
            Target.On = 0;
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win);                                         % Clear stimulus  
        end
        
        %======================= GENERATE CFS TEXTURES ====================
        function R_loadImages(varargin)
            params = varargin{1,1};
%             CFS.FrameSize = params([1,2])*Display.Pixels_per_deg(1);
%             CFS.TexelsPerFrame = params(5);
%             CFS.NoFrames = params(4)*Display.RefreshRate;
%             CFS.Freq = params(3);
%             Image.Radius = params(6)*Display.Pixels_per_deg(1);           
            
            %==================== GENERATE CFS MASK 
            if exist('CFSTextures','var')~= 2
                CFS.FrameRate = 10;                                             % Frame rate (Hz)
                CFS.Duration = 5;                                      
                CFS.TextureSize = [12, 12];
                CFS.Background = Display.Background;

                if exist('CFSTextures','var') 
                    Screen('Close', CFSTextures);                           	% close any previously loaded textures
                end
    %             CFSTextures = GenerateCFS(Display.win, CFS.FrameSize, CFS.TexelsPerFrame, CFS.NoFrames, CFS.Background);
                CFSTextures = GenerateDCFS(DCFS, Display,0);
            end
            
            %==================== LOAD IMAGE(S)        
            if exist('ImageTexture','var')                                  
                Screen('Close', ImageTexture);                           	% close any previously loaded image textures
            end
            AllImageTypes = regexp(genpath(ImageDir),['[^;]*'],'match');    % Get list of image subcategory folders
            CatgeoryDir = AllImageTypes{ImageType+1};                    	% Get path of requested image category
            [x ImageCategory] = fileparts(CatgeoryDir);                   	% Get text string for category
            Image.Dir = ImageCategory;
            Image.No = 1:params(7);
            Image.Mask = 1;
            [ImageTextures, Image] = LoadImages(Image);
            
            %==================== CREATE ALPHA MASK
            if exist('MaskTex','var')~= 2
                Mask.Dim = Dots.Window+4;                                           % Mask will be 2 pixels larger than the stimulus window
                Mask.ApRadius = (min(Mask.Dim)/2)-8;
                Mask.Colour = Display.Background(1);
                Mask.Edge = 0;                                                      % Set edge type: 0 = hard, 1 = gaussian, 2 = cosine
                MaskTex = GenerateAlphaMask(Mask, Display);
                Mask.Rect = [0 0 Mask.Dim];                                      	% Mask size
            end
            
            %===================== CREATE BACKGROUND TEXTURE
            if exist('BorderTexture','var')~= 2
                Stim.Window = Mask.DestRect;
                Stim.Background = Display.Background;
                BorderTexture = BackgroundSquares(Display, Stim);
            end
            
        end
         
        %========================= DISPLAY STIMULI ========================
        function R_stimOn(varargin)
            params=varargin{1,1};
            Eye = params(1);                                % Which eye to present CFS to
            Target.Image = params(2);                     	% Which image to present
            CFS.On = 1;                                     % CFS on?
            Target.On = 1;                                  % Target image on?
            Target.Texture = ImageTextures(Target.Image);   
            Target.Direction = 90;                          % Direction of target motion
         	CFS.Rect = [0 0 CFS.FrameSize];
            CFS.DestRect = CenterRect(CFS.Rect, Display.Rect);

            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            while 1
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    fprintf(strcat(commandIn,'\n'));
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn,'R_stimOff')            % Exits loop if trial is quit or ended
                    Exit=1;
                    break
                elseif strcmp(commandIn,'R_CFSOn') 
                    CFS.On = 1;
             	elseif strcmp(commandIn,'R_ImageOn') 
                    Target.On = 1;
                end
                
                %======== CFS eye
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);
                Screen('DrawTexture', Display.win, BorderTexture);                                  % Draw background texture
                if CFS.On == 1
                    Screen('DrawTexture', Display.win, CFSTextures(f), CFS.Rect, CFS.DestRect);     % Draw CFS stimulus
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);          % Apply Gaussian aperture mask
                    Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{E+1});
                else
                     Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{E+1});
                end
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);            % Draw fixation marker
                
                %======= Image eye
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, double(~Eye));
                Screen('DrawTexture', Display.win, BorderTexture);                                  % Draw background texture
                if Target.On == 1
                	Screen('DrawTexture', Display.win, Target.Texture, Target.Rect, Target.DestRect);
                    Screen('DrawTexture', Display.win, MaskTex, Mask.Rect, Mask.DestRect);          % Apply Gaussian aperture mask
                    Screen('FillOval', Display.win, Photodiode.OnColour, Photodiode.Rect{E+1});
                else
                	Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{E+1});
                end
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);        % Draw fixation marker
                [VBL FrameOnset] = Screen('Flip',Display.win);                              	% Flip window to screen
                [keyIsDown,secs,keyCode] = KbCheck;                                             % Check keyboard for 'escape' press        
                if keyIsDown && keyCode(KbName('Escape'))                                       % Press Esc for abort
                    break;
                end
            end
            R_stimOff;
        end
        
        %=========================== CLEAR TEXTURE ========================
        function R_stimOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background);
                Screen('DrawTexture', Display.win, Fix.Texture); 
                Screen('FillOval', Display.win, Photodiode.OffColour, Photodiode.Rect{Eye});
            end
            [VBL ImageOff] = Screen('Flip', Display.win);                                       % Clear stimulus  
        end
        
    end


%% ============================ apm_Retinotopy ============================
% Functional localizer stimuli for retinotopic mapping of visual areas.
% Polar and eccentricity maps are acheived using rotating wedge and
% expanding annulus stimuli respectively.
%
% INPUTS:
%   MapType:    1 = Polar mapping
%               2 = Eccentricity mapping
%   Direction:  1 = clockwise/expanding; 
%              -1 = anticlockwise/contracting
%   StimType:   1 = chromatic checkerboard
%               2 = 'buffy-otopy'
%               3 = binocular disparity
%   FixTask:    1 = respond to changes in fixation colour
%               2 = respond to changes in stimulus luminance
% 
% REFERENCES:
%   Swisher JD, Halko MA, Merabet LB, McMains SA & Somers DC (2007). Visual 
%       topography of human intraparietal sulcus. J.Neurosci.,
%       27(20):5326-5337.
%   Arcaro MJ, Pinsk MA, Li X & Kastner S (2011). Visuotopic organization
%       of macaque posterior parietal cortex: a functional magnetic resonance
%       imaging study. J.Neurosci., 31(6):2065-2078.
%==========================================================================

    function RetHandle = apm_Retinotopy()

        RetHandle = ...
            {...
                @R_punishStim;... 
                @R_fixOn;...
                @R_setFixColor;...
                @R_fixOff;...
                @R_stimOn;...
                @R_stimOff;...
                @R_makeStim;...
                @R_getBlock;...
            };
        
        persistent Wedge;
        persistent Annulus;
        persistent Stim;
        persistent Mask;
        persistent Fix;
        persistent StimTexture;
        persistent WinTexture;
        persistent MaskTexture;
        
        Fix.Type = 1;
        
     	%======================== PRESENT FIXATION ========================
        function R_fixOn (varargin)
            params=varargin{1,1};
            Fix.Offset = [params(1) params(2)]*Display.Pixels_per_deg(1);
            Fix.Size = params(3)*Display.Pixels_per_deg(1);
            Fix.Eye = 2;
            Fix = fixOn(Fix);
          	[VBL FixOn] = Screen('Flip', Display.win);
        end
        
        %====================== GET FIXATION COLOUR =======================
        function R_setFixColor (varargin)
            params=varargin{1,1};
            Fix.Colour = [params(1) params(2) params(3)]*255;
        end
        
        %========================== CLEAR FIXATION ========================
        function R_fixOff
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background, Fix.DestRect);
            end 
            [VBL FixOff] = Screen('Flip', Display.win);                                         % Clear stimulus  
        end
        
        %===================== PREPARE STIMULUS TEXTURES ==================
        function R_makeStim(varargin)
            params=varargin{1,1};
            Stim.Type = params(1);
        	Stim.Direction = params(2);
            Stim.OuterEcc = params(3);
            Wedge.Angle = params(4);
            
            %====== RADIAL CHECKERBOARD SETTINGS
            Stim.Chromatic = 1;                                     % 1 = chromatic; 0 = monochrome
            Stim.InnerEcc = round(0.5*Display.Pixels_per_deg(1)); 	% Innermost eccentricity
            Stim.OuterEcc = floor(Display.Rect(4)/Display.Pixels_per_deg(1)/2)*Display.Pixels_per_deg(1);    	% Outermost eccentricity (degrees)
            Stim.FlickerFrequency = 4;                              % Set flicker frequency (Hz) (Swisher et al. 2007 use 4Hz)
            Stim.PolarAngle = 10;                                   % Set polar angle (degrees) per division
            Stim.EccDivisions = 10;                                 % Set number of eccentricity divisions

            OuterDim = 2*Stim.OuterEcc;
            InnerDim = 2*Stim.InnerEcc;
            Stim.Rect = [0 0 OuterDim OuterDim];
            Stim.DestRect = CenterRect(Stim.Rect, Display.Rect);
            Logscale = abs(diff(log(Stim.EccDivisions:-1:1)));
            Stim.Eccentricities = ((Logscale*((Stim.OuterEcc-Stim.InnerEcc)/Logscale(end)))+Stim.InnerEcc);
            Stim.Eccentricities = Display.Pixels_per_deg(1):Display.Pixels_per_deg(1):Stim.OuterEcc;
            
            %============================= MAKE ALPHA MASK TEXTURES ===================
            WinTexture = repmat(Display.Background(1), ceil([OuterDim,OuterDim, 4]));   % create an RGBA matrix
            WinTexture(:,:,4) = 255;                                                 	% Set alpha to opaque
            Mask.Rect = [0 0 size(WinTexture,1),size(WinTexture,2)];
            Mask.DestRect = CenterRect(Mask.Rect, Display.Rect);
            Mask.OuterRect = CenterRect([0 0 OuterDim OuterDim],Mask.Rect);
            Mask.InnerRect = CenterRect([0 0 InnerDim InnerDim],Mask.Rect);

            if Stim.Type == 1
            	%====== POLAR MAPPING SETTINGS
                Wedge.Angle = 45;                                                       % Angular width of wedge (degrees)
                Wedge.AngularVelocity = 9;                                              % rotation speed (degrees per second) (Swisher et al. 2007 use 9)
                Wedge.CurrentRotation = 0;%Stim.Direction*Wedge.Angle;
                Wedge.RotationPerFrame = Wedge.AngularVelocity/Display.RefreshRate;
                MaskTexture = Screen('MakeTexture', Display.win, WinTexture, Wedge.CurrentRotation);            % make an opaque texture
                Screen('FillArc',MaskTexture,[0 0 0 0],Mask.OuterRect,Wedge.CurrentRotation,Wedge.Angle);       % draw a transparent wedge
                Screen('FillOval',MaskTexture,[Display.Background,255],Mask.InnerRect);                         % fill the centre
                
            elseif Stim.Type == 2
                %====== ECCENTRICITY MAPPING SETTINGS
                Wedge.CurrentRotation = 0;
                Annulus.Width = 1*Display.Pixels_per_deg(1);                            % width of annulus (degrees)
                Annulus.WidthPerFrame = 0.02*Display.Pixels_per_deg(1);
                Annulus.WidthMin = 1*Display.Pixels_per_deg(1);
                Annulus.WidthMax = 2*Display.Pixels_per_deg(1);
                Annulus.Acceleration = 9;                                               % expansion speed (degrees per second^2)
                Annulus.Speed = 2*Display.Pixels_per_deg(1);
                Annulus.CurrentEcc = Stim.InnerEcc;
                Annulus.EccPerFrame = Stim.Direction*Annulus.Speed/Display.RefreshRate;
            end

            %==================== MAKE RADIAL CHECKERBOARD TEXTURES ===================
            Stim.Angles = (0:Stim.PolarAngle:360);%-(Stim.PolarAngle/2);
            SegmentLum = {-1,1};
            Checker{1} = [1 2;2 1];
            Checker{2} = [2 1;1 2];
            for t = 1:2
                StimTexture(t) = Screen('MakeTexture', Display.win, WinTexture);% Create a blank texture
                for e = numel(Stim.Eccentricities):-1:1                         % For each eccentricity (in descending order)
                    for p = 1:1:numel(Stim.Angles)                          	% For each polar angle
                        s = Checker{t}(mod(e,2)+1,mod(p,2)+1);                  % Find segment brightness ('dark' or 'light')
                        if Stim.Chromatic == 1
                            SegmentColor = ceil(rand(1,3)*255);                 % Generate random RGB
                            SegmentColor = SegmentColor+(50*SegmentLum{s});     % Change brightness ('dark' or 'light')
                        else
                            SegmentColor = max([0 0 0],SegmentLum{s}*255);      % Black or white
                        end
                        Xdim = 2*Stim.Eccentricities(e);
                        StimEccRect = CenterRect([0 0 Xdim Stim.Eccentricities(e)*2],Mask.Rect);   
                        Screen('FillArc',StimTexture(t),SegmentColor,StimEccRect,Stim.Angles(p),360-(p*Stim.PolarAngle));
                    end
                end
            end

            
        end
        
        %======================= PRESENT STIMULUS =========================
        function R_stimOn(varargin)
            params=varargin{1,1};
            Eye = params(1);
            Stim.Cycles = 4;
            
            Cycle = 1;
            StartTime = GetSecs;
            pnet(con,'setreadtimeout',.001);
            pnet(con,'printf', reply);
            while Cycle <= Stim.Cycles
                commandIn = pnet(con,'readline');           % Check to make sure trial is still in progress
                if ~isempty(commandIn)
                    fprintf(strcat(commandIn,'\n'));
                    pnet(con,'printf', reply);              % return control to the socket
                end
                if strcmp(commandIn,'R_stimOff')            % Exits loop if trial is quit or ended
                    break
                end

                %================== DRAW STIMULUS
                t = mod(floor((GetSecs-StartTime)/(1/Stim.FlickerFrequency)),2)+1;
                
                %================== PREPARE APPROPRIATE ALPHA MASK
                if Stim.Type == 1
                    Wedge.CurrentRotation = Wedge.CurrentRotation+(Wedge.RotationPerFrame*Stim.Direction);                     % rotate wedge
                    Cycle = ceil(Wedge.CurrentRotation/360);
                elseif Stim.Type == 2
                    Annulus.Width = Annulus.Width+Annulus.WidthPerFrame;                                      % Set current annulus width
                    Annulus.OuterRect = [0 0 Mask.InnerRect(3)+(Annulus.Width*2) Mask.InnerRect(4)+(Annulus.Width*2)];  %
                    Annulus.InnerRect = [0 0 Mask.InnerRect([3,4])];
                    AnnulusOuterRect = Annulus.OuterRect+[0 0 Annulus.CurrentEcc Annulus.CurrentEcc];
                    AnnulusOuterRect = CenterRect(AnnulusOuterRect,Stim.Rect);
                    AnnulusInnerRect = Annulus.InnerRect+[0 0 Annulus.CurrentEcc Annulus.CurrentEcc];
                    AnnulusInnerRect = CenterRect(AnnulusInnerRect,Stim.Rect);

                    MaskTexture = Screen('MakeTexture', Display.win, WinTexture, 0);
                    Screen('FillOval',MaskTexture,[0 0 0 0],AnnulusOuterRect);
                    Screen('FillOval',MaskTexture,[Display.Background,255],AnnulusInnerRect);
                    Annulus.CurrentEcc = Annulus.CurrentEcc+Annulus.EccPerFrame;
                    if Annulus.CurrentEcc <= Stim.InnerEcc
                        Annulus.CurrentEcc =  Stim.OuterEcc;
                        Annulus.Width = Annulus.WidthMax;
                        Cycle = Cycle+1;
                    elseif Annulus.CurrentEcc >= Stim.OuterEcc
                        Annulus.CurrentEcc =  Stim.InnerEcc;
                        Annulus.Width = Annulus.WidthMin;
                        Cycle = Cycle+1;
                    end
                end
                
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                    Screen('DrawTexture', Display.win, StimTexture(t), Stim.Rect, Stim.DestRect);   % Draw stimulus
                    Screen('DrawTexture', Display.win, MaskTexture, Mask.Rect, Mask.DestRect, Wedge.CurrentRotation);    	% Draw alpha channel mask
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);        % Draw fixation marker
                end
                [VBL FrameOnset] = Screen('Flip',Display.win);                              	% Flip window to screen
                [keyIsDown,secs,keyCode] = KbCheck;                                             % Check keyboard for 'escape' press        
                if keyIsDown && keyCode(KbName('Escape'))                                       % Press Esc for abort
                    break;
                end
            end
            R_stimOff;
        end
        
        %======================= CLEAR STIMULUS =========================
        function R_stimOff(varargin)
            for Eye = 1:2
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);
                Screen('FillRect', Display.win, Display.Background);
            end 
            [VBL FixOff] = Screen('Flip', Display.win); 
        end
        
    end


%% ======================== QPCS SUBFUNCTIONS =============================
% The following functions can be called by all experiment functions above,
% and use and modify persistent variables specified at the top of this
% script.

%% ========================== PENALIZE ERROR ==============================
% If the monkey makes a task error (e.g. breaks fixation or fails to 
% saccade to the correct location on time, then s/he will be penalized 
% by terminating the current trial with an annoying flash (colour
% specified by Display.FlashBackground), optional tone or white noise,
% and a delay before the next trial commences.
%==========================================================================
    function R_punishStim()
        if Penalty.On == 1
            if Audio.On == 1
                if Audio.PsychSound == 1
                    PsychPortAudio('Start', Audio.Penalty, 1);
                else
                    Snd('Play',Audio.Penalty);
                end
            end
            for n = 1:Penalty.NoFlashes
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);     	% Select which eye to draw to
                    Screen('FillRect', Display.win, Display.FlashBackground);
                    Screen('DrawTexture',Display.win, Penalty.BadMonkey);
                end
                [VBL PenaltyOn] = Screen('Flip', Display.win);
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);     	% Select which eye to draw to
                    Screen('FillRect', Display.win, Display.Background(1));
                    Screen('DrawTexture',Display.win, Penalty.BadMonkey);
                end
                Screen('Flip', Display.win);
            end
            Screen('Flip', Display.win);
        end
    end

%% ============================ FIXATION ON ===============================
% Turn on fixation marker in the absence of a current onscreen stimulus.
% The type of fixation marker and other fixation marker parameters are 
% specified by fields belonging to the global structure 'Fix'.
% INPUT:
%       Fix.Colour:  RGB values
%       Fix.Size:    Diameter of fixation marker (pixels)
%       Fix.Offset:  Position of fixation from centre of screen (pixels)
%       Fix.Type:    0 = dot; 1 = cross; 2 = square; 3 = crosshairs
%       Fix.Eye:     0 = Right; 1 = Left; 2 = both
%==========================================================================
    function [Fix] = fixOn(Fix)
        if Fix.Size > 0
            Fix.Rect       	= ceil([0 0 Fix.Size Fix.Size]);
            Fix.DestRect   	= CenterRect(Fix.Rect, Display.Rect)+[Fix.Offset Fix.Offset];
            Fix.Background 	= ones(ceil(Fix.Size), ceil(Fix.Size), 4)*Display.Background(1);
            Fix.Background(:,:,4) = 0;
            if ~isfield(Fix, 'LineWidth'), Fix.LineWidth = 3; end
            Fix.Texture = Screen('MakeTexture', Display.win, Fix.Background);
            if Fix.Type==3
                Fix.Rotation = [0, 180];                                        % Binocular presentation of crosshairs
            else
                Fix.Rotation = [0, 0];                 
            end                             
            switch Fix.Type
                case 0                                                          %============== FILLED CIRCLE
                    Screen('FillOval', Fix.Texture, [Fix.Colour, 255], Fix.Rect);
                case 1                                                          %============== CROSS
                    Fix.Pos = [Fix.Size/2, -Fix.Size/2, 0, 0; 0, 0, Fix.Size/2, -Fix.Size/2];
                    Screen('FillOval', Fix.Texture, [Display.Background 255], Fix.Rect);
                    Screen('DrawLines', Fix.Texture, Fix.Pos, Fix.LineWidth, Fix.Colour, [Fix.Size/2, Fix.Size/2]);
                case 2                                                          %============== SOLID SQUARE
                    Screen('FillRect', Fix.Texture, Fix.Colour, Fix.Rect);
                case 3                                                          %============== BINOCULAR CROSSHAIRS
                    Fix.Pos = [Fix.Size/2, Fix.Size/4, 0, 0; 0, 0, Fix.Size/2, Fix.Size/4];
                    Screen('DrawLines', Fix.Texture, Fix.Pos, Fix.LineWidth, Fix.Colour, [Fix.Size/2, Fix.Size/2]);
                    Screen('FrameRect', Fix.Texture, Fix.Colour, CenterRect([0 0 Fix.Size,Fix.Size]/2,Fix.Rect), Fix.LineWidth);
            end
            %===================== DRAW TO WINDOW
            if Fix.Eye < 2                                                                      % For monocular fixation...
                Eye = abs(floor(Fix.Eye-0.5));                                                  % 0 = Right, 1 = Left
                currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye);             % Select which eye to draw to
                Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect);         
            elseif Fix.Eye == 2                                                                 % For binocular fixation...
                for Eye = 1:2
                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);     	% Select which eye to draw to
                    Screen('DrawTexture', Display.win, Fix.Texture, Fix.Rect, Fix.DestRect, Fix.Rotation(Eye)); 
                end
            end
        end
    end

%% ============================ FIXATION OFF ==============================
% Turn off fixation marker in the absence of a current onscreen stimulus.
% The type of fixation marker is specified by the global variable Fix.type.
%==========================================================================


%% ============================ LOAD IMAGES ===============================
% Loads all images from specified directory and returns a vector containing
% their texture handles ready for presentation.
% INPUTS:
%   Image.Dir:      full path of folder containing images
%   Image.No:       scalar or vector of image index/indices to load
%   Image.Radius:   radius in pixels
%==========================================================================
    function [ImageTextures, Image] = LoadImages(Image)
        if ~isfield(Image,'Radius')
            Image.Radius = 5*Display.Pixels_per_deg(1);                               	% Default radius of image size BEFORE drawing based on size input
        end
        Image.TextureWindow = [2*ImageRadius, 2*ImageRadius];                           % Set size of the stimulus before mask is applied (degrees)
        cd(Image.Dir);                                                                	% Change to image folder
        Images = dir('*.jpg');                                                          % Get all jpg images
        reverseStr = '';
        for n = Image.No
            LoadingText = sprintf('Loading image %d of %d...',n,numel(Image.No));   	% Print loading progress to command window
            fprintf([reverseStr, LoadingText]);
            reverseStr = repmat(sprintf('\b'), 1, length(LoadingText));
            if Debug.On == 1
                DrawFormattedText(Display.win, LoadingText, 20, Display.Rect(4)-40, [0 0 0], []);% Draw text
                Screen('Flip',Display.win);
            end
            [Img ImgCmap] = imread(fullfile(Image.Dir, Images(n).name));                % Load image file
            ImgDim = size(Img);                                                         % Get image dimensions
            if min(ImgDim([1,2])) ~= 2*Image.Radius                                    	% if smallest image dimesnion is not requested size...
                scale = 2*Image.Radius/min(ImgDim([1 2]));                              % Resize image so smallest dimension fits
                Img = imresize(Img, scale);
                ImgDim = size(Img);                                                     % Get new Img image dimensions
            end         
            if max(ImgDim([1,2])) > 2*Image.Radius                                   	% If largest image dimension is too large...
                Crop = (max(ImgDim([1,2]))-(2*Image.Radius))/2;
                if find(ImgDim==max(ImgDim([1,2])))==1
                    Img(end-Crop:end,:,:) = [];
                    Img(1:Crop,:,:) = [];
                elseif find(ImgDim==max(ImgDim([1,2])))==2
                    Img(:,end-Crop:end,:) = [];
                    Img(:,1:Crop,:) = [];
                end
            end
            Image.ImageRect = [0,0,size(Img,2),size(Img,1)];
            ImageTextures(n) = Screen('MakeTexture', Display.win, double(Img));          % Convert image to PTB texture
        end
    end

%% ============================= PLAY AUDIO ===============================
% Check whether audio is on, and which sound driver is in use.
%==========================================================================
    function [] = PlayAudio(Audio, Sound)
        if Audio.On == 1                            
            if Audio.PsychSound == 1
                PsychPortAudio('Start', Sound, 1);
            else
                Snd('Play',Sound);
            end
        end
    end

%% ============================= PLAY EPI ===============================
% Play MRI scanner sound recordings
%==========================================================================
    function [] = PlayEPI(EPIOn)
        if EPIOn == 1         
            EPI.Duration = 5*60;                                                                % Specify maximum playback duration (seconds)
            EPI.SoundDir = 'Z:\UCNI_Sounds\magnet sounds';                        
            [y, fs, nb] = wavread(fullfile(EPI.SoundDir, '4_EPI_1_11025.wav'));                 % Load EPI wave sounds files
            [y2, fs2, nb2] = wavread(fullfile(EPI.SoundDir, '5_EPI_2_11025.wav'));
            EPI.loop = [repmat(y',[1,3]),repmat(y2',[1,3])];                                    % Composite loop from fast and slow EPI waves
            EPI.Reps = ceil(EPI.Duration*fs/numel(EPI.loop));                                   % Find out how many times to play the loops
            EPI.PlayerObj = audioplayer(repmat(EPI.loop,[1,10]),fs);                            % Create an audio player object
            play(EPI.PlayerObj);                                                                % Begin playback

        elseif EPIOn == 0
            if isfield(EPI,'PlayerObj')
                stop(EPI.PlayerObj);                                                            % Stop playback
            end
        end
    end

%% ============================= PLAY METRONOME ===========================
% Play metronome audio
%==========================================================================
    function PlayMetronome(MetronomeOn)
        if MetronomeOn == 1         
            Metronome.Duration = 5*60;                                                                % Specify maximum playback duration (seconds)
            Metronome.SoundDir = 'Z:\UCNI_Sounds\magnet sounds';                        
            [y, fs, nb] = wavread(fullfile(Metronome.SoundDir, '4_EPI_1_11025.wav'));                 % Load Metronome wave sounds files
            [y2, fs2, nb2] = wavread(fullfile(Metronome.SoundDir, '5_EPI_2_11025.wav'));
            Metronome.loop = [repmat(y',[1,3]),repmat(y2',[1,3])];                                    % Composite loop from fast and slow EPI waves
            Metronome.Reps = ceil(Metronome.Duration*fs/numel(Metronome.loop));                                   % Find out how many times to play the loops
            Metronome.PlayerObj = audioplayer(repmat(Metronome.loop,[1,10]),fs);                            % Create an audio player object
            play(Metronome.PlayerObj);                                                                % Begin playback

        elseif MetronomeOn == 0
            if isfield(Metronome,'PlayerObj')
                stop(Metronome.PlayerObj);                                                            % Stop playback
            end
        end
    end


%% ======================= GET VBL ============================
    function [monitorFlipInterval] = GetDisplayTiming(win)
        nrSamples = 50;             % How many samples to take
        stddev = 0.0005;            % Standard deviation should be below this (s)
        timeout = 10;               % Timeout if not completed in this many seconds
        [monitorFlipInterval nrValidSamples stddev] =Screen('GetFlipInterval', win, nrSamples, stddev, timeout);
        
    end

%% ======================= GLOBAL OPTIONS MENU ============================
% Show icons for current settings, and allows changes to global settings
% via mouse pointer.
%==========================================================================
    function [] = ShowMenu
        for b = 1:numel(Button.IconFunctions)
            eval(sprintf('Button.On(b) = %s.On;',Button.IconFunctions{b}));
        end
        if ~isfield(Button,'LastPress')
            Button.LastPress = GetSecs;
        end
        ShowCursor;                                                                                             % make mouse pointer visible
        SetMouse(Button.DestRect(1,1), Button.DestRect(1,4), Display.win, 1);                                   % Set cursor location to near menu
        Button.Hover = zeros(1,numel(Button.IconFunctions));
        while 1
            Button.Change = 0;
         	currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, 0);                                   % darw to the left screen only
            for b = 1:numel(Button.DestRect(:,1))                                                               % For however many buttons...
                Screen('DrawTexture',Display.win,Button.Icon(Button.On(b)+1),[],Button.DestRect(b,:));          % Draw buttons
                Screen('DrawTexture',Display.win,Button.Highlight(Button.Hover(b)+1),[],Button.DestRect(b,:)); 	% Draw highlight if cursor is hovering
            end
            Screen('DrawTexture',Display.win,Audio.SpeakerIcon(Button.On(1)+1),[],Button.DestRect(1,:));        % Display audio settings icon
            for b = 2:numel(Button.DestRect(:,1))
                eval(sprintf('Screen(''DrawTexture'',Display.win,%s.Icon,[],Button.DestRect(b,:));', Button.IconFunctions{b}));
            end
            Screen('Flip', Display.win);                                                                        % Update the menu appearance
            [x y buttons] = GetMouse(Display.win, 1);                                                           % Get current cursor location
            for b = 1:numel(Button.DestRect(:,1))  
                Button.Hover(b) = IsInRect(x,y,Button.DestRect(b,:));
                if Button.Hover(b) && buttons(1) && GetSecs > Button.LastPress+Key.MinInterval
                    Button.LastPress = GetSecs;
                    Button.On(b) = double(~Button.On(b));
                    Button.Change = b;
                    eval(sprintf('%s.On = Button.On(b);',Button.IconFunctions{b}));
                end
            end
            
            %===================== UPDATE FUNCTIONS
            if Button.Change > 0
                switch Button.Change
                    case 1
                        if Audio.On == 1
                            Audio.a.Speak('Audio on.');
                        elseif Audio.On == 0
                           Audio.a.Speak('Audio off.'); 
                        end
                    case 3
                        if Photodiode.On == 0
                            Photodiode.OnColour =  Display.Background;
                            Photodiode.OffColour = Display.Background;              
                        elseif Photodiode.On == 1
                            Photodiode.OffColour = [128 128 128];
                            Photodiode.OnColour = [0 0 0];   
                        end
                    case 4
                        if GammaCorrect.On == 1
                            try
                                inverseCLUT = load(Display.CLUT, 'inverseCLUT');                                	% Load gamma lookup tables
                                Display.OriginalGamma = Screen('ReadNormalizedGammaTable', Display.ScreenID);
                                Screen('LoadNormalizedGammaTable', Display.ScreenID, inverseCLUT.inverseCLUT{2});	% Apply gamma table
                            catch
                                GammaCorrect.On = 0;
                            end
                        elseif GammaCorrect.On == 0
                            if isfield(Display, 'OriginalGamma')
                                Screen('LoadNormalizedGammaTable', Display.ScreenID, Display.OriginalGamma);
                            end
                        end
                    case 5
                        if Sleep.On == 0                                                                         	% Turn screen black so monkey can rest
                            for p = Display.Background:-1:0
                                for e = 1:2
                                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, e-1);         	% Draw to screen
                                    Screen('FillRect', Display.win, p);
                                end
                                Screen('Flip', Display.win);
                            end
                        elseif Sleep.On == 1                                                                     	% Turn screen grey ready for experiment
                            for p = 1:Display.Background
                                for e = 1:2
                                    currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, e-1);           	% Draw to screen
                                    Screen('FillRect', Display.win, p);
                                end
                                Screen('Flip', Display.win);
                            end
                        end
                    case 6
                        
                    case 7
                        PlayEPI(EPI.On);
                    case 8
                        PlayMetronome(Metronome.On);
                    case 9
                        break;
                end
            end
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();                                                      % Check keyboard input for command to toggle menu off
            if keyIsDown && keyCode(Key.Menu) && secs > Key.LastPress+Key.MinInterval+0.4 || keyCode(Key.Exit)
                Key.LastPress = secs;
                break;
            end
        end                                                                             
     	Screen('Flip', Display.win);                                                                                % Clear menu
        HideCursor;                                                                                                 % Hide mouse pointer
    end

end
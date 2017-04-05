

%==================== MF3D_AnimatedInteractionDemo.m ======================
% This script runs a simple demo of how multiple short movie clips/ image 
% sequences can be presented in an interactive closed loop using standard
% PsychToolbox functions.
%
%
%==========================================================================



%============ Prepare PTB window
Stereo              = 0;                            % Stereoscopic presentation?
Display             = DisplaySettings(Stereo);      % Get display settings
Display.Rect        = Display.Rect/2;               % For debugging...
Display.Background  = [0 0 0];                      % Set background color (RGB)
Display.Imagingmode = [];                           
%HideCursor;
KbName('UnifyKeyNames');
Screen('Preference', 'VisualDebugLevel', 1);     
[Display.win, Display.Rect] = Screen('OpenWindow', Display.ScreenID, Display.Background,Display.Rect,[],[], Display.Stereomode, [], Display.Imagingmode);


%============= Load animation data
AnimationDir    = '/Volumes/Seagate Backup 1/NIH_PhD_nonthesis/7. 3DMacaqueFaces/AnimationDemos';
FileFormat      = 'avi';
ClipFiles       = wildcardsearch(AnimationDir, sprintf('*.%s',FileFormat));

Clips           = MF3D_LoadClips(Display.win, ClipFiles);

SourceRect      = Display.Rect;
DestRect        = Display.Rect;


MaskFile = '/Volumes/Seagate Backup 1/NIH_PhD_nonthesis/7. 3DMacaqueFaces/AnimationDemos/TargetMask.png';
Im = imread(MaskFile);


%============= Begin demo
while 1
    %============== Check gaze/ mouse location
    [x,y] = GetMouse(Display.win);
    
    
    Screen('PlayMovie',mov,1,[],Movie.Volume);
    Screen('SetmovieTimeIndex',mov,StartTime,1); 
    MovieTex = Screen('GetMovieImage', Display.win, mov, 1);
    
    
    %============== Draw frame
    for Eye = 1:2                                                                   % For each eye...
        currentbuffer = Screen('SelectStereoDrawBuffer', Display.win, Eye-1);       % Select buffer for each eye
        Screen('DrawTexture', Display.win, MovieTex, SourceRect, DestRect);         % Draw frame for current eye
        if DrawMousePos == 1
            
        end
    end
    [VBL FrameOnset(end+1)] = Screen('Flip', Display.win);                          % Flip frame to screen
    
    %============== Check user input
    [keyIsDown,secs,keyCode] = KbCheck;                                             % Check keyboard for 'escape' press        
    if keyIsDown && keyCode(KbName('Escape')) == 1                                  % Press Esc for abort
        break
    end
    
    
end

%============== Perform clean up
Screen('CloseAll');      
ShowCursor;

%============== Show some stats about software performance
Frametimes      = diff(FrameOnset);
meanFrameRate   = mean(Frametimes(2:end))*1000;
semFrameRate    = (std(Frametimes(2:end))*1000)/sqrt(numel(Frametimes(2:end)));
fprintf('Frames shown............%.0f\n', numel(Frametimes));
fprintf('Mean frame duration.....%.0f ms +/- %.0f ms\n', meanFrameRate, semFrameRate);
fprintf('Max frame duration......%.0f ms\n', max(Frametimes)*1000);



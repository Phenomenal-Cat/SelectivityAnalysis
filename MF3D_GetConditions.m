function Params = MF3D_GetConditions(ExpType)

%======================== MF3D_GetConditions.m ============================
% Reconstruct which stimulus numbers correspond to which conditions for
% each of the different experiment types ('ExpType') in the apm_StereoFaces
% pilot experiments. Input experiment types are:
%       1) Macaque face, single identity, neutral, 3x7 angles
%       2) Macaque face, single identity, fear grimace, 3x7 angles
%       3) Macaque face, single identity, 4 expressions, 3x3 angles, convex & concave
%       4) Macaque face, 5 identities, neutral, 3x3 angles
%       5) Human face, 2 identities, 1 expression, 3x5 angles
%       
%
%==========================================================================

if ExpType == 5
    FaceImageDir    = 'P:\murphya\Stimuli\BlenderRenders\HumanFaces\';
elseif ExpType == 7
    FaceImageDir    = 'P:\murphya\MacaqueFace3D\GazeExperiments\Renders\Experiment7_v2\';
else
    FaceImageDir    = 'P:\murphya\MacaqueFace3D\BlenderFiles\Renders\';

end


if ExpType == 4
    Params.MonkeyIDs 	= [1,2,3,4,5];
else
    Params.MonkeyIDs   = 1;
end

Params.Species     = 'Macaque';
Params.Elevations  = [-30, 0, 30];
Params.Distances   = [-20, 0, 20];
    
switch ExpType
    case 1
        Params.Expressions = {'neutral'};
        Params.Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
        Params.Scales      = [10, 12, 20];

    case 2
        Params.Expressions = {'fear'};
        Params.Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
        Params.Scales      = [10, 12, 20];
        
    case 3
        Params.Expressions = {'neutral','fear','lipsmack','threat'};
        Params.Azimuths    = [-30, 0, 30];
        Params.Scales      = [12];
        
    case 4
        Params.Expressions = {'neutral'};
        Params.Azimuths    = [-30, 0, 30];
        Params.Scales      = [12];
        
    case 5
        Params.Species     = 'Human';
        Params.Expressions = {'neutral'};
        Params.Azimuths    = [-60,-30, 0, 30,60];
        Params.Scales      = [60, 80];                 % Only 2 scales were shown, although a third scale (120%) was generated
        Params.MonkeyIDs   = {'Debbie','Ewan'};
        
    case 6                                                      % Added 04/22/17 for new Spice data
       	Params.Expressions = {'neutral','fear','lipsmack','threat'}; 
        Params.Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
        Params.Scales      = [10, 12, 20];
        Params.Distances   = 0;                             	% Only show 1 distance (since presentation is mono)
        
  	case 7
        Params.Expressions = {'neutral'};
        Params.Azimuths    = [-30,-15, 0, -15, 30];
        Params.Elevations  = [-15, 0, 15];
        Params.Distances   = 0; 
        
    otherwise
        error('ExpType %d is not recognized!', ExpType);
end

Indx        = 1;
for exp = 1:numel(Params.Expressions)
    for az = 1:numel(Params.Azimuths)
        for el = 1:numel(Params.Elevations)
            if ExpType == 7
                for gaz = 1:numel(Params.Azimuths)
                    for gel = 1:numel(Params.Elevations)
                        Pic.ImgFilenames{Indx} = fullfile(FaceImageDir, sprintf('MacaqueGaze_Neutral_Haz%d_Hel%d_Gaz%d_Gel%d_dist0.png', Params.Azimuths(az), Params.Elevations(el), Params.Azimuths(gaz), Params.Elevations(gel)));
                      	Pic.ConditionMatrix(Indx,:) = [az, el, gaz, gel];
                        Indx = Indx+1;
                    end
                end
            elseif ExpType <7
                for d = 1:numel(Params.Distances)
                    for s = 1:numel(Params.Scales)
                        for m = 1:numel(Params.MonkeyIDs)
                            switch Params.Species
                                case 'Macaque'
                                    Params.Filenames{Indx} = sprintf('Macaque_Id%d_%s_az%d_el%d_dist%d_sc%d.png',  Params.MonkeyIDs(m), Params.Expressions{exp}, Params.Azimuths(az), Params.Elevations(el), Params.Distances(d), Params.Scales(s));
                                    Params.FullFilenames{Indx} = fullfile(FaceImageDir, sprintf('Monkey_%d', Params.MonkeyIDs(m)), Params.Filenames{Indx});
                                case 'Human'
                                    Params.Filenames{Indx} = sprintf('Human_%s_az%d_el%d_dist%d_sc%d.png',  Params.MonkeyIDs{m}, Params.Azimuths(az), Params.Elevations(el), Params.Distances(d), Params.Scales(s));
                                    Params.FullFilenames{Indx} = fullfile(FaceImageDir, Params.Filenames{Indx});
                            end
                            if numel(Params.MonkeyIDs) == 1
                                Params.ConditionMatrix(Indx,:) = [exp, az, el, d, s];
                            elseif numel(Params.MonkeyIDs)>1
                                Params.ConditionMatrix(Indx,:) = [m, az, el, d, s];
                            end
                            Indx = Indx+1;
                        end
                    end
                end
            end
        end
    end
end
% Include stereoscopic depth profile as a variable? (1) stereo congruent; 0) mono; -1) incongruent:
if ExpType == 3
    Params.ConditionMatrix = [Params.ConditionMatrix, ones(size(Params.ConditionMatrix,1),1); Params.ConditionMatrix, zeros(size(Params.ConditionMatrix,1),1); Params.ConditionMatrix, -1*ones(size(Params.ConditionMatrix,1),1)];
end
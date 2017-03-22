function Params = GetConditions(ExpType)

%=========================== GetConditions.m ==============================
% Reconstruct which stimulus numbers correspond to which conditions for
% each of the different experiment types ('ExpType') in the apm_StereoFaces
% pilot experiments.
%==========================================================================

if nargin ==0
    ExpType = 1;
end
Params.ExpType = ExpType;

[~,CompName] = system('hostname');  
if strcmpi(CompName(1:end-1), 'Aidans-MacBook-Pro.local')
    FaceImageDir    = 'P:\murphya\MacaqueFace3D\BlenderFiles\Renders\';
else
    FaceImageDir    = 'P:\murphya\MacaqueFace3D\BlenderFiles\Renders\';
end

if ExpType < 4
    Params.MonkeyIDs    = 1;
    Params.Factors      = {'Elevations','Azimuths','Distances','Scales','Expressions'};     % All factors tested
    Params.CondMatCol 	= [3, 2, 4, 5, 1];                                                 	% Which column is each factor coded in?
elseif ExpType >= 4
    Params.MonkeyIDs 	= [1,2,3,4,5];
    Params.Factors      = {'Elevations','Azimuths','Distances','Scales','MonkeyID'};        % All factors tested
    Params.CondMatCol 	= [3, 2, 4, 5, 1];                                                 	% Which column is each factor coded in?
end

Params.Elevations  = [-30, 0, 30];
Params.Distances   = [-20, 0, 20];
if ExpType == 1
    Params.Expressions  = {'neutral'};
    Params.Azimuths     = [-90, -60, -30, 0, 30, 60, 90];
    Params.Scales       = [10, 12, 20];
elseif ExpType == 2
    Params.Expressions = {'fear'};
    Params.Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
    Params.Scales      = [10, 12, 20];
elseif ExpType == 3
    Params.Expressions = {'neutral','fear','lipsmack','threat'};
    Params.Azimuths    = [-30, 0, 30];
    Params.Scales      = [12];
elseif ExpType == 4
    Params.Expressions = {'neutral'};
    Params.Azimuths    = [-30, 0, 30];
    Params.Scales      = [12];
end

Indx        = 1;
for exp = 1:numel(Params.Expressions)
    for az = 1:numel(Params.Azimuths)
        for el = 1:numel(Params.Elevations)
            for d = 1:numel(Params.Distances)
                for s = 1:numel(Params.Scales)
                    for m = 1:numel(Params.MonkeyIDs)
                        Params.Filenames{Indx} = sprintf('Macaque_Id%d_%s_az%d_el%d_dist%d_sc%d.png',  Params.MonkeyIDs(m), Params.Expressions{exp}, Params.Azimuths(az), Params.Elevations(el), Params.Distances(d), Params.Scales(s));
                        Params.FullFilenames{Indx} = fullfile(FaceImageDir, sprintf('Monkey_%d', Params.MonkeyIDs(m)), Params.Filenames{Indx});
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
% Include stereoscopic depth profile as a variable? (1) stereo congruent; 0) mono; -1) inconruent:
if ExpType == 3
    Params.ConditionMatrix = [Params.ConditionMatrix, ones(size(Params.ConditionMatrix,1),1); Params.ConditionMatrix, zeros(size(Params.ConditionMatrix,1),1); Params.ConditionMatrix, -1*ones(size(Params.ConditionMatrix,1),1)];
end
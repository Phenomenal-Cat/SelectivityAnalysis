%===================== StereoFaceConditions.m =============================
% Generate a 'Params' structure to store information about the stimuli
% presented in different variants of the 'StereoFaces' experiment.

Params.ImageDir    = '/Volumes/PROJECTS/murphya/MacaqueFace3D/BlenderFiles/Renders/Monkey_1';
Params.Expressions = {'neutral','fear','lipsmack','threat'};
% Params.Expressions = {'Fear'};
% Params.Azimuths    = [-90, -60, -30, 0, 30, 60, 90];
Params.Azimuths    = [-30, 0, 30];
Params.Elevations  = [-30, 0, 30];
Params.Distances   = [-20, 0, 20];
% Params.Scales      = [10, 12, 20];
Params.Scales      = [12];
Params.Stereo     	= [1, 0, -1];

Params.Header     	= {'Expression','Azimuth','Elevation','Distance','Scale','Stereo'};
Indx                = 1;
for exp = 1:numel(Params.Expressions)
    for az = 1:numel(Params.Azimuths)
        for el = 1:numel(Params.Elevations)
            for d = 1:numel(Params.Distances)
                for s = 1:numel(Params.Scales)
                    Params.ImgFilenames{Indx} = fullfile(Params.ImageDir, sprintf('Macaque_%s_az%d_el%d_dist%d_sc%d.png', Params.Expressions{exp}, Params.Azimuths(az), Params.Elevations(el), Params.Distances(d), Params.Scales(s)));
                    Params.ConditionMatrix(Indx,:) = [exp, az, el, d, s];
                    Indx = Indx+1;
                end
            end
        end
    end
end

%<<<<<<<<< To include stereoscopic depth profile as a
%variable (1) stereo congruent; 0) mono; -1) inconruent:
StereoColumn = [];
for s = 1:numel(Params.Stereo)
    StereoColumn = [StereoColumn; repmat(Params.Stereo(s),[size(Params.ConditionMatrix,1),1])];
end
Params.ConditionMatrix = [repmat(Params.ConditionMatrix, [numel(Params.Stereo),1]), StereoColumn];
% Params.ConditionMatrix = [Params.ConditionMatrix, ones(size(Params.ConditionMatrix,1),1); Params.ConditionMatrix, zeros(size(Params.ConditionMatrix,1),1); Params.ConditionMatrix, -1*ones(size(Params.ConditionMatrix,1),1)];

OutputDir   = '/Volumes/procdata/murphya/Physio/StereoFaces';
Filename    = fullfile(OutputDir, 'StereoFaces_Conditions_20160617');
save(Filename, 'Params');
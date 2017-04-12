
if ismac
    Append = '/Volumes';
end
StimDir = fullfile(Append, 'projects/murphya/MacaqueFace3D/BlenderFiles/Renders/Monkey_1/');

%========== Plot stimulus images for Exp. 1/ 2
Subject         = 'Matcha';
Date            = '20160615';
TimingData      = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/',sprintf('StimTimes_%s_%s.mat', Subject, Date));
load(TimingData)
ImCrop = [650, 200];                                                        % Number of pixels to crop from each edge in X and Y
figure('position', [-1911, 168, 1799, 840]);
Axh = tight_subplot(3,7,0, 0,0);
AxIndx  = 1;
Dist    = 1;
Scale   = 1;
for el = 1:numel(Params.Elevations)
    for az = 1:numel(Params.Azimuths)
        Cond = find(ismember(Params.ConditionMatrix, [1, az, el, Dist, Scale],'rows'));
        axes(Axh(AxIndx));
        [Im, cm, Alpha] = imread(fullfile(StimDir, Params.Filenames{Cond}));
        Imsize      = size(Im);
        Im          = Im(ImCrop(2):(Imsize(1)-ImCrop(2)), ImCrop(1):((Imsize(2)/2)-ImCrop(1)),:);
        Alpha       = Alpha(ImCrop(2):(Imsize(1)-ImCrop(2)), ImCrop(1):((Imsize(2)/2)-ImCrop(1)));
        imh(Cond)   = image(Im);
        alpha(imh(Cond), Alpha);
        axis tight equal off
        AxIndx = AxIndx+1;
    end
end
% export_fig('/Volumes/projects/murphya/MacaqueFace3D/PilotData/PNGs/OrientStim_Fear.png', '-png','-transparent');


%========== Plot stimulus images for Exp. 3
Subject         = 'Avalanche';
Date            = '20160630';
TimingData      = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/',sprintf('StimTimes_%s_%s.mat', Subject, Date));
load(TimingData)
ImCrop = [650, 100];                                                        % Number of pixels to crop from each edge in X and Y
Dist    = 1;
Scale   = 1;

for exp = 1:numel(Params.Expressions)
    figure('position', [-1911, 168, 786, 840]);
    Axh = tight_subplot(3,3,0, 0,0);
    AxIndx  = 1;
    for el = 1:numel(Params.Elevations)
        for az = 1:numel(Params.Azimuths)
            Cond = find(ismember(Params.ConditionMatrix, [exp, az, el, Dist, Scale, 1],'rows'));
            axes(Axh(AxIndx));
            [Im, cm, Alpha] = imread(fullfile(StimDir, Params.Filenames{Cond}));
            Imsize      = size(Im);
            Im          = Im(ImCrop(2):(Imsize(1)-ImCrop(2)), ImCrop(1):((Imsize(2)/2)-ImCrop(1)),:);
            Alpha       = Alpha(ImCrop(2):(Imsize(1)-ImCrop(2)), ImCrop(1):((Imsize(2)/2)-ImCrop(1)));
            imh(Cond)   = image(Im);
            alpha(imh(Cond), Alpha);
            axis tight equal off
            AxIndx = AxIndx+1;
        end
    end
    export_fig(fullfile('/Volumes/projects/murphya/MacaqueFace3D/PilotData/PNGs/',sprintf('OrientStim3x3_%s.png',Params.Expressions{exp})), '-png','-transparent');
end


%========== Plot stimulus images for Exp. 4
Subject         = 'Avalanche';
Date            = '20160712';
StimDir         = fullfile(Append, 'projects/murphya/MacaqueFace3D/BlenderFiles/Renders/');
TimingData      = fullfile(Append, '/procdata/murphya/Physio/StereoFaces/Timing/StereoFaces/',sprintf('StimTimes_%s_%s.mat', Subject, Date));
load(TimingData)
ImCrop = [650, 100];                                                        % Number of pixels to crop from each edge in X and Y
Scale   = 1;
Dist = 1;

for id = 1:numel(Params.MonkeyIDs)
    figure('position', [-1911, 168, 786, 840]);
    Axh = tight_subplot(3,3,0, 0,0);
    AxIndx  = 1;
    for el = 1:numel(Params.Elevations)
        for az = 1:numel(Params.Azimuths)
            Cond = find(ismember(Params.ConditionMatrix, [id, az, el, Dist, Scale],'rows'));
            axes(Axh(AxIndx));
            [Im, cm, Alpha] = imread(fullfile(StimDir, sprintf('Monkey_%d', id), Params.Filenames{Cond}));
            Imsize      = size(Im);
            Im          = Im(ImCrop(2):(Imsize(1)-ImCrop(2)), ImCrop(1):((Imsize(2)/2)-ImCrop(1)),:);
            Alpha       = Alpha(ImCrop(2):(Imsize(1)-ImCrop(2)), ImCrop(1):((Imsize(2)/2)-ImCrop(1)));
            imh(Cond)   = image(Im);
            alpha(imh(Cond), Alpha);
            axis tight equal off
            AxIndx = AxIndx+1;
        end
    end
    export_fig(fullfile('/Volumes/projects/murphya/MacaqueFace3D/PilotData/PNGs/',sprintf('MonkeyID_%d.png',id)), '-png','-transparent');
end


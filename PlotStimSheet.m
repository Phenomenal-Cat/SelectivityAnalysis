
%========================= PlotStimSheet.m ================================
% Plot a subset of stimulus images in a single figure for illustration
% purposes.
%
%==========================================================================


CondFile    = '/Volumes/projects/murphya/MacaqueFace3D/StereoFaceAnalysis/StereoFace_AllConditions_20160613.mat';
load(CondFile);

figure('position', [ -1915, 100, 1913, 890]);
ImIndices   = find(ismember(Pic.ConditionMatrix(:,4:5),[2,3],'rows'));
axh         = tight_subplot(3,7, 0, 0, 0);

Order = [3:3:21, 2:3:21, 1:3:21]+42;
% X = 1:size(im,2)/2;
X = 500:1500;
for p = 1:numel(Order)
    axes(axh(p));
    [im, map, alph] = imread(Pic.ImgFilenames{ImIndices(Order(p))});
    imh(p) = image(im(:,X,:));
    set(imh(p), 'AlphaData', alph(:,X));
    axis equal tight off
end
export_fig(Filename, '-png', '-transparent','-nocrop')




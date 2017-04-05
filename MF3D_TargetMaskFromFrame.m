function Mask = MF3D_TargetCoordsFromFrame(Frame, TargetVals)

%====================== MF3D_TargetCoordsFromFrame.m ======================
% If the screen coordinates (x,y pixels) of onscreen target objects are not
% available then this function can estimate the coordinates of each target
% center based on a single input frame that has been pre-thresholded.
%==========================================================================

if numel(size(Frame))>2
    Frame = double(rgb2gray(Frame));
end

axh(1) = subplot(2,2,1);
imagesc(Frame);
axis equal tight off

axh(2) = subplot(2,2,2);
hist(reshape(Frame(:,:), [1, numel(Frame(:,:))]), 100);
grid on;
xlabel('Hue value','fontsize', 16);
ylabel('Frequency','fontsize', 16);


[centers,radii] = imfindcircles(Frame,[30, 60])


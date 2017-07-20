




StimDir = '/Volumes/PROJECTS/murphya/MacaqueFace3D/GazeExperiments/Renders/Experiment7_v2';


Az  = -30:15:30;
El  = -15:15:15; 
Ypix = 250:800;
Xpix = 300:980;

% figure('position',get(0,'ScreenSize'))
% axh = tight_subplot(numel(El), numel(Az), 0, 0, 0);
% i = 1;
% for el = numel(El):-1:1
%     for az = 1:numel(Az)
%         Filename = fullfile(StimDir, sprintf('MacaqueGaze_Neutral_Haz%d_Hel%d_Gaz%d_Gel%d_dist0.png', 0, 0, Az(az), El(el)));
%         [im, cm, alph] = imread(Filename);
%         
%         im = im(Ypix, Xpix, :);
%         alph = alph(Ypix, Xpix);
%         axes(axh(i));
%         imh(i) = image(im);
%         axis equal tight off
%         alpha(imh(i), alph);
%         i = i+1;
%     end
% end
% % export_fig('GazeStim_EyePos.png','-png','-transparent')


Ypix = [220:770; 250:800; 280:830];
Xpix = 300:980;
Thresh = [0.3, 0.34];

figure('position',get(0,'ScreenSize'))
axh = tight_subplot(numel(El), numel(Az), 0, 0, 0);
i = 1;
for el = 1:numel(El)
    for az = 1:numel(Az)
        Filename = fullfile(StimDir, sprintf('MacaqueGaze_EyesOnly_Neutral_Haz%d_Hel%d_Gaz%d_Gel%d_dist0.png', Az(az), El(el), 0, 0));
        [im, cm, alph] = imread(Filename);
%         im2     = rgb2hsl(double(im)/255);
%         mask    = alph/max(alph(:));
%         mask(im2(:,:,1) > Thresh(1) & im2(:,:,1) < Thresh(2)) = 0;
%         mask    = imfill(mask);
%         alph    = mask*255;       
        im = im(Ypix(el,:), Xpix, :);
        alph = alph(Ypix(el,:), Xpix);
        axes(axh(i));
        imh(i) = image(im);
        axis equal tight off
        alpha(imh(i), alph);
        i = i+1;
    end
end
set(gcf,'position', [-1920 0 1920 1080])
% export_fig('GazeStim_HeadPos.png','-png','-transparent')



% %=================  Gaze directions
% figure('position',get(0,'ScreenSize'))
% axh = tight_subplot(numel(El), numel(Az), 0, 0, 0);
% i = 1;
% for el = 1:numel(El)
%     for az = 1:numel(Az)
%         Filename = fullfile(StimDir, sprintf('MacaqueGaze_Neutral_Haz%d_Hel%d_Gaz%d_Gel%d_dist0.png', Az(az), El(el), -Az(az), El(el)));
%         [im, cm, alph] = imread(Filename);
%         
%         im = im(Ypix(el,:), Xpix, :);
%         alph = alph(Ypix(el,:), Xpix);
%         axes(axh(i));
%         imh(i) = image(im);
%         axis equal tight off
%         alpha(imh(i), alph);
%         i = i+1;
%     end
% end
% set(gcf,'position', [-1920 0 1920 1080])

% %================= Eyes closed
% figure('position',get(0,'ScreenSize'))
% axh = tight_subplot(numel(El), numel(Az), 0, 0, 0);
% i = 1;
% for el = 1:numel(El)
%     for az = 1:numel(Az)
%         Filename = fullfile(StimDir, sprintf('MacaqueGaze_Eyes_Neutral_Haz%d_Hel%d_GazNaN_GelNaN_dist0.png', Az(az), El(el)));
%         [im, cm, alph] = imread(Filename);
%         
%         im = im(Ypix(el,:), Xpix, :);
%         alph = alph(Ypix(el,:), Xpix);
%         axes(axh(i));
%         imh(i) = image(im);
%         axis equal tight off
%         alpha(imh(i), alph);
%         i = i+1;
%     end
% end
% set(gcf,'position', [-1920 0 1920 1080])

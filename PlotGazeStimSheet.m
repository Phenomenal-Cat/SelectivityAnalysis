




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

BorderSize  = 200;
Offsets     = [0,0; BorderSize,BorderSize; -BorderSize, BorderSize; -BorderSize, -BorderSize; BorderSize, -BorderSize];



for o = 1%:size(Offsets, 1)
    figure('position',get(0,'ScreenSize'))
    axh = tight_subplot(numel(El), numel(Az), 0.01, 0.01, 0.01);
    i = 1;
    for el = 1:numel(El)
        for az = 1:numel(Az)
            Filename = fullfile(StimDir, sprintf('MacaqueGaze_EyesOnly_Neutral_Haz%d_Hel%d_Gaz%d_Gel%d_dist0.png', Az(az), El(el), 0, 0));
            %Filename = fullfile(StimDir, sprintf('MacaqueGaze_Neutral_Haz%d_Hel%d_Gaz%d_Gel%d_dist0.png', Az(az), El(el), 0, 0));
            [im, cm, alph] = imread(Filename);
    % %         im2     = rgb2hsl(double(im)/255);
    % %         mask    = alph/max(alph(:));
    % %         mask(im2(:,:,1) > Thresh(1) & im2(:,:,1) < Thresh(2)) = 0;
    % %         mask    = imfill(mask);
    % %         alph    = mask*255;       
    %         im = im(Ypix(el,:), Xpix, :);
    %         alph = alph(Ypix(el,:), Xpix);
            ImSize = size(im);
            im = im(BorderSize:(size(im,1)-BorderSize), BorderSize:(size(im,2)/2-BorderSize), :);
    
%             im = im(:,1:size(im,2)/2, :);
%             alph = alph(:,1:size(alph,2)/2);
            axes(axh(i));
            imh(i) = image(BorderSize+Offsets(o,1), BorderSize+Offsets(o,2), im);
%             alpha(imh(i), alph);
            grid on
            axis equal tight %off
            hold on;
            plot(repmat(ImSize(2)/4,[1,2]),repmat(ImSize(1)/2,[1,2]), '.r','markersize', 20);
            set(gca, 'xlim', [0, ImSize(2)/2], 'ylim', [0, ImSize(1)],'color', [0.5, 0.5, 0.5]);
            rectangle('position',[BorderSize,BorderSize,ImSize(2)/2-BorderSize*2,ImSize(1)-BorderSize*2]);
            drawnow
            i = i+1;
        end
    end
    set(gcf,'position', [-1920 0 1920 1080])
% export_fig('GazeStim_HeadPos.png','-png','-transparent')

end


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


%==================== PlotFingerPrintStim

StimDir = '/Volumes/PROJECTS/koyanok/Stimuli/Categories_Stim10k_small';
NoStim = 60;

figure('position',get(0,'Screensize'));
axh = tight_subplot(6,10,0.02,0.02,0.02);
ImageRadius = 300;
for s = 1:NoStim
    StimFile = fullfile(StimDir, sprintf('Cat%03d.png', s));
    [img] = imread(StimFile);
    ImgDim = size(img);                                                         % Get image dimensions
    if min(ImgDim([1,2])) ~= 2*ImageRadius                                    	% if smallest image dimesnion is not requested size...
        scale = 2*ImageRadius/min(ImgDim([1 2]));                               % Resize image so smallest dimension fits
        img = imresize(img, scale);
        ImgDim = size(img);                                                     % Get new Img image dimensions
    end         
    if max(ImgDim([1,2])) > 2*ImageRadius                                   	% If largest image dimension is too large...
        Crop = (max(ImgDim([1,2]))-(2*ImageRadius))/2;
        if find(ImgDim==max(ImgDim([1,2])))==1
            img(end-Crop:end,:,:) = [];
            img(1:Crop,:,:) = [];
        elseif find(ImgDim==max(ImgDim([1,2])))==2
            img(:,end-Crop:end,:) = [];
            img(:,1:Crop,:) = [];
        end
    end
    
    axes(axh(s));
    image(img);
    axis equal tight off
    title(sprintf('%d', s), 'fontsize', 14);
    drawnow
end
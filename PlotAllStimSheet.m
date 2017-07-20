
%========================= PlotAllStimSheet.m ================================
% Plot a subset of stimulus images in a single figure for illustration
% purposes.
%
%==========================================================================


ImageDir    = '/Volumes/PROJECTS/murphya/MacaqueFace3D/BlenderFiles/Renders/Monkey_0';
SaveDir     = '/Volumes/PROJECTS/murphya/MacaqueFace3D/Figures/';

SqueezeFramed = 1;
Exp = {'Neutral','Threat','Fear','Coo','Yawn'};
Az  = [-90, -60, -30, 0, 30, 60, 90];
El  = [-30, 0, 30];
OffsetPix   = 200;
X           = OffsetPix:((1920/2)-OffsetPix);
Y           = 50:1000;
FigPos      = [-1921, 379, 1804, 629];
%FigPos = [ -1915, 100, 1913, 890];

for exp = 1:numel(Exp)
    fh(exp)     = figure('position', FigPos);
    axh         = tight_subplot(numel(El), numel(Az), 0, 0, 0);

    Order = [3:3:21, 2:3:21, 1:3:21]+42;
    % X = 1:size(im,2)/2;
    
    for e = 1:numel(El)
        for a = 1:numel(Az)
            AxIndx = ((e-1)*numel(Az)) + a;
            ImgFilename = fullfile(ImageDir, sprintf('Macaque_HeadRot_%s_az%d_el%d_dist-20_sc100.png', Exp{exp}, 0-Az(a), El(e)));
            axes(axh(AxIndx));
            [im, map, alph] = imread(ImgFilename);
            if SqueezeFramed == 1
                Image   = imresize(im(Y,X,:), [round(numel(Y)/2), numel(X)]);
                Alph    = imresize(alph(Y,X), [round(numel(Y)/2), numel(X)]);
            end
            imh(AxIndx) = image(Image);
            set(imh(AxIndx), 'AlphaData', Alph);
            axis equal tight off
        end
    end
    Filename = fullfile(SaveDir, sprintf('MacaqueOrient_Head_%s.png', Exp{exp}));
	export_fig(Filename, '-png', '-transparent','-nocrop','-m2');
    
end



%===================== Draw spherical representation of 3D orientations
figure;
Ndivs = 360;
[X,Y,Z] = sphere(Ndivs);
sphh = surf(X,Y,Z, 'edgealpha', 0, 'facealpha', 0.2, 'facecolor', repmat(0.5,[1,3]));
lighting phong
camlight left
hold on;
for e = 1:numel(El)
    RowIndx = round(Ndivs/2 + (El(e)/90)*(Ndivs/2))+1;
    lh = plot3(X(RowIndx,:), Y(RowIndx,:), Z(RowIndx,:), '-r');
end
for a = 1:numel(Az)
    ColIndx = round(Ndivs/2 + (Az(a)/90)*(Ndivs/4))+1;
    lh = plot3(X(:,ColIndx), Y(:,ColIndx), Z(:,ColIndx), '-b');
end
n = 1;
for e = 1:numel(El)
    for a = 1:numel(Az)
        xyz(n,1) = cosd(Az(a))*cosd(El(e));
        xyz(n,2) = sind(Az(a))*cosd(El(e));
        xyz(n,3) = sind(El(e));
        ph(n) = plot3(xyz(n,1),xyz(n,2),xyz(n,3), '.g', 'markersize', 20);
        n = n+1;
    end
end
axis equal tight
view(90, 0);

EyePoint = [1,0,0];
HeadPoint = xyz(19,:);
BodyPoint = xyz(12,:);
EyeArrow = mArrow3([0 0 0], EyePoint, 'color', [1,1 0]);
HeadArrow = mArrow3([0 0 0], HeadPoint, 'color', [0,1 1]);
BodyArrow = mArrow3([0 0 0], BodyPoint, 'color', [1,0 1]);
legend([EyeArrow,HeadArrow,BodyArrow], {'Gaze direction','Head direction','Body direction'},'fontsize', 20)


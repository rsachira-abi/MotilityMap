% AnalyzeVideoSegment: 
% Analyze a video file (avi) given the start and end frame number.
%
% INPUT:
% 1. Video file name
% 2. Start and end frame
% 
% OUTPUT:
% 1. Video sequence with the mesh
% 2. Cache MAT file
%
% Author: Sachira Kuruppu


%% Initialization & Configuration

clear;
clc;

LookAhead = input('LookAhead: '); % Number of frames to look ahead for optimization.
NumElements = input('Enter number of elements: ');

Optimize = (LookAhead > 1);

fem = FiniteElementStrainCalculator(NumElements, LookAhead);

%% Process video

[vFile, vPath] = uigetfile('*.avi', 'Select the video file');

StartFrameNum = input('Enter start frame number: ');
EndFrameNum = input('Enter end frame number: ');

% Find the frames
v = VideoReader(fullfile(vPath, vFile));

FrameNum = 1;
while (FrameNum < StartFrameNum)
    frame = readFrame(v);
    FrameNum = FrameNum + 1;
end

i = 1;
while(FrameNum <= EndFrameNum)
    FrameList{i} = readFrame(v);
    i = i + 1;
    FrameNum = FrameNum + 1;
end

% Crop frames
[~, rect] = imcrop(FrameList{1});
imshow(FrameList{1})

for i = 1:length(FrameList)
    FrameList{i} = im2double(rgb2gray(imcrop(FrameList{i}, rect)));
end

%% Generate displacement fields

Width = size(FrameList{1}, 2);
Height = size(FrameList{1}, 1);

NumFrames = length(FrameList);

[MeshX, MeshY] = meshgrid(128:60:(Width - 128), 128:60:(Height - 128));
PointGrid = [MeshX(:), MeshY(:)];

fem = fem.GenerateDisplacementFields(PointGrid, @(i)FrameList{i}, NumFrames);

%% Specify boundary

figure, imshow(FrameList{1});

TopBoundary = FiniteElementStrainCalculator.SelectBoundary(2000);
RightBoundary = FiniteElementStrainCalculator.SelectBoundary(2000, TopBoundary(end,1), TopBoundary(end,2));
BottomBoundary = FiniteElementStrainCalculator.SelectBoundary(2000, RightBoundary(end,1), RightBoundary(end,2));
LeftBoundary = FiniteElementStrainCalculator.SelectBoundary(2000, BottomBoundary(end,1), BottomBoundary(end,2), TopBoundary(1,1), TopBoundary(1,2));
%%
fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

%% Fit to deformation

tic
[Px, Py] = fem.OptimizeMeshDeformation(NumFrames, [], [], Optimize, LookAhead);
toc

% Save cache
save([vFile, 'cache_2.mat'], '-v7.3');

%% Generate deformation video

disp([newline, newline]);
disp('1: Longitudinal strain');
disp('2: Transverse strain');
strain_indx = input('What strain do you want to plot? ');

vout = VideoWriter([vFile, '.output.avi']);
vout.FrameRate = 15;
open(vout);

ValRangeX = 1:0.2:(fem.m - 1);
ValRangeY = 1:0.5:(fem.m - 1);
[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

MinVal = -0.2;
MaxVal = 0.2;

NN = [];

%cl = [-2.7321, 1.2898];
cl = [-0.4, 0.4];
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);

fig = figure;

for i = 1:NumFrames
    I = FrameList{i};
    
    Px2 = Px(:,i);
    Py2 = Py(:,i);
    
    [X, Y] = fem.ExtractBoundary(Px2, Py2);
    
    [E, lambda] = fem.CalculateStrain(Px2, Py2, MaterialPoints);

    if isempty(NN)
        [X_, Y_, NN] = fem.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2);
    else
        [X_, Y_, ~] = fem.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2, NN);
    end
    
    I = insertMarker(I, [X, Y], '*', 'Color', 'yellow');
    
    imshow(I);
    hold on;
    scatter(X_, Y_, 35, lambda(:,2), 'filled');
    colormap(mycolormap);
    colorbar;
    caxis(cl);
    drawnow;
    f = getframe(fig);
    
    hold off;
    
    writeVideo(vout, f);
end

close(vout);

%% Generate strain map

keyboard; % Are you sure?

clearvars -except fem Px Py StartFrameNum EndFrameNum NumFrames
%load([vFile, 'cache.mat'], 'Px', 'Py');

[tFile, tPath] = uigetfile;
load(fullfile(tPath, tFile));

ValRangeX = 1:0.05:(fem.m - 1);
%ValRangeY = 1:0.1:(fem.m - 1);
ValRangeY = (fem.m - 2):0.1:(fem.m - 1);

tvec = timevec(StartFrameNum:EndFrameNum);

StrainMapX = zeros(NumFrames, length(ValRangeX));
StrainMapY = zeros(NumFrames, length(ValRangeX));

[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

for t = 1:NumFrames
    disp(['Processing Frame: ', num2str(t)]);
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    [E, lambda] = fem.CalculateStrain(Px2, Py2, MaterialPoints);
    
    StrainX_2D = reshape(lambda(:,1), length(ValRangeY), length(ValRangeX));
    StrainY_2D = reshape(lambda(:,2), length(ValRangeY), length(ValRangeX));
    
    StrainMapX(t,:) = nanmean(StrainX_2D, 1);
    StrainMapY(t,:) = nanmean(StrainY_2D, 1);
end

%% Display strain map

StrainMap = StrainMapY;
Ref = tvec(1);
tvec = tvec - Ref;

figure;
h = imagesc(StrainMap);
set(h, 'YData', tvec);
ylim([0, max(tvec)]);
xlabel('Length along the small intestine');
ylabel('Time (s)');

cl = [-0.4, 0.4];%caxis; %[-0.0955, 0.1744];
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);
colormap(mycolormap);
colorbar;
caxis(cl);

%% Pause

keyboard;

%% Display strain-rate map

[~, StrainRateMap] = gradient(StrainMap);

figure;
h = imagesc(StrainRateMap);
set(h, 'YData', tvec);
ylim([0, max(tvec)]);
xlabel('Length along the small intestine');
ylabel('Time (s)');

cl = caxis; %[-0.0955, 0.1744];
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);
colormap(mycolormap);
colorbar;
caxis(cl);
title('Strain rate map');

%% Display composite map

f = figure;
ax1 = axes(f);

% Strain map
h1 = imagesc(ax1, StrainMap);
set(h1, 'YData', tvec);
ylim(ax1, [0, max(tvec)]);
xlabel(ax1, 'Length along the small intestine');
ylabel(ax1, 'Time (s)');

cl = caxis(ax1);
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);
colormap(ax1, mycolormap);
caxis(ax1, cl);

colorbar(ax1, 'westoutside');

% Strain rate map
ax2 = axes(f, 'Position', ax1.Position);

AlphaMap = abs(StrainRateMap);
AlphaMap = AlphaMap ./ max(AlphaMap(:)); % normalize

h2 = imagesc(ax2, StrainRateMap, 'AlphaData', AlphaMap);
set(h2, 'YData', tvec);
ylim(ax2, [0, max(tvec)]);
% xlabel(ax2, 'Length along the small intestine');
% ylabel(ax2, 'Time (s)');

cl = caxis(ax2);
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);
colormap(ax2, mycolormap);
caxis(ax2, cl);
ax2.Color = 'none';

colorbar(ax2, 'eastoutside');

% Align again?
ax2.Position = ax1.Position;

%% Still Image

frame_num = 294;

I = FrameList{frame_num};

Px2 = Px(:,frame_num);
Py2 = Py(:,frame_num);

[X, Y] = fem.ExtractBoundary(Px2, Py2);

I = insertMarker(I, [X, Y], '*', 'Color', 'yellow');

figure, imshow(I);

%% Change strain reference

NewReference = StrainMap(451,:);

for i = 1:size(StrainMap, 1)
    StrainMap(i,:) = StrainMap(i,:) - NewReference;
end

%% Calculate velocity

h = imline;
pos = getPosition(h);

velocity = (pos(2,1) - pos(1,1)) / (pos(2,2) - pos(1,2))

%% Simplify map

LowestMean = Inf;
BestRow = 1;

for row = 1:size(StrainMap,1)
    disp(['Checking row = ', num2str(row)]);
    SubMap = repmat(StrainMap(row,:), size(StrainMap,1), 1);
    NewMap = StrainMap - SubMap;
    
    m = sum(NewMap(:) .^ 2);
    
    if (m < LowestMean)
        LowestSum = m;
        BestRow = row;
    end
end

%% Generate deformation video 2

keyboard; % You probably don't need to run this part

vout = VideoWriter([vFile, '.output.avi']);
vout.FrameRate = 15;
open(vout);

% X boundary lines
[MeshX, MeshY] = meshgrid(1:0.01:(fem.m - 1), 1:(fem.m - 1));
BoundaryMaterial = [MeshX(:), MeshY(:)];

% Y boundary lines
[MeshX, MeshY] = meshgrid(1:(fem.m - 1), 1:0.01:(fem.m - 1));
BoundaryMaterial = [BoundaryMaterial; MeshX(:), MeshY(:)];

%clearvars -except fem Px Py StartFrameNum EndFrameNum NumFrames BoundaryMaterial
%load([vFile, 'cache.mat'], 'Px', 'Py');

[tFile, tPath] = uigetfile;
load(fullfile(tPath, tFile));

NN = [];

cl = [-0.4, 0.4];
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);

fig = figure;

for t = 1:NumFrames
    disp(['Processing Frame: ', num2str(t)]);
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    [E, lambda] = fem.CalculateStrain(Px2, Py2, BoundaryMaterial);
    if isempty(NN)
        [X_, Y_, NN] = fem.CalculatePoint(BoundaryMaterial(:,1), BoundaryMaterial(:,2), Px2, Py2);
    else
        [X_, Y_, ~] = fem.CalculatePoint(BoundaryMaterial(:,1), BoundaryMaterial(:,2), Px2, Py2, NN);
    end
    
    I = cat(3, FrameList{t}, FrameList{t}, FrameList{t});
    imshow(I);
    hold on;
    scatter(X_, Y_, 35, lambda(:,2), 'filled');
    colormap(mycolormap);
    colorbar;
    caxis(cl);
    drawnow;
    f = getframe(fig);
    
    hold off;
    
    writeVideo(vout, f);
end

close(vout);

%% Boundary still figure

frame_num = 105;

Px2 = Px(:,frame_num);
Py2 = Py(:,frame_num);

% X boundary lines
[MeshX, MeshY] = meshgrid(1:0.01:(fem.m - 1), 1:(fem.m - 1));
BoundaryMaterial = [MeshX(:), MeshY(:)];

% Y boundary lines
[MeshX, MeshY] = meshgrid(1:(fem.m - 1), 1:0.01:(fem.m - 1));
BoundaryMaterial = [BoundaryMaterial; MeshX(:), MeshY(:)];

[E, lambda] = fem.CalculateStrain(Px2, Py2, BoundaryMaterial);
if isempty(NN)
    [X_, Y_, NN] = fem.CalculatePoint(BoundaryMaterial(:,1), BoundaryMaterial(:,2), Px2, Py2);
else
    [X_, Y_, ~] = fem.CalculatePoint(BoundaryMaterial(:,1), BoundaryMaterial(:,2), Px2, Py2, NN);
end



figure;

subplot(3, 1, 1);
I = cat(3, FrameList{1}, FrameList{1}, FrameList{1});
imshow(I)

subplot(3, 1, 2);
I = cat(3, FrameList{frame_num}, FrameList{frame_num}, FrameList{frame_num});
imshow(I);

hold on;
scatter(X_, Y_, 14, lambda(:,1), 'filled');
colormap(mycolormap);
caxis(cl);
hold off;

subplot(3, 1, 3);
imshow(I);

hold on;
scatter(X_, Y_, 14, lambda(:,2), 'filled');
colormap(mycolormap);
caxis(cl);
hold off;





% AnalyzeVideoSegment:
% 
% Analyze a video file (avi) given the start and end frame number.
%
% INPUT:
% 1. Video file name
% 2. Start and end frame
% 
% OUTPUT:
% 1. Video sequence with the mesh
% 2. Cache MAT file

%% Initialization

clear;
clc;

LookAhead = input('LookAhead: '); % Number of frames to look ahead for optimization.
NumElements = input('Enter number of elements: ');

Optimize = (LookAhead > 1);

%% Process video

vPath = uigetdir;
vFile = input('Enter name for the project: ', 's');

all_files = dir(fullfile(vPath, '*.raw'));
all_files = struct2cell(all_files);
all_files = cellfun(@extractFrameNumCellFunc , all_files(1,:));
all_files = sort(all_files);

StartFrameNum = input('Enter start frame number: ');
EndFrameNum = input('Enter end frame number: ');

FirstFrame = ImageRegistrationCUDA.ImportRaw(fullfile(vPath, [num2str(StartFrameNum), '.raw']));
StartFrameIndx = find(all_files == StartFrameNum);
EndFrameIndx = find(all_files == EndFrameNum);
FrameNameList = all_files(StartFrameIndx:EndFrameIndx);

% Crop frames
figure;
[FirstFrame, rect] = imcrop(FirstFrame);
imshow(FirstFrame)

%% Generate displacement fields

Width = size(FirstFrame, 2);
Height = size(FirstFrame, 1);

NumFrames = length(FrameNameList);

[MeshX, MeshY] = meshgrid(128:60:(Width - 128), 128:60:(Height - 128));
PointGrid = int32([MeshX(:), MeshY(:)]);

reg = ImageRegistrationCUDA(FirstFrame, PointGrid, 64);

DisplacementFields = cell(NumFrames, 4);

mkdir('tmp');

sz = reg.getTemplateSize();
Queue = complex(double.empty(sz(1), sz(2), 0));
QueueIndx = int16.empty(1, 0);
qLength = 0;
TemplateIndx = 0;
max_weight = 0;

for i = 1:NumFrames
    FrameName = [num2str(FrameNameList(i)), '.raw'];
    disp(['Frame = ', FrameName]);
    
    cacheFile = fullfile('tmp', [FrameName, '.mat']);
    if ~isfile(cacheFile)
        target = ImageRegistrationCUDA.ImportRaw(fullfile(vPath, FrameName));
        target = imcrop(target, rect);
        %[target, ~] = undistortImage(target, cameraParams, 'OutputView', 'full');
        target_grad = reg.grad(target);
        [padded_target_grad, ~] = reg.padArray(target_grad, []);
        
        save(cacheFile, 'padded_target_grad')
    else
        load(cacheFile, 'padded_target_grad');
    end
    
    if (qLength < LookAhead)
        qLength = qLength + 1;
        Queue(:,:,qLength) = padded_target_grad;
        QueueIndx(qLength) = i;
        disp(['--- Added, length = ', num2str(qLength)]);
    end
    
    if (qLength == LookAhead)
        disp('--- Processing');
        [Fx, Fy, error] = reg.calcDisplacementFields(Queue);
        weight = 1 ./ error;
        weight(weight > 1) = 1;
        max_weight = max(weight);
        
        TemplateIndx = TemplateIndx + 1;
        DisplacementFields{TemplateIndx, 1} = Fx;
        DisplacementFields{TemplateIndx, 2} = Fy;
        DisplacementFields{TemplateIndx, 3} = QueueIndx;
        DisplacementFields{TemplateIndx, 4} = weight;
        
        Queue(:,:,1) = [];
        QueueIndx(1) = [];
        qLength = qLength - 1;
        disp(['--- Removed, length = ', num2str(qLength)]);
    end
end

while (qLength > 0)
    disp('--- Processing');
    [Fx, Fy, error] = reg.calcDisplacementFields(Queue);
    weight = 1 ./ error;
    weight(weight > 1) = 1;
    max_weight = max(weight);
    
    TemplateIndx = TemplateIndx + 1;
    DisplacementFields{TemplateIndx, 1} = Fx;
    DisplacementFields{TemplateIndx, 2} = Fy;
    DisplacementFields{TemplateIndx, 3} = QueueIndx;
    DisplacementFields{TemplateIndx, 4} = weight;
    
    Queue(:,:,1) = [];
    QueueIndx(1) = [];
    qLength = qLength - 1;
    disp(['--- Removed, length = ', num2str(qLength)]);
end

%% Specify boundary

figure, imshow(FirstFrame);

TopBoundary = FreeFormDefStrainCalculator.SelectBoundary(2000);
RightBoundary = FreeFormDefStrainCalculator.SelectBoundary(2000, TopBoundary(end,1), TopBoundary(end,2));
BottomBoundary = FreeFormDefStrainCalculator.SelectBoundary(2000, RightBoundary(end,1), RightBoundary(end,2));
LeftBoundary = FreeFormDefStrainCalculator.SelectBoundary(2000, BottomBoundary(end,1), BottomBoundary(end,2), TopBoundary(1,1), TopBoundary(1,2));

ffd = FreeFormDefStrainCalculator(11, LookAhead, DisplacementFields, max_weight);
ffd = ffd.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

%% Fit to deformation

tic
[Px, Py] = ffd.OptimizeMeshDeformation(NumFrames, [], [], Optimize, LookAhead);
toc

% Save cache
mkdir(vFile);
save(fullfile(vFile, 'cache.mat'), '-v7.3');

%% Generate displacement fields video

vout = VideoWriter(fullfile(vFile, 'displacements.avi'));
vout.FrameRate = 15;
open(vout);

ValRangeX = 1:0.2:(ffd.m - 1);
ValRangeY = 1:0.2:(ffd.m - 1);
[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

NN = [];

cl = [-0.4, 0.4];
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);

fig = figure;

for i = 1:(NumFrames - 1)
    disp(['Writing frame: ', num2str(i)]);
    
    I = ImageRegistrationCUDA.ImportRaw(fullfile(vPath, [num2str(FrameNameList(i)), '.raw']));
    I(I > 1) = 1;
    I = imcrop(I, rect);
    
    Px2 = Px(:,i);
    Py2 = Py(:,i);
    
    if isempty(NN)
        [X_, Y_, NN] = ffd.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2);
    else
        [X_, Y_, ~] = ffd.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2, NN);
    end
    
    Fx = DisplacementFields{i,1}{1};
    Fy = DisplacementFields{i,2}{1};
    
    u = Fx(X_, Y_);
    v = Fy(X_, Y_);
    
    imshow(I);
    hold on;
    scatter(X_, Y_, 35, 'g+');
    quiver(X_, Y_, u, v);
    hold off;
    
    drawnow;
    f = getframe(fig);
    
    writeVideo(vout, f);
end

close(vout);

%% Generate deformation video

disp([newline, newline]);
disp('1: Longitudinal strain');
disp('2: Transverse strain');
strain_indx = input('What strain do you want to plot? ');

vout = VideoWriter(fullfile(vFile, 'strain.avi'));
vout.FrameRate = 15;
open(vout);

ValRangeX = 1:0.2:(ffd.m - 1);
ValRangeY = 1:0.2:(ffd.m - 1);
[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

MinVal = -50;
MaxVal = 50;

NN = [];

cl = [-0.4, 0.4];
center = -1 * cl(2) / (cl(1) - cl(2));
mycolormap = customcolormap([0, center, 1], [1, 0, 0; 1, 1, 1; 0, 0, 1]);

fig = figure;

for i = 1:NumFrames
    disp(['Writing frame: ', num2str(i)]);
    
    I = ImageRegistrationCUDA.ImportRaw(fullfile(vPath, [num2str(FrameNameList(i)), '.raw']));
    I(I > 1) = 1;
    I = imcrop(I, rect);
    
    Px2 = Px(:,i);
    Py2 = Py(:,i);
    
    [X, Y] = ffd.ExtractBoundary(Px2, Py2);
    
    [E, lambda] = ffd.CalculateStrain(Px2, Py2, MaterialPoints);

    if isempty(NN)
        [X_, Y_, NN] = ffd.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2);
    else
        [X_, Y_, ~] = ffd.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2, NN);
    end
    
    I = insertMarker(I, [X, Y], '*', 'Color', 'yellow');
    
    imshow(I);
    hold on;
    scatter(X_, Y_, 35, lambda(:,strain_indx), 'filled');
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

clearvars -except ffd Px Py StartFrameIndx EndFrameIndx NumFrames strain_indx
%load([vFile, 'cache.mat'], 'Px', 'Py');

[tFile, tPath] = uigetfile;
load(fullfile(tPath, tFile));

ValRangeX = 1:0.01:(ffd.m - 1);
ValRangeY = 1:0.1:2; %1:0.1:(ffd.m - 1);

tvec = timevec(StartFrameIndx:EndFrameIndx);

StrainMapX = zeros(NumFrames, length(ValRangeX));
StrainMapY = zeros(NumFrames, length(ValRangeX));

[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

for t = 1:NumFrames
    disp(['Processing Frame: ', num2str(t)]);
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    [E, lambda] = ffd.CalculateStrain(Px2, Py2, MaterialPoints);
    
    StrainX_2D = reshape(lambda(:,1), length(ValRangeY), length(ValRangeX));
    StrainY_2D = reshape(lambda(:,2), length(ValRangeY), length(ValRangeX));
    
    StrainMapX(t,:) = nanmean(StrainX_2D, 1);
    StrainMapY(t,:) = nanmean(StrainY_2D, 1);
end

%% Display strain map

if (strain_indx == 2)
    StrainMap = StrainMapY;
else
    StrainMap = StrainMapX;
end

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

%% Helper functions

%--------------------------------------------------------------------------
% extractFrameNumCellFunc
% Designed to be run with cellfun(). Extract the numerical frame number
% from the raw file listing names obtained with dir().
%
% Input:
%           x = Name of a single raw file. e.g.: '100.raw'
% Output:
%           y = Frame number extracted from the raw file name. e.g.: 100
%--------------------------------------------------------------------------
function y = extractFrameNumCellFunc (x)
    str_parts = split(x, '.');
    y = str2double(str_parts{1});
end


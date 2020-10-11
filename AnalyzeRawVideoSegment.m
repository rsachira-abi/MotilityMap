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
vFile = input('Enter name for the video outputs (eg: exp1.avi): ', 's');

all_files = dir(fullfile(vPath, '*.raw'));
all_files = struct2cell(all_files);

StartFrameNum = input('Enter start frame number: ');
EndFrameNum = input('Enter end frame number: ');

FirstFrame = ImageRegistration.ImportRaw(fullfile(vPath, [num2str(StartFrameNum), '.raw']));
StartFrameIndx = find(contains(all_files(1,:), num2str(StartFrameNum)));
EndFrameIndx = find(contains(all_files(1,:), num2str(EndFrameNum)));
FrameNameList = all_files(1, StartFrameIndx:EndFrameIndx);

% Crop frames
[FirstFrame, rect] = imcrop(FirstFrame);
imshow(FirstFrame)

%% Generate displacement fields

Width = size(FirstFrame, 2);
Height = size(FirstFrame, 1);

NumFrames = length(FrameNameList);

[MeshX, MeshY] = meshgrid(128:60:(Width - 128), 128:60:(Height - 128));
PointGrid = [MeshX(:), MeshY(:)];

reg = ImageRegistration(FirstFrame, PointGrid);

DisplacementFields = cell(NumFrames, 4);

mkdir('tmp');

sz = reg.getTemplateSize();
Queue = complex(double.empty(sz(1), sz(2), 0));
QueueIndx = int16.empty(1, 0);
qLength = 0;
TemplateIndx = 0;
max_weight = 0;

for i = 1:NumFrames
    FrameName = FrameNameList{i};
    disp(['Frame = ', FrameName]);
    
    cacheFile = fullfile('tmp', [FrameName, '.mat']);
    if ~isfile(cacheFile)
        target = ImageRegistration.ImportRaw(fullfile(vPath, FrameName));
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

TopBoundary = FiniteElementStrainCalculator.SelectBoundary(2000);
RightBoundary = FiniteElementStrainCalculator.SelectBoundary(2000, TopBoundary(end,1), TopBoundary(end,2));
BottomBoundary = FiniteElementStrainCalculator.SelectBoundary(2000, RightBoundary(end,1), RightBoundary(end,2));
LeftBoundary = FiniteElementStrainCalculator.SelectBoundary(2000, BottomBoundary(end,1), BottomBoundary(end,2), TopBoundary(1,1), TopBoundary(1,2));

fem = FiniteElementStrainCalculator(11, LookAhead, DisplacementFields, max_weight);
fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

%% Fit to deformation

tic
[Px, Py] = fem.OptimizeMeshDeformation(NumFrames, [], [], Optimize, LookAhead);
toc

% Save cache
mkdir('output');
save(fullfile('output', [vFile, 'cache.mat']), '-v7.3');

%% Generate deformation video

vout = VideoWriter(fullfile('output', [vFile, '.output.avi']));
vout.FrameRate = 15;
open(vout);

ValRangeX = 1:0.2:(fem.m - 1);
ValRangeY = 1:0.2:(fem.m - 1);
[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

MinVal = -50;
MaxVal = 50;

NN = [];

figure;
for i = 1:NumFrames
    I = ImageRegistration.ImportRaw(fullfile(vPath, FrameNameList{i}));
    I(I > 1) = 1;
    I = imcrop(I, rect);
    
    Px2 = Px(:,i);
    Py2 = Py(:,i);
    
    [X, Y] = fem.ExtractBoundary(Px2, Py2);
    
    E = fem.CalculateStrain(Px2, Py2, MaterialPoints);

    if isempty(NN)
        [X_, Y_, NN] = fem.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2);
    else
        [X_, Y_, ~] = fem.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Px2, Py2, NN);
    end
    
    I = insertMarker(I, [X, Y], '*', 'Color', 'yellow');
    V = E(:,1);
%     V(V > MaxVal) = MaxVal;
%     V(V < MinVal) = MinVal;
%     V = (V - MinVal) ./ (MaxVal - MinVal);
    
    imshow(I);
    hold on;
    scatter(X_, Y_, 35, V);
    colormap('jet');
    caxis([MinVal, MaxVal]);
    colorbar;
    hold off;
    
    f = getframe;
    
    writeVideo(vout, f);
end

close(vout);

%% Generate strain map

clearvars -except fem Px Py StartFrameNum EndFrameNum NumFrames
%load([vFile, 'cache.mat'], 'Px', 'Py');

[tFile, tPath] = uigetfile;
load(fullfile(tPath, tFile));

ValRangeX = 1:0.01:(fem.m - 1);
ValRangeY = (fem.m - 2):0.1:(fem.m - 1);

tvec = timevec(StartFrameNum:EndFrameNum);

StrainMap = zeros(NumFrames, length(ValRangeX));

[MatMeshX, MatMeshY] = meshgrid(ValRangeX, ValRangeY);
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

for t = 1:NumFrames
    disp(['Processing Frame: ', num2str(t)]);
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    E = fem.CalculateStrain(Px2, Py2, MaterialPoints);
    
    Stretch = E(:,1); %sqrt(2 .* E(:,1) + 1);
    ZeroCenteredStretch = Stretch;% - 1;
    
    ZeroCenteredStretch2D = reshape(ZeroCenteredStretch, length(ValRangeY), length(ValRangeX));
    
    StrainMap(t,:) = mean(ZeroCenteredStretch2D, 1);
end

%% Display strain map

limit = min(abs(min(StrainMap(:))), max(StrainMap(:)));

figure;
h = imagesc(StrainMap);
set(h, 'YData', tvec);
ylim([tvec(1), tvec(end)]);
xlabel('Length along the small intestine');
ylabel('Time (s)');

colormap('jet');
%caxis([-limit, limit]);
colorbar;

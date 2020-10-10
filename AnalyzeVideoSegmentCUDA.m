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

addpath('..');
%addpath('../SourceFunctions/');

%%
fem = FiniteElementStrainCalculator(11);

%% Process video

vPath = uigetdir;
vFile = 'RenameThis.avi';

StartFrameNum = input('Enter start frame number: ');
EndFrameNum = input('Enter end frame number: ');

FirstFrame = ImageRegistrationCUDA.ImportRaw(fullfile(vPath, [num2str(StartFrameNum), '.raw']));

% Crop frames
[FirstFrame, rect] = imcrop(FirstFrame);
imshow(FirstFrame)

%% Generate displacement fields

Width = size(FrameList{1}, 2);
Height = size(FrameList{1}, 1);

NumFrames = length(FrameList);

[MeshX, MeshY] = meshgrid(128:60:(Width - 128), 128:60:(Height - 128));
PointGrid = [MeshX(:), MeshY(:)];

reg = ImageRegistrationCUDA(template, PointGrid);

DisplacementFields = cell(NumFrames, 4);

mkdir('tmp');

sz = reg.getTemplateSize();
Queue = complex(double.empty(sz(1), sz(2), 0));
QueueIndx = int.empty(1, 0);
qLength = 0;
TemplateIndx = 0;
LookAhead = 10;

for FrameNum = StartFrameNum:EndFrameNum
    disp(['Frame = ', num2str(FrameNum)]);
    
    cacheFile = fullfile('tmp', [num2str(FrameNum), '.mat']);
    if ~isfile(cacheFile)
        target = ImageRegistrationCUDA.ImportRaw(fullfile(vpath, [num2str(FrameNum), '.raw']));
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
        QueueIndx(qLength) = FrameNum;
        disp(['--- Added, length = ', num2str(qLength)]);
    end
    
    if (qLength == LookAhead)
        disp('--- Processing');
        [Fx, Fy, error] = reg.calcDisplacementFields(Queue);
        
        TemplateIndx = TemplateIndx + 1;
        DisplacementFields{TemplateIndx, 1} = Fx;
        DisplacementFields{TemplateIndx, 2} = Fy;
        DisplacementFields{TemplateIndx, 3} = QueueIndx;
        DisplacementFields{TemplateIndx, 4} = error;
        
        Queue(:,:,1) = [];
        QueueIndx(1) = [];
        qLength = qLength - 1;
        disp(['--- Removed, length = ', num2str(qLength)]);
    end
end

while (qLength > 0)
    disp('--- Processing');
    [Fx, Fy, error] = reg.calcDisplacementFields(Queue);
    
    TemplateIndx = TemplateIndx + 1;
    DisplacementFields{TemplateIndx, 1} = Fx;
    DisplacementFields{TemplateIndx, 2} = Fy;
    DisplacementFields{TemplateIndx, 3} = QueueIndx;
    DisplacementFields{TemplateIndx, 4} = error;
    
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

fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

%% Fit to deformation

tic
[Px, Py] = fem.OptimizeMeshDeformation(NumFrames, [], [], false, 10);
toc

% Save cache
save([vFile, 'cache.mat'], '-v7.3');

%% Generate deformation video

vout = VideoWriter([vFile, '.output.avi']);
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
    I = FrameList{i};
    I(I > 1) = 1;
    
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

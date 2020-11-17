
% Unit tests for StrainCalculator class

clear;
clc;

addpath('..');

%% Zero error strain test

fem = FreeFormDefStrainCalculator(11);

[DispFieldX, DispFieldY] = meshgrid(0:10:400, 0:10:400);

[MeshX, MeshY] = meshgrid(100:10:300, 100:10:300);

KnownStretch =...
    [
    1, 1;
    1.01, 1.01;
    1.02, 1.02;
    1.03, 1.03;
    1.04, 1.04;
    1.05, 1.05;
    1.06, 1.06;
    1.07, 1.07;
    1.08, 1.08;
    1.09, 1.09
    ];

for i = 1:size(KnownStretch, 1)
    PointGrid{i} = [DispFieldX(:), DispFieldY(:)];
    
    for j = (i + 1):size(KnownStretch, 1)
        relative_stretch = (KnownStretch(j) - KnownStretch(i)) / KnownStretch(i);
        
        Disp{i,j} = [PointGrid{i}(:,1) .* relative_stretch, PointGrid{i}(:,2) .* relative_stretch];
    end
end

NumFrames = length(PointGrid);

fem = fem.GenerateDummyDisplacementFields3(PointGrid, NumFrames, 0, Disp);

TopBoundary = [MeshX(1,:)', MeshY(1,:)'];
RightBoundary = [MeshX(:,end), MeshY(:,end)];
BottomBoundary = [flipud(MeshX(end,:)'), flipud(MeshY(end,:)')];
LeftBoundary = [flipud(MeshX(:,1)), flipud(MeshY(:,1))];

fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

tic
[Px, Py] = fem.OptimizeMeshDeformation(NumFrames, [], [], false, 10);
toc

total_rms_error = 0;

[MatMeshX, MatMeshY] = meshgrid(1:0.1:(fem.m - 1), 1:0.1:(fem.m - 1));

MaterialPoints = [MatMeshX(:), MatMeshY(:)];
%MaterialPoints = [0, 0];

for t = 1:NumFrames
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    [E, lambda] = fem.CalculateStrain(Px2, Py2, MaterialPoints);
    
    StretchX = lambda(:,1) + 1; %sqrt(2 .* E(:,1) + 1);
    StretchY = lambda(:,2) + 1; %sqrt(2 .* E(:,4) + 1);
    StretchX = round(StretchX, 4);
    StretchY = round(StretchY, 4);
    
    rms_error = sqrt(mean((StretchX - KnownStretch(t,2)).^2 + (StretchY - KnownStretch(t,1)).^2));
    disp(['Test ', num2str(t), ': StretchX = ', num2str(mean(StretchX)), ', StretchY = ', num2str(mean(StretchY)), ', error = ', num2str(rms_error)])
    
    total_rms_error(t) = rms_error;
end

mean_rms = mean(total_rms_error);
disp(['Mean RMS error = ', num2str(mean_rms)]);

try
    assert(mean_rms == 0);
catch
    warning('Test Failed!');
    keyboard;
end

%% Strain test without optimization

fem = FreeFormDefStrainCalculator(11);

[DispFieldX, DispFieldY] = meshgrid(0:10:400, 0:10:400);

[MeshX, MeshY] = meshgrid(100:10:300, 100:10:300);

KnownStretch =...
    [
    1, 1;
    1.01, 1.01;
    1.02, 1.02;
    1.03, 1.03;
    1.04, 1.04;
    1.05, 1.05;
    1.06, 1.06;
    1.07, 1.07;
    1.08, 1.08;
    1.09, 1.09
    ];

for i = 1:size(KnownStretch, 1)
    PointGrid{i} = [DispFieldX(:), DispFieldY(:)];
    
    for j = (i + 1):size(KnownStretch, 1)
        relative_stretch = (KnownStretch(j) - KnownStretch(i)) / KnownStretch(i);
        
        Disp{i,j} = [PointGrid{i}(:,1) .* relative_stretch, PointGrid{i}(:,2) .* relative_stretch];
    end
end

NumFrames = length(PointGrid);

fem = fem.GenerateDummyDisplacementFields3(PointGrid, NumFrames, 0.1, Disp);

TopBoundary = [MeshX(1,:)', MeshY(1,:)'];
RightBoundary = [MeshX(:,end), MeshY(:,end)];
BottomBoundary = [flipud(MeshX(end,:)'), flipud(MeshY(end,:)')];
LeftBoundary = [flipud(MeshX(:,1)), flipud(MeshY(:,1))];

fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

tic
[Px, Py] = fem.OptimizeMeshDeformation(NumFrames, [], [], false, 10);
toc

total_rms_error = 0;

[MatMeshX, MatMeshY] = meshgrid(1:0.1:(fem.m - 1), 1:0.1:(fem.m - 1));

MaterialPoints = [MatMeshX(:), MatMeshY(:)];
%MaterialPoints = [0, 0];

for t = 1:NumFrames
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    [E, lambda] = fem.CalculateStrain(Px2, Py2, MaterialPoints);
    
    StretchX = lambda(:,1) + 1; %sqrt(2 .* E(:,1) + 1);
    StretchY = lambda(:,2) + 1; %sqrt(2 .* E(:,4) + 1);
    StretchX = round(StretchX, 4);
    StretchY = round(StretchY, 4);
    
    rms_error = sqrt(mean((StretchX - KnownStretch(t,2)).^2 + (StretchY - KnownStretch(t,1)).^2));
    disp(['Test ', num2str(t), ': StretchX = ', num2str(mean(StretchX)), ', StretchY = ', num2str(mean(StretchY)), ', error = ', num2str(rms_error)])
    
    total_rms_error(t) = rms_error;
end

mean_rms = mean(total_rms_error);
disp(['Mean RMS error = ', num2str(mean_rms)]);

try
    assert(mean_rms < 0.1);
catch
    warning('Test Failed!');
    keyboard;
end

%% Strain test with optimization

fem = FreeFormDefStrainCalculator(11);

[DispFieldX, DispFieldY] = meshgrid(0:10:400, 0:10:400);

[MeshX, MeshY] = meshgrid(100:10:300, 100:10:300);

KnownStretch =...
    [
    1, 1;
    1.01, 1.01;
    1.02, 1.02;
    1.03, 1.03;
    1.04, 1.04;
    1.05, 1.05;
    1.06, 1.06;
    1.07, 1.07;
    1.08, 1.08;
    1.09, 1.09
    ];

for i = 1:size(KnownStretch, 1)
    PointGrid{i} = [DispFieldX(:), DispFieldY(:)];
    
    for j = (i + 1):size(KnownStretch, 1)
        relative_stretch = (KnownStretch(j) - KnownStretch(i)) / KnownStretch(i);
        
        Disp{i,j} = [PointGrid{i}(:,1) .* relative_stretch, PointGrid{i}(:,2) .* relative_stretch];
    end
end

NumFrames = length(PointGrid);

fem = fem.GenerateDummyDisplacementFields3(PointGrid, NumFrames, 0.1, Disp);

TopBoundary = [MeshX(1,:)', MeshY(1,:)'];
RightBoundary = [MeshX(:,end), MeshY(:,end)];
BottomBoundary = [flipud(MeshX(end,:)'), flipud(MeshY(end,:)')];
LeftBoundary = [flipud(MeshX(:,1)), flipud(MeshY(:,1))];

fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

tic
[Px, Py] = fem.OptimizeMeshDeformation(NumFrames, [], [], true, 10);
toc

total_rms_error = 0;

[MatMeshX, MatMeshY] = meshgrid(1:0.1:(fem.m - 1), 1:0.1:(fem.m - 1));

MaterialPoints = [MatMeshX(:), MatMeshY(:)];
%MaterialPoints = [0, 0];

for t = 1:NumFrames
    Px2 = Px(:,t);
    Py2 = Py(:,t);
    
    [E, lambda] = fem.CalculateStrain(Px2, Py2, MaterialPoints);
    
    StretchX = lambda(:,1) + 1; %sqrt(2 .* E(:,1) + 1);
    StretchY = lambda(:,2) + 1; %sqrt(2 .* E(:,4) + 1);
    StretchX = round(StretchX, 4);
    StretchY = round(StretchY, 4);
    
    rms_error = sqrt(mean((StretchX - KnownStretch(t,2)).^2 + (StretchY - KnownStretch(t,1)).^2));
    disp(['Test ', num2str(t), ': StretchX = ', num2str(mean(StretchX)), ', StretchY = ', num2str(mean(StretchY)), ', error = ', num2str(rms_error)])
    
    total_rms_error(t) = rms_error;
end

mean_rms = mean(total_rms_error);
disp(['Mean RMS error = ', num2str(mean_rms)]);

try
    assert(mean_rms < 0.1);
catch
    warning('Test Failed!');
    keyboard;
end

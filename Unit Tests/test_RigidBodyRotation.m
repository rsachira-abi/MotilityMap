% This script test strain error during rigid body rotation, with our
% method, and spatial differentiated strain-rate.

clear;
clc;

addpath('..');

%% Generate displacement fields

[DispFieldX, DispFieldY] = meshgrid(-400:10:400, -400:10:400);
corner_indices = sub2ind(size(DispFieldX), [1; 1; size(DispFieldX,1); size(DispFieldX,1); 1], [1; size(DispFieldX,2); size(DispFieldX,2); 1; 1]);

PointGrid = [DispFieldX(:), DispFieldY(:), ones(size(DispFieldY(:)))];

rotz_ = @(rad) [cos(rad), -sin(rad), 0; sin(rad), cos(rad), 0; 0, 0, 1];
rotz = @(theta) rotz_(deg2rad(theta));

newPointGrid = (rotz(1) * PointGrid')';
newPointGrid = round(newPointGrid, 1);
%newPointGrid = (rotz(-40) * newPointGrid')';
%newPointGrid = PointGrid;%(rotz(deg2rad(5)) * PointGrid')';
%newPointGrid(:,1) = newPointGrid(:,1) + 5;

DispX = newPointGrid(:,1) - PointGrid(:,1);
DispY = newPointGrid(:,2) - PointGrid(:,2);

figure;

plot(PointGrid(corner_indices,1), PointGrid(corner_indices,2), 'b-'); hold on;
scatter(PointGrid(:,1), PointGrid(:,2), [], 'b');

plot(newPointGrid(corner_indices,1), newPointGrid(corner_indices,2), 'r-');
scatter(newPointGrid(:,1), newPointGrid(:,2), [], 'r+');

quiver(PointGrid(:,1), PointGrid(:,2), DispX, DispY, 'k');

%% Test with our method

[MeshX, MeshY] = meshgrid(-300:10:300, -300:10:300);

Disp{1,2} = [DispX, DispY];

% Generate displacement fields for the rotation.
Fx = scatteredInterpolant(PointGrid(:,1), PointGrid(:,2), DispX, 'natural');
Fy = scatteredInterpolant(PointGrid(:,1), PointGrid(:,2), DispY, 'natural');

DisplacementFields = cell(2, 4);
DisplacementFields{1, 1} = {Fx};
DisplacementFields{1, 2} = {Fy};
DisplacementFields{1, 3} = [2];
DisplacementFields{1, 4} = [1];

fem = FreeFormDefStrainCalculator(4, 1, DisplacementFields, 1);

TopBoundary = [MeshX(1,:)', MeshY(1,:)'];
RightBoundary = [MeshX(:,end), MeshY(:,end)];
BottomBoundary = [flipud(MeshX(end,:)'), flipud(MeshY(end,:)')];
LeftBoundary = [flipud(MeshX(:,1)), flipud(MeshY(:,1))];

fem = fem.FitInitialMesh({TopBoundary, RightBoundary, BottomBoundary, LeftBoundary});

[MatMeshX, MatMeshY] = meshgrid(1:0.1:(fem.m - 1), 1:0.1:(fem.m - 1));
MaterialPoints = [MatMeshX(:), MatMeshY(:)];

[Px, Py] = fem.OptimizeMeshDeformation(2, [], [], false, 2);
Px2 = Px(:,2);
Py2 = Py(:,2);
    
[E, lambda] = fem.CalculateStrain(Px2, Py2, MaterialPoints);

[bX, bY] = fem.ExtractBoundary(Px2, Py2);
scatter(bX, bY, 'y');

mean(lambda, 1)
std(lambda, 1)

%% With spatial differentiation










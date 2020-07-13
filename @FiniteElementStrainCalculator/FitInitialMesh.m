%-------------------------------------------------------------------
% FitInitialMesh:
% Fit the initial mesh to the given boundary.
%
% INPUT:
% image : Image for which the mesh should be defined. Used to show
%         this to the user to allow him to mark the boundary of
%         the mesh.
%
% NumNodesPerRow : Number of nodes in a row of the mesh.
% NumNodesPerCol : Number of nodes in a column of the mesh.
%
% OUTPUT:
% No outputs. But sets Elements, A, Xg_initial, Yg_initial
%-------------------------------------------------------------------
function this = FitInitialMesh (this, Boundary)
    %% Identify boundaries
    
    NumNodesPerRow = 5 * this.n;
    NumNodesPerCol = 5 * this.n;
    
    KnownPoints = [];
    KnownPointMaterial = [];
    
    [NodesC, NodesR] = meshgrid(0:NumNodesPerRow, 0:NumNodesPerCol);
    PredictedMaterialX = normalize(NodesC, 2, 'range');
    PredictedMaterialY = normalize(NodesR, 1, 'range');

    % Select top bounndary
    [TopBoundaryNodes, TopBoundaryPoints, TopBoundaryMaterial] = this.FindBoundaryNodes(NumNodesPerRow, Boundary{1}, PredictedMaterialX(1,:)); %this.SelectBoundary(NumNodesPerRow);
    KnownPointMaterial = [KnownPointMaterial; TopBoundaryMaterial, zeros(size(TopBoundaryPoints,1),1)];
    KnownPoints = [KnownPoints; TopBoundaryPoints];

    % Select right boundary
    [RightBoundaryNodes, RightBoundaryPoints, RightBoundaryMaterial] = this.FindBoundaryNodes(NumNodesPerCol, [TopBoundaryNodes(end,:); Boundary{2}], PredictedMaterialY(:,end));
    KnownPointMaterial = [KnownPointMaterial; ones(size(RightBoundaryPoints,1),1), (RightBoundaryMaterial)];
    KnownPoints = [KnownPoints; RightBoundaryPoints];

    % Select bottom boundary
    [BottomBoundaryNodes, BottomBoundaryPoints, BottomBoundaryMaterial] = this.FindBoundaryNodes(NumNodesPerRow, flipud([RightBoundaryNodes(end,:); Boundary{3}]), PredictedMaterialX(end,:));
    KnownPointMaterial = [KnownPointMaterial; (BottomBoundaryMaterial), ones(size(BottomBoundaryPoints,1),1)];
    KnownPoints = [KnownPoints; BottomBoundaryPoints];

    % Select left boundary
    [LeftBoundaryNodes, LeftBoundaryPoints, LeftBoundaryMaterial] = this.FindBoundaryNodes(NumNodesPerCol, flipud([BottomBoundaryNodes(1,:); Boundary{4}]), PredictedMaterialY(:,1));
    KnownPointMaterial = [KnownPointMaterial; zeros(size(LeftBoundaryPoints,1),1), (LeftBoundaryMaterial)];
    KnownPoints = [KnownPoints; LeftBoundaryPoints];
    
    KnownPointMaterial = KnownPointMaterial .* this.m;
    
    %% Guess internal points
    Fx = scatteredInterpolant(KnownPointMaterial(:,1), KnownPointMaterial(:,2), KnownPoints(:,1));
    Fy = scatteredInterpolant(KnownPointMaterial(:,1), KnownPointMaterial(:,2), KnownPoints(:,2));
    
    KnownPointMaterial = [KnownPointMaterial; PredictedMaterialX(:) * this.m, PredictedMaterialY(:) * this.m];
    KnownPoints = [KnownPoints; Fx([PredictedMaterialX(:) * this.m, PredictedMaterialY(:) * this.m]), Fy([PredictedMaterialX(:) * this.m, PredictedMaterialY(:) * this.m])];

    %%  Fit
    
%     [meshX, meshY] = meshgrid(0:10, 0:10);
%     MatMeshX = normalize(meshX, 2, 'range') .* this.m;
%     MatMeshY = normalize(meshY, 1, 'range') .* this.m;
%     
%     KnownPoints = [meshX(:), meshY(:)];
%     KnownPointMaterial = [MatMeshX(:), MatMeshY(:)];
    
%     KnownPointMaterial = [PredictedMaterialX(:), PredictedMaterialY(:)] .* this.m;
%     KnownPoints = [Fx(KnownPointMaterial(:,1), KnownPointMaterial(:,2)), Fy(KnownPointMaterial(:,1), KnownPointMaterial(:,2))];
    
    [Px, Py] = this.LeastSquaresMeshFit(KnownPoints, KnownPointMaterial);
    
    this.Px_initial = Px;
    this.Py_initial = Py;
end







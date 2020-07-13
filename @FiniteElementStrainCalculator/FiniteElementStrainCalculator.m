%-------------------------------------------------------------------
% FiniteElementStrainCalculator:
%
% Calculates the strain at each timestep given a set of node
% locations at those timesteps. Node locations at the 0th timestep
% is considered as the undeformed state.
%
% This class fits a bicubic Hermite mesh on to the node points and
% use that to calculate the Green strain. Method given in paper:
%
% Malcolm, D., Nielsen, P., Hunter, P. et al. Biomechan Model
% Mechanobiol (2002) 1: 197. https://doi.org/10.1007/s10237-002-0018-8
%-------------------------------------------------------------------
classdef FiniteElementStrainCalculator
    properties (Constant)
        % Configuration of the displacement fields.
        PIECE_SIZE = 10;
        ERROR_CUTOFF = 128 * 128 * 0.001; % 0.1% error
    end

    properties (Access = private)
        LOOK_AHEAD = 1;
        NumElementParameters = 16;

        DisplacementFields = {};
        MaxWeight = 0; % To normalize the weight
        
        MaterialPoints = [];    % Material point locations of the tracking points.
        TrackingPoints = [];    % Tracking points at each timestep [X, Y] x timestep
        AssignedElement = [];
        
        n = 4;
        k = 2;  % quadratic
    end
    
    properties (GetAccess = public, SetAccess = private)
        % Control points
        Px_initial = [];
        Py_initial = [];
        
        m = 0;  % Number of knot segments
    end
    
    methods (Access = public)
        %-------------------------------------------------------------------
        % Constructor:
        % Inititlizes MaterialPoints, NodePoints, Elements and A.
        %
        % INPUT:
        % NodePoints : Node locations at different timesteps
        %              [X, Y] x timestep.
        %-------------------------------------------------------------------
        function this = FiniteElementStrainCalculator (n, look_ahead)
            if (nargin >= 1)
                this.n = n;
            end
            if (nargin == 2)
                this.LOOK_AHEAD = look_ahead;
            end
            this.m = this.n - this.k + 1;
        end
        
        function this = SetN (this, n)
            this.n = n;
            this.m = this.n - this.k + 1;
        end
        
        %-------------------------------------------------------------------
        % GenerateDisplacementFields:
        % Generate the displacement fields for a frame sequence.
        %
        % INPUTS:
        % PointGrid = Grid of points used to generate the displacement
        %             field from.
        % LoadFrameCallback = Callback function that accepts a frame number
        %                     and returns that frame. e.i. 
        %                     frame <- f(frame_number)
        %
        % OUTPUT:
        % No return. Sets the 'DisplacementFields' private member variable.
        %-------------------------------------------------------------------
        this = GenerateDisplacementFields (this, PointGrid, LoadFrameCallback, NumFrames);
        this = GenerateDummyDisplacementFields (this, PointGrid, NumFrames, noise_std, signal_level, LookAhead);
        this = GenerateDummyDisplacementFields2 (this, PointGrid, NumFrames, noise_std, signal_level, LookAhead);
        this = GenerateDummyDisplacementFields3 (this, PointGrid, NumFrames, noise_std, signal_level, LookAhead);
        
        %------------------------------------------------------------------
        % LoadDisplacementFields:
        %
        % Loads previously generated set of displacement fields from a
        % file.
        %
        % INPUT:
        % filepath = Path to dispalcement field file.
        %
        % OUTPUT:
        % No return. Sets the 'DispalcementField' private member variable.
        %------------------------------------------------------------------
        function this = LoadDisplacementFields (this, filepath)
                tmp = importdata(filepath);
                this.DisplacementFields = tmp.displacement_fields;
                this.MaxWeight = tmp.max_weight;
        end
        
        %------------------------------------------------------------------
        % SaveDisplacementFields:
        %
        % Save the displacement fields in the 'DisplacementFields' member
        % variable to a file.
        %
        % INPUT:
        % filepath = Path to the save location, including the file name.
        %------------------------------------------------------------------
        function SaveDisplacementFields (this, filepath)
            displacement_fields = this.DisplacementFields;
            max_weight = this.MaxWeight;
            
            save(filepath, 'displacement_fields', 'max_weight');
        end
    
        %-------------------------------------------------------------------
        % FitInitialMesh:
        % Fit the initial mesh to the given boundary.
        % 
        % INPUT:
        % Boundary : Boundary to which the mesh should be fitted.
        %
        % NumNodesPerRow : Number of nodes in a row of the mesh.
        % NumNodesPerCol : Number of nodes in a column of the mesh.
        %
        % OUTPUT:
        % No outputs. But sets Elements, A, Xg_initial, Yg_initial
        %-------------------------------------------------------------------
        this = FitInitialMesh (this, Boundary);
        
        %-------------------------------------------------------------------
        % OptimizeDeformation:
        % Iteratively fits the mesh to frames from MinFrame to MaxFrames.
        % This is done by searching for the nodal parameters that best
        % explains the displacement fields.
        %
        % INPUT:
        % XgGuess : Initial guess of X nodal parameters 
        %           (num parameters x num frames)
        % YgGuess : Initial guess of Y nodal parameters
        % MinFrame: Starting frame
        % MaxFrame: Last frame
        %
        % OUTPUT:
        % Xg, Yg : Fitted nodal parameters (num parameters x num frames)
        %-------------------------------------------------------------------
        [Px, Py] = OptimizeMeshDeformation (this, NumFrames, Px0, Py0, DoOptimize, Iterations);
        [Px, Py] = OptimizeMeshDeformation2 (this, NumFrames, Px0, Py0, DoOptimize, Iterations);
        
        %-------------------------------------------------------------------
        % CalculateStrain:
        % Calculate the strain tensor for the mesh defined by the global
        % mesh parameters Xg, Yg.
        %
        % INPUT:
        % Xg, Yg: Global mesh parameters.
        %
        % OUTPUT:
        % X, Y         : Global (euler) coordinates of the strain tensor.
        % StrainTensor : Strain tensor for the global (euler) coordinate
        %                location in X, Y.
        %-------------------------------------------------------------------
        [E, lambda] = CalculateStrain (this, Px, Py, MaterialPoints);
        [E, lambda] = CalculateStrainRelative (this, Px, Py, Px0, Py0, MaterialPoints);
        %[X, Y, StrainTensor] = CalculateStrainV2 (this, Xg, Yg);
        
        %-------------------------------------------------------------------
        % ExtractBoundary:
        % Extract the boundary of the mesh given the global X and Y
        % parameters.
        %
        % INPUT:
        % Xg    : Global X parameters of the nodes
        % Yg    : Global Y parameters of the nodes
        %
        % OUTPUT:
        % boundary  : Boundary of the mesh {element number}[X, Y]
        %-------------------------------------------------------------------
        function [X, Y] = ExtractBoundary (this, Px, Py)
            NN = [];
            r = 0;
            
            BoundaryMaterial = [];
            
            % X boundary lines
            [MeshX, MeshY] = meshgrid(1:0.01:(this.m - 1), 1:(this.m - 1));
            BoundaryMaterial = [BoundaryMaterial; MeshX(:), MeshY(:)];
            
            % Y boundary lines
            [MeshX, MeshY] = meshgrid(1:(this.m - 1), 1:0.01:(this.m - 1));
            BoundaryMaterial = [BoundaryMaterial; MeshX(:), MeshY(:)];
            
            c = 1;
            for i = 0:this.n
                for j = 0:this.n
                    NN(:,c) = this.N_(i, this.k, BoundaryMaterial(:,1)) .* this.N_(j, this.k, BoundaryMaterial(:,2));
                    c = c + 1;
                end
            end
            
            X = NN * Px;
            Y = NN * Py;
        end
        
        %-------------------------------------------------------------------
        % CalculatePoint:
        % Calculate the global (euler) coordinates of the material point
        % (x,y), based on its element (assigned element) described by the
        % element parameters Ux, Uy.
        %
        % INPUT:
        % (x, y) : Material coordinates of the point.
        % Ux, Uy : X and Y element properties of the assigned element.
        %
        % OUTPUT:
        % Point : Euler coordinates of the material point.
        %-------------------------------------------------------------------
        function [X, Y, NN] = CalculatePoint (this, E1, E2, Px, Py, NN)
            if (nargin < 6)
                c = 1;
                for i = 0:this.n
                    for j = 0:this.n
                        NN(:,c) = this.N_(i, this.k, E1) .* this.N_(j, this.k, E2);
                        c = c + 1;
                    end
                end
            end
            
            X = NN * Px;
            Y = NN * Py;
        end
    end
    
    methods (Access = private)
        
        [Px, Py] = LeastSquaresMeshFit (this, KnownPoints, MaterialPoints);

        %-------------------------------------------------------------------
        % FitToDeformation:
        % Iteratively fits the mesh to frames from MinFrame to MaxFrames.
        % This is done by searching for the nodal parameters that best
        % explains the displacement fields.
        %
        % INPUT:
        % Xg0,Yg0   : Current guess of guess of X and Ynodal parameters
        %           that describe the deformation.
        %           (num parameters x num frames)
        % MinFrame: Starting frame
        % MaxFrame: Last frame
        %
        % OUTPUT:
        % ssd : Fitness of the current parameters (sum of squared difference).
        %-------------------------------------------------------------------
        [ssd, g] = NonlinearMinimization (this, PxPy, GivenPx, GivenPy, GivenRange, CurrentRange, ConvertionMatrix);

        g = OptimizationDerivatives (this, Ux0, Uy0, Ux1, Uy1, MaterialPoints, Fx, Fy);
        
        %-------------------------------------------------------------------
        % FindBoundaryNodes:
        % Find the boundary nodes to determine elements, based on the
        % boundary points given.
        %
        % INPUT:
        % NumNodes      : Number of nodes to use for fitting.
        % Points        : User selected boundary points.
        % ZeroX, ZeroY  : Hardcoded coordinates of zero length. 
        %
        % OUTPUT:
        % BoundaryNodes : Nodes of the boundary.
        %-------------------------------------------------------------------
        [BoundaryNodes, BoundaryPoints, BoundaryPointsMaterial] = FindBoundaryNodes (~, NumNodes, Points, NodeMaterial);
        
        %-------------------------------------------------------------------
        % ElementBoundary:
        % Output the set of points representing the boundary of the element
        % described by Ux, and Uy.
        %
        % INPUT:
        % Ux, Uy   : Element parameters that describe the element.
        %
        % OUTPUT:
        % (Xs, Ys) : Set of coordinates (X, Y) of the boundary of the 
        %            element. 
        %-------------------------------------------------------------------
        [Xs, Ys]  = ElementBoundary (this, Ux, Uy);
        
        %-------------------------------------------------------------------
        % T:
        % Returns the knot value corresponding to i
        %-------------------------------------------------------------------
        function t = T (this, i)
            if (i < this.k)
                t = 0;
            elseif (i >= this.k && i <= this.n)
                t = i - this.k;
            elseif (i > this.n)
                t = this.n - this.k + 1;
            end
        end
        
        %-------------------------------------------------------------------
        % N_:
        % Basis function
        %-------------------------------------------------------------------
        function V = N_ (this, i, k, t)
            V = zeros(size(t));
            
            if (k == 0)
                V(t >= this.T(i) & t < this.T(i + 1)) = 1;
                return;
            end
            
            N1 = this.N_(i, k - 1, t);
            N2 = this.N_(i + 1, k - 1, t);
            
            N1_indx = find(N1);
            N2_indx = find(N2);
            
            if ~isempty(N1_indx)
                V(N1_indx) = (t(N1_indx) - this.T(i)) ./ (this.T(i + k) - this.T(i)) .* N1(N1_indx);
            end
            
            if ~isempty(N2_indx)
                V(N2_indx) = V(N2_indx) + ((this.T(i + k + 1) - t(N2_indx)) ./ (this.T(i + k + 1) - this.T(i + 1)) .* N2(N2_indx));
            end
        end
        
        function [dx_dE1, dx_dE2, dy_dE1, dy_dE2] = Derivatives (this, Px, Py, MaterialPoints)
            Px_2d = reshape(Px, this.n + 1, this.n + 1)';
            Py_2d = reshape(Py, this.n + 1, this.n + 1)';
            
            dx_dE1 = zeros(size(MaterialPoints, 1), 1);
            dx_dE2 = zeros(size(MaterialPoints, 1), 1);
            dy_dE1 = zeros(size(MaterialPoints, 1), 1);
            dy_dE2 = zeros(size(MaterialPoints, 1), 1);
            
            for i = 0:(this.n - 1)
                for j = 0:this.n
                    Qx = (Px_2d(i + 2, j + 1) - Px_2d(i + 1, j + 1)) * this.k / (this.T(i + this.k + 1) - this.T(i + 1));
                    Qy = (Py_2d(i + 2, j + 1) - Py_2d(i + 1, j + 1)) * this.k / (this.T(i + this.k + 1) - this.T(i + 1));
                    
                    dx_dE1 = dx_dE1 + this.N_(i + 1, this.k - 1, MaterialPoints(:,1)) .* this.N_(j, this.k, MaterialPoints(:,2)) .* Qx;
                    dy_dE1 = dy_dE1 + this.N_(i + 1, this.k - 1, MaterialPoints(:,1)) .* this.N_(j, this.k, MaterialPoints(:,2)) .* Qy;
                end
            end
            
            for i = 0:this.n
                for j = 0:(this.n - 1)
                    Qx = (Px_2d(i + 1, j + 2) - Px_2d(i + 1, j + 1)) * this.k / (this.T(j + this.k + 1) - this.T(j + 1));
                    Qy = (Py_2d(i + 1, j + 2) - Py_2d(i + 1, j + 1)) * this.k / (this.T(j + this.k + 1) - this.T(j + 1));
                    
                    dx_dE2 = dx_dE2 + this.N_(i, this.k, MaterialPoints(:,1)) .* this.N_(j + 1, this.k - 1, MaterialPoints(:,2)) .* Qx;
                    dy_dE2 = dy_dE2 + this.N_(i, this.k, MaterialPoints(:,1)) .* this.N_(j + 1, this.k - 1, MaterialPoints(:,2)) .* Qy;
                end
            end
            
            dx_dE1 = round(dx_dE1, 4);
            dx_dE2 = round(dx_dE2, 4);
            dy_dE1 = round(dy_dE1, 4);
            dy_dE2 = round(dy_dE2, 4);
        end
        
        [E, F] = GenerateConversionMatrices (this, MaterialPoints);
        
    end
    
    methods (Static)     
        %-------------------------------------------------------------------
        % SelectBoundary:
        % Let the user select the boundary of the mesh. Used for initial
        % fitting.
        %
        % INPUT:
        % NumNodes      : Number of nodes to use for fitting.
        % ZeroX, ZeroY  : Hardcoded coordinates of zero length. 
        %
        % OUTPUT:
        % BoundaryNodes : Nodes of the boundary.
        %-------------------------------------------------------------------
        BoundaryNodes = SelectBoundary (NumNodes, ZeroX, ZeroY, LastX, LastY);
    end
end
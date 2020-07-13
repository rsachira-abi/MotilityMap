%-------------------------------------------------------------------
% OptimizeDeformation:
% Iteratively fits the mesh to frames from MinFrame to MaxFrames.
% This is done by searching for the nodal parameters that best
% explains the displacement fields.
%
% OUTPUT:
% Xg, Yg : Fitted nodal parameters (num parameters x num frames)
%-------------------------------------------------------------------
function [Px, Py] = OptimizeMeshDeformation (this, NumFrames, Px0, Py0, DoOptimize, Iterations)
    if isempty(Px0) || isempty(Py0)
        % Initial guesses of the nodal parameters.
        % Number of nodes x Number of frames.
        Px(:,1) = this.Px_initial;
        Py(:,1) = this.Py_initial;
    else
        Px(:,1) = Px0;
        Py(:,1) = Py0;
    end
    
    PieceSize = this.PIECE_SIZE;
    
    StepSize = floor(PieceSize / 2);
    
    [MeshX, MeshY] = meshgrid(0:(this.n * 5), 0:(this.n * 5));
    MaterialX = normalize(MeshX, 2, 'range') .* this.m;
    MaterialY = normalize(MeshY, 1, 'range') .* this.m;
    MaterialPoints = [MaterialX(:),  MaterialY(:)];

    [E, F] = this.GenerateConversionMatrices(MaterialPoints);
    NN = F';
    
    opts = optimoptions('fminunc', 'Display', 'iter-detailed', 'FiniteDifferenceType', 'forward', 'StepTolerance', 1e-10, 'FunctionTolerance', 1e-4,...
            'OptimalityTolerance', 1e-4, 'MaxFunctionEvaluations', 1e10, 'MaxIterations', Iterations, 'UseParallel', true);
    
    for CurrentFrame = 1:StepSize:NumFrames
        disp(['-- Processing frames from ', num2str(CurrentFrame), ' to ', num2str(CurrentFrame + StepSize - 1)]);

        EulerPointsX = NN * Px(:,CurrentFrame);
        EulerPointsY = NN * Py(:,CurrentFrame);
        
        PointsArray = [EulerPointsX, EulerPointsY];
        
        InitializationList = CurrentFrame + 1:min([NumFrames, (CurrentFrame + StepSize)]);
        
        if ~isempty(InitializationList)
            for i = InitializationList
                Fx = this.DisplacementFields{i - 1,1}{1};
                Fy = this.DisplacementFields{i - 1,2}{1};
                
                PointsArray(:,1) = PointsArray(:,1) + Fx(PointsArray(:,1), PointsArray(:,2));
                PointsArray(:,2) = PointsArray(:,2) + Fy(PointsArray(:,1), PointsArray(:,2));
                
                fx = F * PointsArray(:,1);
                fy = F * PointsArray(:,2);
                
                Px(:,i) = (E' * E) \ (E' * fx);
                Py(:,i) = (E' * E) \ (E' * fy);
                %Px(:,i) = inv(E) * fx;
                %Py(:,i) = inv(E) * fy;
            end
            i = min([NumFrames, (CurrentFrame + StepSize)]);
        else
            assert(size(Px,2) == NumFrames);
            disp('--- Skipping');
            i = NumFrames;
        end
        
        GivenRange = CurrentFrame - StepSize + 1:CurrentFrame;
        GivenRange = GivenRange(GivenRange > 0);
        
        GivenPx = Px(:,GivenRange);
        GivenPy = Py(:,GivenRange);
        
        CurrentRange = CurrentFrame + 1:i;
        
        if (DoOptimize)
            %options = optimset('PlotFcns', @optimplotfval, 'TolFun', 1e-6, 'TolX', 1e-4, 'MaxIter', Iterations);
            fun = @(PxPy) this.NonlinearMinimization(PxPy, GivenPx, GivenPy, GivenRange, CurrentRange, NN);
            
            PxPy0 = cat(3, Px(:,CurrentRange), Py(:,CurrentRange));
            PxPy = fminunc(fun, PxPy0(:), opts);
            %PxPy = fminsearch(fun, PxPy0(:), options);
            
            PxPy = reshape(PxPy, [], length(CurrentRange), 2);
            
            Px(:,CurrentRange) = PxPy(:,:,1);
            Py(:,CurrentRange) = PxPy(:,:,2);
        end
    end
end

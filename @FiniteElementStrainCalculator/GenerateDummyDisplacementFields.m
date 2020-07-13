%--------------------------------------------------------------------------
% GenerateDummyDisplacementFields:
%
% Generate dummy displacement fields for a frame sequence.
%
% INPUTS:
% PointGrid = Grid of points used to generate the displacement
%             field from.
% LoadFrameCallback = Callback function that accepts a frame number
%                     and returns that frame. e.i. frame <- f(frame_number)
%
% OUTPUT:
% No return. Sets the 'DisplacementFields' private member variable.
%--------------------------------------------------------------------------
function this = GenerateDummyDisplacementFields (this, PointGrid, NumFrames, noise_std, signal_level, LookAhead)
    P = PointGrid;

    for i = 1:NumFrames
        disp(['Processing frame ', num2str(i), ' of ', num2str(NumFrames)]);
        
        this.DisplacementFields{i,1} = [];
        this.DisplacementFields{i,2} = [];
        this.DisplacementFields{i,3} = [];
        
        
        Dx = ones(size(P,1), 1) .* signal_level;
        Dy = ones(size(P,1), 1) .* signal_level;
        
        Dx = Dx + normrnd(0, noise_std, size(P,1), 1);
        Dy = Dy + normrnd(0, noise_std, size(P,1), 1);
        
        Fx = scatteredInterpolant(P(:,1), P(:,2), Dx, 'natural');
        Fy = scatteredInterpolant(P(:,1), P(:,2), Dy, 'natural');
        
        for j = (i + 1):min([NumFrames, i + LookAhead])
            rms_error = 1; %sqrt(mean(error.^2));
            if (rms_error > this.ERROR_CUTOFF)
                disp(['-- Frame ', num2str(i), ' Exceeded error @ ', 'Frame ', num2str(j), '. (error = ', num2str(rms_error), ')']);
                break;
            end
            
            weight = 1 / rms_error;
            
            if (this.MaxWeight < weight)
                this.MaxWeight = weight; % keep max weight, to normalize the weight later.
            end

            if isempty(this.DisplacementFields{i,1})
                this.DisplacementFields{i,1}{1} = Fx;
                this.DisplacementFields{i,2}{1} = Fy;
                this.DisplacementFields{i,3}(1) = j;
                this.DisplacementFields{i,4}(1) = weight;
            else
                this.DisplacementFields{i,1}{end + 1} = Fx;
                this.DisplacementFields{i,2}{end + 1} = Fy;
                this.DisplacementFields{i,3}(end + 1) = j;
                this.DisplacementFields{i,4}(end + 1) = weight;
            end
        end
    end
end
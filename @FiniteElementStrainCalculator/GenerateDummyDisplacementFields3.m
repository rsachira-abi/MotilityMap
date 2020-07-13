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
function this = GenerateDummyDisplacementFields3 (this, PointGrid, NumFrames, noise_std, Disp)
    SNR = [];
    
    for i = 1:(NumFrames - 1)
        disp(['Processing frame ', num2str(i), ' of ', num2str(NumFrames)]);
        
        this.DisplacementFields{i,1} = [];
        this.DisplacementFields{i,2} = [];
        this.DisplacementFields{i,3} = [];
        
        P = PointGrid{i};
        
        for j = (i + 1):NumFrames
            Dx = Disp{i,j}(:,1);
            Dy = Disp{i,j}(:,2);
            
            NoiseX = normrnd(0, noise_std, size(P,1), 1);
            NoiseY = normrnd(0, noise_std, size(P,1), 1);
            
            Dx = Dx + NoiseX;
            Dy = Dy + NoiseY;
            
            SNR_ = 10 * log (sum(Dx.^2 + Dy.^2) / sum(NoiseX.^2 + NoiseY.^2));
            SNR = [SNR; SNR_];
            
            Fx = scatteredInterpolant(P(:,1), P(:,2), Dx, 'natural');
            Fy = scatteredInterpolant(P(:,1), P(:,2), Dy, 'natural');
        
            rms_error = 1; %sqrt(mean(error.^2));
            
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
    
    SNR = SNR(:);
    disp(['Mean SNR = ', num2str(mean(SNR)), ', std = ', num2str(std(SNR))]);
end
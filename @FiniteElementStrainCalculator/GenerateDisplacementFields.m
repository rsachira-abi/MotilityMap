%--------------------------------------------------------------------------
% GenerateDisplacementFields:
%
% Generate the displacement fields for a frame sequence.
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
function this = GenerateDisplacementFields (this, PointGrid, LoadFrameCallback, NumFrames)
    PointsArray = PointGrid;
    ImageSize = size(LoadFrameCallback(1));

    DownsampleRatio = pow2(floor(log2(round(max(ImageSize)/SavitskyGolayRegistration.CCSize))));
    for i=1:log2(DownsampleRatio)+1
        DownSamples(i) = DownsampleRatio/2^(i-1);
    end

    this.MaxWeight = 0;

    for i = 1:NumFrames
        disp(['Processing frame ', num2str(i), ' of ', num2str(NumFrames)]);

        for num = 1:length(DownSamples)
                DownSampledTargetImg = imresize(LoadFrameCallback(i), 1/DownSamples(num), 'bicubic');
                DownSampledPoints = ceil(PointsArray ./ DownSamples(num));
                ReferenceTemplates{num} = SavitskyGolayRegistration.GenerateTemplates(DownSampledTargetImg, DownSampledPoints, zeros(size(PointsArray)));
        end

        this.DisplacementFields{i,1} = [];
        this.DisplacementFields{i,2} = [];
        this.DisplacementFields{i,3} = [];

        for j = (i + 1):min([NumFrames, i + this.LOOK_AHEAD])
            [P, Dx, Dy, ~, error] = SavitskyGolayRegistration.RegisterImage(ReferenceTemplates, LoadFrameCallback(j), PointsArray, DownSamples);

            Fx = scatteredInterpolant(P(:,1), P(:,2), Dx, 'natural');
            Fy = scatteredInterpolant(P(:,1), P(:,2), Dy, 'natural');
            
            rms_error = sqrt(mean(error.^2));
            if (rms_error > this.ERROR_CUTOFF)
                disp(['-- Frame ', num2str(i), ' Exceeded error @ ', 'Frame ', num2str(j), '. (error = ', num2str(rms_error), ')']);
                %break;
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
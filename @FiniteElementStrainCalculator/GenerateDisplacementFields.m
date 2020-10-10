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
function this = GenerateDisplacementFields (this, PointsArray, LoadFrameCallback, NumFrames)
    this.MaxWeight = 0;
    
    template = LoadFrameCallback(1);
    reg = ImageRegistration(template, PointsArray);
    
    sz = reg.getTemplateSize();
    Queue = complex(double.empty(sz(1), sz(2), 0));
    QueueJ = int16.empty(1,0);
    qLength = 0;
    TemplateIndx = 0;
    for i = 1:NumFrames
        disp(['Frame: ', num2str(i)]);
        
        target = LoadFrameCallback(i);
        target_grad = reg.grad(target);
        [padded_target_grad, ~] = reg.padArray(target_grad, []);
        
        if (qLength < this.LOOK_AHEAD)
            qLength = qLength + 1;
            Queue(:,:,qLength) = padded_target_grad;
            QueueJ(1,qLength) = i;
            disp(['--- Added, length = ', num2str(qLength)]);
        end
        
        if (qLength == this.LOOK_AHEAD)
            disp('--- Processing');
            [Fx, Fy, error] = reg.calcDisplacementFields(Queue);
            weight = 1 ./ error;
            weight(weight > 1) = 1;
            max_weight = max(weight);
            
            this.DisplacementFields{TemplateIndx + 1, 1} = Fx;
            this.DisplacementFields{TemplateIndx + 1, 2} = Fy;
            this.DisplacementFields{TemplateIndx + 1, 3} = QueueJ;
            this.DisplacementFields{TemplateIndx + 1, 4} = weight;
            TemplateIndx = TemplateIndx + 1;
            
            if (this.MaxWeight < max_weight)
                this.MaxWeight = max_weight; % keep max weight, to normalize the weight later.
            end
            
            Queue(:,:,1) = [];
            QueueJ(:,1) = [];
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
		
		this.DisplacementFields{TemplateIndx + 1, 1} = Fx;
        this.DisplacementFields{TemplateIndx + 1, 2} = Fy;
        this.DisplacementFields{TemplateIndx + 1, 3} = QueueJ;
        this.DisplacementFields{TemplateIndx + 1, 4} = weight;
        TemplateIndx = TemplateIndx + 1;
		
		Queue(:,:,1) = [];
		QueueJ(:,1) = [];
		qLength = qLength - 1;
		disp(['--- Removed, length = ', num2str(qLength)]);
	end
end
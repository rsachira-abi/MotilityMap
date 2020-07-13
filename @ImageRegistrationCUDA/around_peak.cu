__global__ void around_peak (double *CC, const int *maxX, const int *maxY, const int SquareSize, const double peakThreshold, const int HalfCCSize, const int PointLength) {

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int i_ = blockIdx.z * blockDim.z + threadIdx.z;
	
	int CCSize = (2 * HalfCCSize) + 1;
	
	if (i_ >= PointLength || i >= CCSize || j >= CCSize)
		return;
	
	int j_min = maxX[i_] - SquareSize;
	int i_min = maxY[i_] - SquareSize;
	int j_max = maxX[i_] + SquareSize;
	int i_max = maxY[i_] + SquareSize;
	
	int indx = (i_ * CCSize * CCSize) + (j * CCSize) + i;
	
	if (i >= i_min && i <= i_max && j >= j_min && j <= j_max) {
		
		double CCval = CC[indx];
		
		if (CCval > peakThreshold)
			CC[indx] = 1;
		else
			CC[indx] = 0;
	}
	else {
		CC[indx] = 0;
	}
}
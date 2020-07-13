__global__ void windowROI (double2 *roi_1, double2 *roi_2, const double2 *padded_template,
								const double2 *padded_target, const double *window2D, const int *points,
								const int *dispX, const int *dispY, const int HalfCCSize, const int PointLength, const int ImageWidth, const int ImageHeight) {

	
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int i_ = blockIdx.z * blockDim.z + threadIdx.z;
	
	int CCSize = (2 * HalfCCSize) + 1;
	
	if (i_ >= PointLength)
		return;
	
	int x1 = points[i_] - HalfCCSize - 1 + j;
	int y1 = points[i_ + PointLength] - HalfCCSize - 1 + i;
	int x2 = points[i_] - HalfCCSize - 1 + j + dispX[i_];
	int y2 = points[i_ + PointLength] - HalfCCSize - 1 + i + dispY[i_];
	
	if (i >= CCSize || j >= CCSize)
		return;
	
	double2 r1, r2;
	
	double2 I1 = padded_template[x1*ImageHeight + y1];
	double2 I2 = padded_target[x2*ImageHeight + y2];
	
	r1.x = window2D[j*CCSize + i] * I1.x;
	r1.y = window2D[j*CCSize + i] * I1.y;
	r2.x = window2D[j*CCSize + i] * I2.x;
	r2.y = window2D[j*CCSize + i] * I2.y;
	
	int indx = (i_ * CCSize * CCSize) + (CCSize * j) + i; //(i * CCSize * PointLength) + (j * PointLength) + i_;
	roi_1[indx] = r1;
	roi_2[indx] = r2;
}
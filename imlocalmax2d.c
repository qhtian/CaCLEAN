/***********************************************************************************
 * imlocalmax2d.c
 *
 * Compilation:
 *   mex -v -g -largeArrayDims imlocalmax2d.c
 * Matlab Calling:
 *   p=imlocalmax2d(double(I),logical(h));
 * By Qinghai Tian.
 **********************************************************************************/

 /***********************************************************************************/
#include <matrix.h>
#include <mex.h>
#include <math.h>
// #include <omp.h>

// #ifdef NAN
	/* NAN is supported */
// #endif
#ifdef INFINITY
/* INFINITY is supported */
#endif

// Main funciton is here.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *I, *II, *pOutputXY, Imax, *ImaxArray;
	mxLogical *H;
	size_t NumOfElements, mRow, nCol, NumOfFramePixelNew, NumOfFramePixel;
	size_t mRowH, nColH, NumOfElementsOfH, NumOfMask = 0, NumOfFramePixelH;
	size_t i, j, m, n, mRowHhalf, nColHhalf, currPos, currPosH, mRowNew;
	size_t NumOfDim, NumOfDimH, mRowXnCol;

	// Check inputs.
	if (nrhs != 2) mexErrMsgIdAndTxt("MATLAB:imlocalmax2d:NumOfInput", "Number of input argument must be 2. Im=imlocalmax2d(I,H);");
	if (nlhs > 1) mexErrMsgIdAndTxt("MATLAB:imlocalmax2d:NumOfOutput", "Number of output argument must be <= 1.");
	if (!(mxIsDouble(prhs[0])))	mexErrMsgIdAndTxt("MATLAB:imlocalmax2d:DataType", "First input image must be double.");

	if (!(mxIsLogical(prhs[1]))) mexErrMsgIdAndTxt("MATLAB:imlocalmax2d:NumOfInput", "Second input image must be logical.");

	NumOfDim = mxGetNumberOfDimensions(prhs[0]);
	NumOfDimH = mxGetNumberOfDimensions(prhs[1]);

	I = mxGetPr(prhs[0]);
	H = mxGetLogicals(prhs[1]);

	NumOfElements = mxGetNumberOfElements(prhs[0]);
	NumOfElementsOfH = mxGetNumberOfElements(prhs[1]);

	mRow = mxGetM(prhs[0]);
	nCol = mxGetN(prhs[0]);

	mRowH = mxGetM(prhs[1]);
	nColH = mxGetN(prhs[1]);

	NumOfFramePixelH = mRowH*nColH;

	if ((mRowH < 3) || (nColH < 3)) mexErrMsgIdAndTxt("MATLAB:imlocalmax2d:NumOfInput", "Second input image must be logical.");

	mRowHhalf = mRowH / 2;
	nColHhalf = nColH / 2;

	for (i = 0; i < NumOfElementsOfH; i++) {
		if (H[i]) {
			NumOfMask = NumOfMask + 1;
		}
	}
	if (NumOfMask < 1) mexErrMsgIdAndTxt("MATLAB:imlocalmax2d:ValidPixelMask", "Second logical mask must contains one pixel.");

	II = (double*) mxMalloc((mRow + mRowH * 2) * (nCol + nColH * 2) * sizeof(double));
	if (II == NULL) return;
	for (i = 0; i < (mRow + mRowH * 2) * (nCol + nColH * 2); i++) 	II[i] = -INFINITY; // Initialize array.

	ImaxArray = (double*)mxMalloc(mRow * nCol * sizeof(double));
	if (ImaxArray == NULL) {
		mxFree(II);
		return;
	}

	// Pad array.
	mRowNew = mRow + mRowH * 2;
	for (i = mRowH; i < (mRow + mRowH); i++) {
		for (j = nColH; j < (nCol + nColH); j++) {
			II[j*mRowNew + i] = I[(j - nColH)*mRow + (i - mRowH)];
		}
	}


	// Cycle every elements.
	for (j = nColH; j < (nCol + nColH); j++) {
		for (i = mRowH; i < (mRow + mRowH); i++) {
			// Current Pixel's surroundings.
			Imax = -INFINITY;
		 	for (n = 0; n < nColH; n++) {
		 		for (m = 0; m < mRowH; m++) {
					currPos = (j + n - nColHhalf)*mRowNew + (i + m - mRowHhalf);
					currPosH = n*mRowH + m;
					if((H[currPosH])&&(II[currPos]>Imax))	Imax=II[currPos];
		 		}
		 	}
			currPos = (j -nColH)*mRow + (i - mRowH);
			ImaxArray[currPos] = Imax;
		}
	}

	// Output.
	if (nlhs == 1) {
		plhs[0] = mxCreateNumericMatrix(mRow, nCol, mxDOUBLE_CLASS, mxREAL);
		pOutputXY = mxGetPr(plhs[0]); // create a C pointer to a copy of the output matrix.
		for (n = 0; n < NumOfElements; n++) pOutputXY[n] = ImaxArray[n];
	}

	// Clear allowcated memories.
	mxFree(II);
	mxFree(ImaxArray);
}

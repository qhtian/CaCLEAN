/* 	CaCLEANXYKernel, find global maxima, find its center and CLEAN it. Do this until
		global maxima is below the given threshold.

	Calling from Matlab, e.g.:
		xy=CaCLEANXYKernel(double(I),bw,th,double(cleanPSF(:)),[maxIter]);
            I: 			will be modified and becomes the residual.
			th:			1 element, double;
			cleanPSF: 	a CLEAN PSF in a column format.
			maxIter:	optional, 1 element, double; default value = 1e6;

			xy:			output of centers of CLEAN cleanPSF.

	Compile with this command: 	mex -largeArrayDims CaCLEANXYKernel.c

    Author: Qinghai Tian (tian_qhcn@icloud.com).
    www.lipplab.de
    For research use only.
*/

#include <math.h>
#include <time.h>
#include <matrix.h>

// main mex entrance function here.
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *I, *Im, *threshold, *PSFInfo, I_Max, xyCenter[2]={0.0,0.0}, *maxIterNo, maxIterDouble;
	double *pOutputXY, *pOutputResidual;
	unsigned long long	pixelNum, currPos, *pixelIndexArray, i, j, PSFpieceX0, PSFpieceY0,ii,
						NumberOfElements=0, PSFRadius=0,PSFsize=0, mrows, ncols, i_MaxPos, j_MaxPos,
						pixelNumStep, pixelNumDone, pixelNumNew, IterationElapsed, maxIteration=0,currMaxPos;
	unsigned long long *CaRelease;
	signed long long currMaxPosIndex;
	mxLogical *bw;

    /////////////////////////    Check inputs.   /////////////////////////
	if ((mxGetM(prhs[1]) != mxGetM(prhs[0])) || (mxGetN(prhs[1]) != mxGetN(prhs[0])))
    	return;

    // Create a pointer to the input matrix. prhs: parameter of right hand side.
    I         = mxGetPr(prhs[0]);
    bw  = mxGetLogicals(prhs[1]);
    threshold = mxGetPr(prhs[2]);
    PSFInfo   = mxGetPr(prhs[3]);
    maxIteration=(unsigned long long)1e7;
    if(nrhs==5){
		maxIterNo=mxGetPr(prhs[4]);
        maxIterDouble=*maxIterNo;
		maxIteration=(unsigned long long)(maxIterDouble);
	}

	// Size of the input image.
    mrows = mxGetM(prhs[0]);  // First dimensions of the matrix input matrix.
    ncols = mxGetN(prhs[0]);  // Second dimensions of the matrix input matrix.

	// Copy the input image for later modification.
	Im = (double *) mxMalloc(sizeof(double) * mrows * ncols);
	if (Im==NULL) return;
	for(i=0; i<mrows * ncols; i++)
		Im[i] = I[i];

	CaRelease=(unsigned long long *) mxMalloc(sizeof(unsigned long long) * mrows * ncols);
	if (CaRelease==NULL) {
		mxFree(Im);
		return;
	}
	for(i=0; i<mrows * ncols; i++)
		CaRelease[i] = 0;

	// Size of PSF.
    NumberOfElements = mxGetNumberOfElements(prhs[3]);
    PSFsize = (unsigned int) sqrt((int)NumberOfElements);
    PSFRadius = (PSFsize - 1)/2;

    /* Create one column to store index of every single pixel in the image. */
    pixelIndexArray = (unsigned long long *) mxMalloc(sizeof(unsigned long long) * mrows * ncols);
    if (pixelIndexArray==NULL) {
    	mxFree(Im);
		mxFree(CaRelease);
    	return;
	}

    // Index all positive pixels to speed up.
	pixelNum=0;
    for(j=0;j<ncols;j++) {       // Dim 2.
        for(i=0;i<mrows;i++) {   // Dim 1.
			currPos = j * mrows + i;
            if (bw[currPos]) {
                if (Im[currPos]>*threshold) {
                    pixelIndexArray[pixelNum]=currPos;
                    pixelNum=pixelNum+1;
                }
            }
        }
    }


	/* ========== the main while loop of The CLEAN Algorithm is here. =========== */
	pixelNumStep = pixelNum / 50;
	pixelNumDone = 0;
	IterationElapsed = 0;

	while(IterationElapsed<maxIteration)
	{
        //=============== Search for the maxima and its indexing. ===============
        I_Max=*threshold;//Im[pixelIndexArray[0]];;
        currMaxPosIndex = -1;
        for(ii=0; ii<pixelNum; ii++)
        {
            if(Im[pixelIndexArray[ii]]>I_Max)
            {
                I_Max=Im[pixelIndexArray[ii]];
                currMaxPosIndex = ii;
            }
		}

		//====== Clean one small piece of ACO (like a PSF) from the sample. =====
		if(currMaxPosIndex>=0) //(I_Max>*threshold)
		{
			currMaxPos = pixelIndexArray[currMaxPosIndex];
			CaRelease[currMaxPos]=CaRelease[currMaxPos]+1;

			i_MaxPos = currMaxPos % mrows;  // Dim1
        	j_MaxPos = currMaxPos / mrows;  // Dim2

			// Find out the area of the maxima position and subtract one copy of ACO.
			PSFpieceX0 = i_MaxPos - PSFRadius;
			PSFpieceY0 = j_MaxPos - PSFRadius;

			for(j=0; j<PSFsize; j++) {
				for(i=0; i<PSFsize; i++) {
					currPos = (PSFpieceY0 + j) * mrows + (PSFpieceX0 + i);
                    Im[currPos] = Im[currPos] - PSFInfo[i+j*PSFsize];
				}
			}

			// Check whether it is possible to reduce the cycling number of pixels
			if(Im[currMaxPos]<=*threshold)
			{
				pixelNumDone = pixelNumDone + 1;
				if (pixelNumDone > pixelNumStep)
				{
					pixelNumNew=0;
					pixelNumDone=0;
					for(j=0;j<pixelNum;j++)
					{
						if(Im[pixelIndexArray[j]]>*threshold)
						{
							pixelIndexArray[pixelNumNew] = pixelIndexArray[j];
							pixelNumNew = pixelNumNew + 1;
						}
					}
					pixelNum = pixelNumNew;
				}
			}
		}
		else {
			maxIteration = IterationElapsed;  // Make sure that next iteration will be stopped.
		}

		IterationElapsed = IterationElapsed + 1;
	}

    // Create Output Matrix.  // plhs: parameter of left hand side.
    plhs[0] = mxCreateDoubleMatrix((unsigned long long) mrows, ncols, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((unsigned long long) mrows, ncols, mxREAL);

    pOutputXY = mxGetPr(plhs[0]); // create a C pointer to a copy of the output matrix.
	pOutputResidual =  mxGetPr(plhs[1]);

	for(i=0; i<mrows*ncols; i++) {
		pOutputXY[i]=(double)CaRelease[i];
		pOutputResidual[i]=Im[i];
	}

    // clear new allocated memory
    mxFree(pixelIndexArray);
	mxFree(CaRelease);
	mxFree(Im);
}

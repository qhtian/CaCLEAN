/***********************************************************************************
 * isodata.c
 *
 * Compilation:
 *   mex -v -largeArrayDims isodata.c
 * Matlab Calling:
 *   p=isodata(double(I));
 * IEEE Transactions on Systems, Man, and Cybernetics, vol. SMC-8, No. 8, 630, 1978
 * Picture Thresholding using an iterative selection method. T.W.Ridler and S. Calvard
 **********************************************************************************/

// MATLAB Primitive Types
// mxClassID Value     MATLAB Type     MEX Type	C Primitive Type
// mxINT8_CLASS        int8            int8_T      char, byte
// mxUINT8_CLASS       uint8           uint8_T     unsigned char, byte
// mxINT16_CLASS       int16           int16_T     short
// mxUINT16_CLASS      uint16          uint16_T    unsigned short
// mxINT32_CLASS       int32           int32_T     int
// mxUINT32_CLASS      uint32          uint32_T    unsigned int
// mxINT64_CLASS       int64           int64_T     long long
// mxUINT64_CLASS      uint64          uint64_T    unsigned long long
// mxSINGLE_CLASS      single          float       float
// mxDOUBLE_CLASS      double          double      double



/***********************************************************************************/
#include <matrix.h>
#include <mex.h>
#include <math.h>

#ifdef NAN
/* NAN is supported */
#endif
#ifdef INFINITY
/* INFINITY is supported */
#endif



// Main funciton is here.
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double    *I, th1, th2, th_old, th_new, tol = 1e-6;
    char      *Ic;
    unsigned char *Iuc;
    unsigned short *Ius;
    unsigned int *Iui;
    short     *Is;
    int       *Ii;
    long long *Ill;
    unsigned long long *Iull;
    float     *If;

    long double sum1, sum2;
    size_t    i, maxiter = 100, iter = 0, NumOfElments, bwNum, memClear=0, elementSize=0;
    
    mxClassID   category;
    
    if(nrhs != 1){
        mexErrMsgIdAndTxt( "MATLAB:isodata:NumOfInput",
                "Number of input argument must be 1.");
    }
    
    if(nlhs > 1){
        mexErrMsgIdAndTxt( "MATLAB:isodata:NumOfOutput","Number of output argument must be <=1");
    }
    
    NumOfElments = mxGetNumberOfElements(prhs[0]);
    
    // Convert to double for processing.
    if(!mxIsDouble(prhs[0])){
        I = (double *) mxMalloc(sizeof(double) * NumOfElments);
        if(I==NULL) return;
        memClear=1;
        category = mxGetClassID(prhs[0]);
        switch (category){
            case mxUINT8_CLASS:  Iuc=mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Iuc[i];  break;
            case mxUINT16_CLASS: Ius=mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Ius[i];  break;
            case mxUINT32_CLASS: Iui=mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Iui[i];  break;
            case mxINT8_CLASS:   Ic =mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Ic[i];   break;
            case mxINT16_CLASS:  Is =mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Is[i];   break;
            case mxINT32_CLASS:  Ii =mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Ii[i];   break;
            case mxINT64_CLASS:  Ill=mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)Ill[i];  break;
            case mxUINT64_CLASS: Iull=mxGetData(prhs[0]); for (i=0; i<NumOfElments; i++) I[i] = (double)Iull[i]; break;
            case mxSINGLE_CLASS: If =mxGetData(prhs[0]);  for (i=0; i<NumOfElments; i++) I[i] = (double)If[i];   break;
            default: mexErrMsgIdAndTxt( "MATLAB:isodata:DataTypeInput","Input data type is not supported."); 
                if (memClear==1){ mxFree(I); return;} break;
        }
    }
	else{
		I = mxGetPr(prhs[0]);
	}
        

    // Main body of cycling.
    sum1 = 0;
    for (i=0; i<NumOfElments; i++){
        if ((!isnan(I[i])) & (!isinf(I[i])))
            sum1 = sum1 + I[i];
    }
    th_old = sum1 / (long double) NumOfElments;

    while (iter < maxiter){
        iter = iter + 1;
        sum1 = 0;
        sum2 = 0;
        bwNum = 0;
        for (i=0; i<NumOfElments; i++){
            if ((!isnan(I[i])) & (!isinf(I[i]))) {
                if(I[i]<=th_old){
                    sum1 = sum1 + I[i];
                    bwNum = bwNum + 1;
                }
                else{
                    sum2 = sum2 + I[i];
                }
            }
        }
        th1 = sum1 / (long double) bwNum;
        th2 = sum2 / (long double) (NumOfElments - bwNum);
        th_new = (th1 + th2) / 2;
        if(fabs(th_new-th_old) < tol){
            break;
        }
        else{
            th_old = th_new;
        }
    }

    if (memClear==1){
        mxFree(I);
    }
    
    if(nlhs == 1)	{
		plhs[0] = mxCreateDoubleScalar(th_new);
	}
	else{
		mexPrintf("   isopoint = %-4.5f\n", th_new);
	}
}
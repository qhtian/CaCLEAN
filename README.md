CaCLEAN algorithm

In beating cardiomyocytes, synchronized localized Ca2+ release from clusters of ryanodine receptors (RyR) is amongst the most important determinants of cellular contractility and the heart´s output. Fundamental properties including the release reliability of such individual clusters can be particularly important to understand the development of cardiac diseases, e.g. arrhythmia and heart failure. Novel revolution analysis methods, the CaCLEAN algorithm, were developed to understand the behavior of Ca2+ releasing RyR clusters in beating cardiac myocytes by combining the astronomical CLEAN algorithm with known properties of Ca2+ diffusion. 

Function Lists and notes:
  1. CICRcleanSimp: CICRcleanSimp calculates the Calcium release map, and it is the kernel funciton of the CaCLEAN algorithm.
  2. CICRrebuildSimp: CICRrebuildSimp calculates the upstroke of a calcium transient from the calculated calcium release map that is derived with CaCLEAN algorithm.
  3. CICRsimulation: CICRsimulation simulates confocal recordings of cardiac calcium transient.
  4. CRUProps: CRUProps segments the calcium release map and calculates the properties of single Calcium Release Units (CRU).
  5. generatePSF: Internal function to generate Analytical CLEAN Object (ACO). An ACO is similar to a Point Spread Function.
  6. SampleCell: SampleCell genarates a sample ventricular myocyte.
  7. sCMOSdynamicOffsetCorrection: sCMOSdynamicOffsetCorrection removes dynamic offset of sCMOS camera.
  8. CICRcleanMultiCaTBatchSimp: the batch mode for the CaCLEAN, which reads the raw data into matlab, calculates the Ca2+ release maps and saves all the results automatically. This batch mode was customized to fit our workflow as described in the eLife paper. Because the source codes uses functions from other sources, and due to legal reasons this function will NOT be published here. Upon request, a stripped version can be provided.
  
Demonstration of the CaCLEAN algorithm:
  Demo_with_SampleData.m is a short demonstration script to run the CaCLEAN progam. Sample data is included the SampleData.mat file.

Before the MatLab functions can be used, please be noted:

  1. CRUProps function needs fminsearchbnd function by John D'Errico (woodchips@rochester.rr.com). It is available in the File Exchange of MathWorks. Please download it and put it into current folder containing all the functions.
  
  2. CICRcleanSimp needs smoothn function: Damien Garcia, Robust smoothing of gridded data in one and higher dimensions with missing values. Comput Stat Data Anal. 2010 Apr 1; 54(4): 1167–1178. Please download it and put it into current folder containing all the functions.
  
  3. Some kernel functions are programed in C programming language with MatLab/MEX interfaces. Before you can use them, please compile these functions in MatLab Command Window first:
  
	>> mex -largeArrayDims CaCLEANXYKernel.c  % For CICRcleanSimp function. This is the kernel loop funciton.
	>> mex -largeArrayDims imlocalmax2d.c     % For CRUProps funciton.
	>> mex -largeArrayDims isodata.c          % Isodata algorithm, useful image segmentation utility.
	
  4. The program package here does NOT contain denosing program CANDLE.
  
  5. Reading data, e.g. tiff files into MatLab, Loci Bio-Formats plugin for ImageJ and toolbox for MatLab (http://www.openmicroscopy.org/bio-formats/downloads/) were used. Please refer to the documentation for detailed information.


Please refer to the original paper for detailed information:
Tian Q, et. al., An adaptation of astronomical image processing enables characterization and functional 3D mapping of individual sites of excitation-contraction coupling in rat cardiac muscle, eLife, 2017, in press.


Qinghai Tian (tian_qhcn@icloud.com)
Center for Molecular Signaling (PZMS)
Institute for Molecular Cell Biology
Research Center for Molecular Imaging and Screening
Building 61, Medical Faculty
Saarland University Hospital
Saarland University
66421 Homburg/Saar, Germany
www.pzms.uni-saarland.de   -or-   www.lipplab.de

%% Load sample data.
load('SampleData.mat')
% RawUpstroke is the raw data.
% Idenoised is the denoised version of the raw upstroke from Figure 1.
% Bgr is the denoised basal fluorescence.
% Membrane is the signal recorded with CellMask DeepRed and deconvolved with AutoDeblur (25 Cycles).
% Mask is the cell mask that labels out the entire cell.
% xyt_dim is the physical dimension of a single pixel in [x(um), y(um), t(ms)];

%% Do the CaCLEAN here.
CleanObj=CICRcleanSimp(Idenoised,Bgr,Mask,xyt_dim,'ApparentDiffusionK',60,'CleanDiffusionK',30,'CaCleanThreshold',10);

%% Display the sample data.
subplot(1,3,1)
imagesc(CleanObj.CaRelease2D); axis image off; caxis([0,1500])
title('CaCLEANed CRU map')

subplot(1,3,2)
imagesc(Membrane); axis image off; caxis([0,2000])
title('Deconvolved Cell Membrane')

G=CleanObj.CaRelease2D;
G=G/1500; G(G<0)=0; G(G>1)=1;

R=Membrane; R=R/2000;
R(R<0)=0; R(R>1)=1;


subplot(1,3,3)
imagesc(cat(3,R,G,zeros(size(G)))); axis image off
title('Merged Ca2+ Release Map and Cell Membrane')


%% Rebuild the upstroke from the CaCLEAN results and put it in the field "CICRrebuilt".
CleanObj=CICRrebuildSimp(CleanObj);
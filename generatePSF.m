function [PSF,CaDiffuse]=generatePSF(xyt_dim,ApparentDiffusionK)
    %% [PSF,CaDiffuse]=generatePSF(xyt_dim,ApparentDiffusionK);
    % Internal function to generate Analytical CLEAN Object (ACO). An ACO
    % is similar to a Point Spread Function.
    %
    % [PSF,CaDiffuse]=generatePSF(xyt_dim,ApparentDiffusionK);
    %
    % Qinghai Tian
    % Institute for Molecular Cellbiology
    % Medical Facalty of University of
    % Saarland.
    % Homburg, Germany.
    % tian_qhcn@icloud.com
    
    CaDiffuse=@(t,D,x,y) exp(-(x.^2+y.^2)/(4*D*t));
    CaDiffuseSize=xyt_dim(1):xyt_dim(1):30;
    CaDiffuseSize=[sort(-CaDiffuseSize,'ascend'),0,CaDiffuseSize];
    
    if mod(numel(CaDiffuseSize),2)==0
        CaDiffuseSize=numel(CaDiffuseSize)/2;
        CaDiffuseSize=0:xyt_dim(1):xyt_dim(1)*CaDiffuseSize;
        CaDiffuseSize=cat(2,-fliplr(CaDiffuseSize(2:end)),CaDiffuseSize);
    end
    [x,y]=ndgrid(CaDiffuseSize,CaDiffuseSize);
    PSF=CaDiffuse(xyt_dim(3)/1000,ApparentDiffusionK,x,y);
    PSF_bw=PSF>max(PSF(:))*0.001;
    PSF_bw=sum(PSF_bw,2);
    PSF_bw=ceil(sum(PSF_bw>0)/2);
    PSF_siz=(size(PSF,1)+1)/2;
    PSF=PSF((PSF_siz-PSF_bw):(PSF_siz+PSF_bw),(PSF_siz-PSF_bw):(PSF_siz+PSF_bw));
    PSF=PSF/sum(PSF(:));
end
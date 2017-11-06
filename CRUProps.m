function CleanObj=CRUProps(CleanObj)
    %% CleanObj=CRUProps(CleanObj);
    % CRUProps segments the calcium release map and calculates the
    % properties of single Calcium Release Units (CRU).
    %
    % Qinghai Tian
    % Institute for Molecular Cellbiology
    % Medical Facalty of University of
    % Saarland.
    % Homburg, Germany.
    % tian_qhcn@icloud.com
    
    %% Some parameters.
    xyt_dim=CleanObj(1).xyt_dim;
    
    % Some temporary parameters;
    minimumPixelNumPerCluster=0.2/xyt_dim(1)/xyt_dim(2); % um^2, 2.8 Pixels
    
    clusterTerritoryRadius=ceil(0.5/xyt_dim(1));          % 1.5 um.
    localThrsholdNormalizedToLocalMax=0.05;
    
    
    %% Calculation start here.
    hLocalMax=fspecial('disk',clusterTerritoryRadius);
    hLocalMax=hLocalMax==max(hLocalMax(:));
        

    for k=1:numel(CleanObj)
        % Calculate the dF/F0
        currMask=CleanObj(k).Mask;
        I=double(CleanObj(k).CaRelease2D)./double(CleanObj(k).DataBgr);
        I=I.*(I>0).*currMask;
        I(isnan(I))=0;
        I=I.*(I>CleanObj(k).CaCleanThreshold/mean(CleanObj(k).DataBgr(CleanObj(k).Mask))); % Just to remove values close to zero.
        
        
        % Local max mask.
        localMaxMask=imlocalmax2d(I,hLocalMax);
        localMaxMask=(I>(localMaxMask*localThrsholdNormalizedToLocalMax));
        
        % Watershed
        CRULabel=watershed(-I);
        CRULabel=single(CRULabel).*currMask.*localMaxMask;

        CRULabelProps=regionprops(CRULabel,I,'BoundingBox','WeightedCentroid','Area');
        CRUnum=numel(CRULabelProps);
        
        FWHM=nan(numel(CRULabelProps),1);
        Amp=nan(numel(CRULabelProps),1);
        Isiz=size(I);
        parfor j=1:CRUnum
            if CRULabelProps(j).Area<minimumPixelNumPerCluster; continue; end
            
            currProps=CRULabelProps(j).BoundingBox;
            x1=ceil(currProps(2));     if x1<1; x1=1; end
            x2=floor(x1+currProps(4)); if x2>Isiz(1);x2=Isiz(1); end %#ok<PFBNS>
            y1=ceil(currProps(1));     if y1<1; y1=1; end
            y2=floor(y1+currProps(3)); if y2>Isiz(2);y2=Isiz(2); end
            
            currBW=CRULabel(x1:x2,y1:y2)==j; %#ok<PFBNS>

            currCRU=I(x1:x2,y1:y2); %#ok<PFBNS>
            currCRU=currCRU.*currBW;
            [FWHM(j),Amp(j)]=PeakFitting(currCRU,currBW);
        end
        
        
        for j=CRUnum:-1:1
            if isnan(FWHM(j))
                CRULabel(CRULabel==j)=0;
                CRULabelProps(j)=[];
                Amp(j)=[];
                FWHM(j)=[];
            else
                CRULabelProps(j).FWHM=FWHM(j)*xyt_dim(1);
                CRULabelProps(j).Amp=Amp(j);
            end
        end

        CleanObj(k).CaRelease2D_dFF0=I;
        CleanObj(k).CRUProps=CRULabelProps;
        CleanObj(k).CRULabel=CRULabel;
        
        
        CleanObj(k).CV_IntraCaTClusterFiring=std(Amp)/mean(Amp)*100;
        fprintf('\t%d\tNumCluster=%0.0f\n',k,numel(CRULabelProps));
    end
    
    CV_InterCaTClusterFiring=calculateCVinterCaT(cat(3,CleanObj.CRULabel),cat(3,CleanObj.CaRelease2D_dFF0));
    for k=1:numel(CleanObj)
        CleanObj(k).CV_InterCaTClusterFiring=CV_InterCaTClusterFiring(k);
    end
end


function [FWHM,Amp]=PeakFitting(I,mask)
    I=double(I);
    [sizey,sizex] = size(I);
    
%% Get center of mass, amplitude, and sigma.
    [X,Y]=ndgrid(1:sizey,1:sizex);
    bwMaxLoc=(I==max(I(:)));
    cx=X(bwMaxLoc); if numel(cx)>1; cx=cx(1); end
    cy=Y(bwMaxLoc); if numel(cy)>1; cy=cy(1); end
    
    distance=sqrt((X-cx).^2+(Y-cy).^2); 
    sigma=sqrt(sum(mask(:))/pi)/2;

    distance=distance(mask);
    I=I(mask);

    I_max=max(I);
    Dis_max=max(distance);

%% Do a Gaussian fitting with mu=0 and Bgr=0.
    Gau=@(x,p)p(1)*exp(-x.^2/2/p(2)^2);
    fun_dev=@(x,y,p)sum((Gau(x,p)-y).^2);
    options=optimset('MaxIter',100000000,'Display','off');
    p=fminsearchbnd(@(p)fun_dev(distance,I,p),[I_max sigma],...
        [I_max/10 sigma/10], [I_max*20 Dis_max*3],options);
    FWHM=p(2)*2;
    Amp=p(1);
end

function CV=calculateCVinterCaT(ClusterLabel,I)
    Isiz=size(ClusterLabel,3);
    CV=nan(Isiz,1);
    bw=ClusterLabel>0;
    for k=1:Isiz
        currBW=bw(:,:,k);
        ImProps=regionprops(currBW,I(:,:,k),'Area','MeanIntensity');
        ImProps=[ImProps(:).MeanIntensity]'.*[ImProps(:).Area]';
        ImProps=repmat(ImProps,[1,Isiz]);
        for j=1:Isiz
            if j==k; continue; end
            currImProps=regionprops(currBW,I(:,:,j),'Area','MeanIntensity');
            currImProps=[currImProps(:).MeanIntensity]'.*[currImProps(:).Area]';
            ImProps(:,j)=currImProps;
        end
        currCV=std(ImProps,0,2)./mean(ImProps,2)*100;
        CV(k)=mean(currCV);
    end 
end
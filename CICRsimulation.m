function varargout=CICRsimulation(varargin)
    %% varargout=CICRsimulation(varargin);
    % CICRsimulation simulates confocal recordings of cardiac calcium transient.
    %
    % Name-Value parameters:
    %   Please check the parameter defininaitons in the first section of
    %   the program body.
    % Output:
    %   SparkWithNoise:    Spark recordings with noise.
    %   SparkNoiseFree:    Spark recordings without noise.
    %   SparkPosition:     (1)ID (2)xc (3)yc (4)t_onset (5)Bgr (6)dF/F0
    %                      (7)FWHM/2 (8)Decay.
    % Qinghai Tian
    % Institute for Molecular Cellbiology
    % Medical Facalty of University of
    % Saarland.
    % Homburg, Germany.
    % tian_qhcn@icloud.com
    
    %% Input.
    CPU_NumCores=feature('NumCores');
    p=inputParser;
    
    p.addParameter('LaserIntensity',30,@(x)isscalar(x) && x>0 && x<=100);  % Percentage.
    
    p.addParameter('DetectorOffset',400,@(x)isscalar(x) && x>0);           % a.u.
    p.addParameter('DetectorGain',2.7,@(x)isscalar(x) && x>0);
    p.addParameter('DetectorGaussNoise',6.2,@(x)isscalar(x) && x>0);       % a.u.
    p.addParameter('xyt_dim',[0.215,0.215,6.85],@(x)numel(x)==3);          % In micrometer.
    
    p.addParameter('CaReleaseSigma',1,@(x)isscalar(x) && x>0);             % In millisecond.
    p.addParameter('CaReleaseTau',  3,@(x)isscalar(x) && x>0);             % In millisecond.
    p.addParameter('CaReleaseFWHMMax',0.2,@(x)isscalar(x) && x>0);         % In micrometer.
    p.addParameter('CaReleaseTauFWHM',1,@(x)isscalar(x) && x>0);           % In millisecond.
    p.addParameter('CaReleaseDurationSD',3,@(x)isscalar(x) && x>0);
    
    p.addParameter('Xdim',[],@(x)isscalar(x) && x>0);                      % In pixel.
    p.addParameter('Ydim',[],@(x)isscalar(x) && x>0);                      % In pixel.
    p.addParameter('Tdim',600,@(x)isscalar(x) && x>0);                     % In ms.
    p.addParameter('ReleaseAmplitude',[],@(x)isempty(x) || (isscalar(x) && x>0));
    
    p.addParameter('ApparentDiffusionK',50,@(x)isscalar(x) && x>0);        % In millisecond.
    p.addParameter('DecayFF0Rate',0.03,@(x)isscalar(x) && x>0);            % In F = F - dF/F0 * K * t.
    
    p.addParameter('CaReleaseNum',35000,@(x)isscalar(x) && x>=0);
    
    p.addParameter('ExpAmp', false, @(x)islogical(x));                     % Distribution.
    p.addParameter('RampAmp',false, @(x)islogical(x));                     % Distribution.
    p.addParameter('CaReleaseAmp',1.0,@(x)isscalar(x) && x>0);               % In dF/F0.
    
    p.addParameter('CPUNumCores',CPU_NumCores,@(x)isscalar(x) && x>0 && x<CPU_NumCores);
    
    % Parameters for the FWHM function.
    parse(p, varargin{:});
    p=p.Results;
    clear('varargin')
    
    %% Parameters and sample Ca release.
    Gain=p.DetectorGain;
    xyt_dim=p.xyt_dim;
    Xdim=p.Xdim;
    Ydim=p.Ydim;
    Tdim=ceil(p.Tdim/xyt_dim(3));
    ReleaseTrigger=Tdim*0.3;
    LaserIntensity=p.LaserIntensity;
    CaReleaseAmp=p.CaReleaseAmp;
    CaReleaseDurationSD=p.CaReleaseDurationSD;
    CaReleaseNum=p.CaReleaseNum;
    CaReleaseSigma=p.CaReleaseSigma;
    CaReleaseTau=p.CaReleaseTau;
    CaReleaseFWHMMax=p.CaReleaseFWHMMax;
    CaReleaseTauFWHM=p.CaReleaseTauFWHM;
    ApparentDiffusionK=p.ApparentDiffusionK;
    DetectorOffset=p.DetectorOffset;
    DetectorGaussNoise=p.DetectorGaussNoise;
    ReleaseAmplitude=p.ReleaseAmplitude;
    DecayFF0Rate=p.DecayFF0Rate;

    
    
    SampleCaReleaseSpark=SparkDsptSimAllInOneV1_2('Mu',xyt_dim(3),'SpR',xyt_dim(1),'TpR',xyt_dim(3),...
        'FWHMMax',CaReleaseFWHMMax,'TauFWHM',CaReleaseTauFWHM,...
        'Sigma',CaReleaseSigma,'Tau',CaReleaseTau);
    clear('CaReleaseFWHMMax','CaReleaseTauFWHM','CaReleaseSigma','CaReleaseTau');
    
    %% Load sample bgr.
    if isempty(Xdim)
        S=SampleCell;
        Xdim=size(S.Bgr,1);
        Ydim=size(S.Bgr,2);
    else
        S=SampleCell([Xdim,Ydim]);
    end
    Bgr=single(S.Bgr);
    Bgr=Bgr/S.LaserIntensity*LaserIntensity;
    BgrBW=S.BgrBW;
    clear('S');


    %% Pararmeters.
    fprintf('    ==================================================================\n')
    fprintf('    Spark movie settings:\n')
    fprintf('    %-38s%0.2f um X %0.2f um, %0.1f ms\n','Resolution (x,y,t):',xyt_dim(1),xyt_dim(2),xyt_dim(3));
    fprintf('    %-38s%0.0f%%\n','Laser Intensity:',LaserIntensity);
    fprintf('    %-38s%0.2f\n','Device Gain:',Gain);
    fprintf('    %-38s%0.0f\n','CaR Number:',CaReleaseNum);
    fprintf('    %-38s%0.3f\n','Spark amplitude (dF/F0):',CaReleaseAmp);
    fprintf('    %-38s%0.3f\n','Apparent diffusion K (um^2/s):',ApparentDiffusionK);
    fprintf('    %-38s%0.2f / %0.2f a.u.\n','Detector offset/noise:',DetectorOffset,DetectorGaussNoise);
    if p.ExpAmp
        fprintf('    %-38sExponential, center = %0.2f\n','Amplitude distribution:',CaReleaseAmp);
    elseif p.RampAmp
        fprintf('    %-38sRamp from 0 to %0.2f\n','Amplitude distribution:',CaReleaseAmp);
    else
        fprintf('    %-38sConstant, center = %0.2f\n','Amplitude distribution:',CaReleaseAmp);
    end
    fprintf('    ==================================================================\n')
    
    %% coordinates template.
    fprintf('    %-38s','Generating coordinates:');
    P=generateCaReleaseProbabilityMap(Xdim,1.075/xyt_dim(1),Ydim,1.935/xyt_dim(1));
    [xPos,yPos]=generateCaReleaseProbabilityReleasePos(P,CaReleaseNum);
    CaReleaseNum=numel(xPos);
    tPos=round((randn(CaReleaseNum,1)*CaReleaseDurationSD)/xyt_dim(3)+ReleaseTrigger);
    
    xytTemplate=cat(2,xPos,yPos,tPos);
    
    clear('tPos','xPos','yPos','P')
    fprintf('done\n    %-38s%0.0f\n','CaR Number New:',CaReleaseNum); %fprintf('done\n');
    
    
    %% Amplitude distribution.
    if p.ExpAmp
        AmpDistribution=exprnd(CaReleaseAmp,[CaReleaseNum,1]);
    elseif p.RampAmp
        AmpDistribution=rand(CaReleaseNum,1)*CaReleaseAmp;
    else
        AmpDistribution=ones(CaReleaseNum,1)*CaReleaseAmp;
    end
    
    %% Ca release recording coordinates.
    xyt=round(rand(CaReleaseNum,1)*size(xytTemplate,1));
    for k=1:numel(xyt)
        if xyt(k)==0
            continue;
        end
        if ~BgrBW(xytTemplate(xyt(k),1),xytTemplate(xyt(k),2));
            xyt(k)=0;
        end
    end
    
    while min(xyt)<=0
        bw=xyt<=0;
        xytTemp=round(rand(sum(bw(:)),1)*size(xytTemplate,1));
        xyt(bw)=xytTemp;
        for k=1:numel(xyt)
            if xyt(k)==0; continue; end
            if ~BgrBW(xytTemplate(xyt(k),1),xytTemplate(xyt(k),2)); xyt(k)=0; continue; end
        end
    end
    
    xyt=xytTemplate(xyt,:);
    clear('xytTemplate','bw','k');
    
    %% Ca release recording. Put sparks in.
    
    CICRNoiseFree=zeros(Xdim,Ydim,Tdim,'single');
    SparkPosition=zeros(CaReleaseNum,6);
    fprintf('%-38s      0%%','    Putting in Ca Release:');
    
    CurrSarkSize=[size(SampleCaReleaseSpark,1),size(SampleCaReleaseSpark,2),size(SampleCaReleaseSpark,3)];
    
    for k=1:CaReleaseNum
        % Define x positions to put in.
        x1=round(xyt(k,1)-CurrSarkSize(1)/2);
        x2=x1+CurrSarkSize(1)-1;
        x1Spark=1;
        x2Spark=CurrSarkSize(1);
        if x1<1
            xyshift=1-x1;
            x1=1;
            x1Spark=xyshift+1;
        end
        if x2>Xdim
            xyshift=x2-Xdim;
            x2=Xdim;
            x2Spark=CurrSarkSize(1)-xyshift;
        end
        
        % Define y positions to put in.
        y1=round(xyt(k,2)-CurrSarkSize(2)/2);
        y2=y1+CurrSarkSize(2)-1;
        y1Spark=1;
        y2Spark=CurrSarkSize(2);
        if y1<1
            xyshift=1-y1;
            y1=1;
            y1Spark=xyshift+1;
        end
        if y2>Ydim
            xyshift=y2-Ydim;
            y2=Ydim;
            y2Spark=CurrSarkSize(2)-xyshift;
        end
        
        % Define t positions to put in.
        t1=xyt(k,3);
        t2=t1+CurrSarkSize(3)-1;
        t1Spark=1;
        t2Spark=CurrSarkSize(3);
        if t1<1
            xyshift=1-t1;
            t1=1;
            t1Spark=xyshift+1;
        end
        if t2>Tdim
            xyshift=t2-Tdim;
            t2=Tdim;
            t2Spark=CurrSarkSize(3)-xyshift;
        end

        % Calculate Bgr
        CurrBgr=Bgr(x1:x2,y1:y2);
        CurrBgr=median(CurrBgr(:));
        
        CurSparkTemp=SampleCaReleaseSpark*AmpDistribution(k)*CurrBgr;
        
        
        CICRNoiseFree(x1:x2,y1:y2,t1:t2)=CICRNoiseFree(x1:x2,y1:y2,t1:t2)+...
            CurSparkTemp(x1Spark:x2Spark,y1Spark:y2Spark,t1Spark:t2Spark);

        SparkPosition(k,:)=[k,xyt(k,1),xyt(k,2),xyt(k,3),CurrBgr,...
            AmpDistribution(k)];
        
        fprintf('\b\b\b\b%3.0f%%',floor(k/CaReleaseNum*100));
        
    end

    clear('j','t1','t2','x1','x2','y1','y2','xyt','CurrBgr','CurrSarkSize','SampleSpark','xytTemplate')
    fprintf('\n');
    
    
    %% Check the amplitude.
    if ~isempty(ReleaseAmplitude)
        fprintf('    %-38s','Adjusting amplitude:');
        dFF0=sum(sum(sum(CICRNoiseFree,3),2),1);
        dFF0=dFF0/sum(Bgr(BgrBW));
        
        FF0Ratio=max(dFF0)/ReleaseAmplitude;
        CICRNoiseFree=CICRNoiseFree/FF0Ratio;
        fprintf('done\n');
    else
        FF0Ratio=1;
    end
    
    
    %% Apply Ca diffusion.
    fprintf('%-38s      0%%','    Calcium diffusing:');
    p.PSFdiffusion=generatePSF(xyt_dim,ApparentDiffusionK);
    PSFdiffusion=p.PSFdiffusion;
    p.PSFdiffusionHalf=generatePSF([xyt_dim(1:2),xyt_dim(3)/2],ApparentDiffusionK);
    PSFdiffusionHalf=p.PSFdiffusionHalf;
    CICRBeforeDiffusion=CICRNoiseFree;

    BgrBW=imfill(BgrBW,'holes');
    for k=2:Tdim
        % Diffusion.
        CICRLastTimePieceDiffuse=ImDiffuseParallel(CICRNoiseFree(:,:,k-1),BgrBW,PSFdiffusion,CPU_NumCores);
        CICRCurrTimePieceDiffuse=ImDiffuseParallel(CICRNoiseFree(:,:,k),BgrBW,PSFdiffusionHalf,CPU_NumCores);
        
        CICRFF0Removal=(CICRLastTimePieceDiffuse+CICRCurrTimePieceDiffuse)*DecayFF0Rate*xyt_dim(3)/2;
        
        CICRNoiseFree(:,:,k)=CICRLastTimePieceDiffuse+CICRCurrTimePieceDiffuse-CICRFF0Removal;
        
        fprintf('\b\b\b\b%3.0f%%',floor(k/Tdim*100));
    end
    clear('CICRLastTimePieceDiffuse','CICRCurrTimePieceDiffuse','CICRFF0Removal');
    fprintf('\n');
    
    
    %% Record the Release Ca Pos.
    SparkPositionMovie=zeros(size(CICRNoiseFree),'single');
    for k=1:size(SparkPosition,1)
        SparkPosition(k,6)=SparkPosition(k,6)/FF0Ratio;
        SparkPositionMovie(SparkPosition(k,2),SparkPosition(k,3),SparkPosition(k,4))=...
            SparkPositionMovie(SparkPosition(k,2),SparkPosition(k,3),SparkPosition(k,4))+SparkPosition(k,6);
    end
    
    %% Add Poissonian noise.
    fprintf('    %-38s','Generating noise:');
    SparkOnCellWithPoissonNoise=zeros(size(CICRNoiseFree),'single');
    parfor k=1:size(CICRNoiseFree,3)
        SparkOnCellWithPoissonNoise(:,:,k)=Gain*poissrnd((CICRNoiseFree(:,:,k)+Bgr)/Gain);
    end
    
    % Final movie.
    CICRRecording=SparkOnCellWithPoissonNoise + DetectorOffset+randn(size(SparkOnCellWithPoissonNoise))*DetectorGaussNoise;
    fprintf('done\n')
    
    %% Output
    if nargout==1
        output.CaRelease=CICRBeforeDiffusion;
        output.CaT=CICRRecording;
        output.CaTNoiseFree=CICRNoiseFree;
        output.CaReleasePos=SparkPositionMovie;
        output.CaReleaseInfo=SparkPosition;
        output.CaReleaseInfoFormat='(1)ID (2)dF/F0 (3)Bgr (4)xc (5)yc (6)t_onset';
        output.SampleCaRelease=SampleCaReleaseSpark;
        output.Bgr=Bgr;
        output.BgrBW=BgrBW;
        output.InputParameter=p;
        varargout{1}=output;
    end
end



function Icon=ImDiffuseParallel(I,bw,PSF,CPUCore_Num)
    if size(I,1)>size(I,2); I=I'; end
    step=floor(size(I,2)/CPUCore_Num);
    stepStart=1:step:size(I,2);
    stepEnd=stepStart+step-1;
    stepEnd(end)=size(I,2);
    
    Icon=zeros(size(I));
    parfor k=1:numel(stepStart)
        Icon=Icon+ImDiffuse(I,bw,PSF,1,size(I,1),stepStart(k),stepEnd(k));
    end
    if size(I,1)>size(I,2); Icon=Icon'; end
end

function Icon=ImDiffuse(I,bw,PSF,k1,k2,j1,j2)
    Xdim=size(I,1);
    Ydim=size(I,2);
    Icon=zeros([Xdim,Ydim]);
    
    PSFsiz=size(PSF);
    for k=k1:k2
        for j=j1:j2
            x1=k-(PSFsiz(1)-1)/2;
            x2=k+(PSFsiz(1)-1)/2;
            x1Spark=1;
            x2Spark=PSFsiz(1);
            if x1<1
                xyshift=1-x1;
                x1=1;
                x1Spark=xyshift+1;
            end
            if x2>Xdim
                xyshift=x2-Xdim;
                x2=Xdim;
                x2Spark=PSFsiz(1)-xyshift;
            end
            
            y1=j-(PSFsiz(2)-1)/2;
            y2=j+(PSFsiz(2)-1)/2;
            y1Spark=1;
            y2Spark=PSFsiz(2);
            if y1<1
                xyshift=1-y1;
                y1=1;
                y1Spark=xyshift+1;
            end
            if y2>Ydim
                xyshift=y2-Ydim;
                y2=Ydim;
                y2Spark=PSFsiz(2)-xyshift;
            end
            
            currPSF=PSF(x1Spark:x2Spark,y1Spark:y2Spark).*bw(x1:x2,y1:y2);
            currPSF_sum=sum(currPSF(:));
            if currPSF_sum==0;
                continue;
            else
                currPSF=currPSF/currPSF_sum;
            end
            
            currPSF=currPSF*I(k,j);
            Icon(x1:x2,y1:y2)=Icon(x1:x2,y1:y2)+currPSF;
        end
    end
end


function varargout=SparkDsptSimAllInOneV1_2(varargin)
    %% Spark simulation with decriptive model.
    %% Input.
    p=inputParser;
    % Parameters for the HillExp function.
    p.addParameter('Mu',30,@(x)isscalar(x));                             % In millisecond.
    p.addParameter('Sigma',6,@(x)isscalar(x));                           % In millisecond.
    p.addParameter('Tau',40,@(x)isscalar(x));                            % In millisecond.
    
    % Parameters for the FWHM function.
    p.addParameter('FWHMMax',0.56,@(x)isscalar(x));                       % In micrometer.
    p.addParameter('TauFWHM',15,@(x)isscalar(x));                        % In millisecond.
    
    p.addParameter('SpR',0.28,@(x)isscalar(x));                         % In micrometer.
    p.addParameter('TpR',1,@(x)isscalar(x));                             % In millisecond.
    
    p.addParameter('DecayLimit',0.001,@(x)isscalar(x) && x<=0.01 && x>0);
    parse(p, varargin{:});
    p=p.Results;
    clear('varargin')
    
    %% 2D Gaussian surface on X and Y, and ExpExp on temporal dimension.
    % p(1), amplitude; p(2), mu; p(3), exponential tau; p(4), gaussian sigma.
    GauConvExpFWHM=@(FWHMMax,T,Mu,TauFWHM) FWHMMax*(1-exp((-(T-Mu)./TauFWHM))).*((1-exp((...
        -(T-Mu)./TauFWHM)))>0)+1e-9;
    GauConvExp=@(x,AmpMax,Mu,Tau,Sigma) AmpMax/2.*exp((Sigma^2+2*Tau*(Mu-x))/2/Tau^2).*(1-(Sigma^2+...
        Tau*(Mu-x))./abs(Sigma^2+Tau*(Mu-x)).*erf(abs(Sigma^2+Tau*(Mu-x))/sqrt(2)/Sigma/Tau));
    GauConvExpSpark=@(X,Y,T,FWHMMax,Mu,TauFWHM,AmpMax,Tau,Sigma) exp(-(X.^2+Y.^2)/2./(...
        GauConvExpFWHM(FWHMMax,T,Mu,TauFWHM)).^2).*GauConvExp(T,AmpMax,Mu+Sigma*2,Tau,Sigma);
    
    %% Dimensions.
    tlen=-p.Tau*log(p.DecayLimit)+round(p.Mu);
    T=0:p.TpR:tlen;
    
    PositivePart=0:p.SpR:p.FWHMMax*3;   NegativePart=sort(PositivePart(2:end),'descend')*(-1);
    XYrange=[NegativePart PositivePart];
    
    [X,Y,T]=ndgrid(XYrange,XYrange,T);
    clear('tlen','PositivePart','NegativePart','XYrange')
    
    %% Generate spark.
    spark=GauConvExpSpark(X,Y,T,p.FWHMMax,p.Mu,p.TauFWHM,1,p.Tau,p.Sigma);
    clear('X','Y','T','HillexpampSpark','Amp','FWHM')
    
    %% Normalization.
    spark=spark/max(spark(:));
    %% Remove zeros in the matrix / Shrink the matrix.
    spark_positive_index=find(spark>0);
    [spark_x,spark_y,spark_t]=ind2sub(size(spark),spark_positive_index);
    spark=spark(min(spark_x):max(spark_x),min(spark_y):max(spark_y),min(spark_t):max(spark_t));

    %% Output.
    if nargout==1
        varargout(1) = {spark};
    end
    if nargout==2
        varargout(1) = {spark};
        varargout(2) = {linescan};
    end
end

function I=generateCaReleaseProbabilityMap(Xdim,xspace,Ydim,yspace)
    sizsigma=0.5; % pixel; 0.107um.
    Gau2D=@(x,y,x0,y0,sigma) exp(-((x-x0).^2+(y-y0).^2)/2/sigma^2);
    
    I=zeros(Xdim,Ydim);
    
    [x,y]=ndgrid(1:Xdim,1:Ydim);
    x0=1:xspace:Xdim;
    y0=1:yspace:Ydim;
    
    for k=1:numel(x0)
        for j=1:numel(y0)
            I=I+Gau2D(x,y,x0(k),y0(j),sizsigma);
        end
    end
    I=I/max(I(:));
end

function [xPos,yPos]=generateCaReleaseProbabilityReleasePos(P,NumCaR)
    [x,y]=ndgrid(1:size(P,1),1:size(P,2));
    P=round(P*NumCaR/sum(P(:)));
    
    xPos=[];    yPos=[];
    for k=1:size(P,1)
        for j=1:size(P,2)
            if P(k,j)>0
                currX=ones(P(k,j),1)*x(k,j); xPos=cat(1,xPos,currX);
                currY=ones(P(k,j),1)*y(k,j); yPos=cat(1,yPos,currY);
            end
        end
    end
end


function [PSF,CaDiffuse]=generatePSF(xyt_dim,ApparentDiffusionK)
    CaDiffuse=@(t,D,x,y) exp(-(x.^2+y.^2)/(4*D*t));
    CaDiffuseSize=-30:xyt_dim(1):30;
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
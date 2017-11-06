function CleanObj=CICRcleanSimp(I,Bgr,Mask,xyt_dim,varargin)
    %% CleanObj=CICRcleanSimp(I,Bgr,Mask,xyt_dim,varargin);
    % CICRcleanSimp calculates the Calcium release map, and it is the kernel
    % funciton of the CaCLEAN algorithm.
    % 
    % Inputs:
    %   I: Image matrix containing the upstroke of a calcium transient.
    %      Usually it should be denoised before use.
    %   Bgr: Image matrix containing the basal fluorescence of a calcium
    %      transient. Usually it should be denoised simultaneously with I
    %      and it should contain more than one frame.
    %   Mask: the binary cell mask.
    %   xyt_dim: pixel dimension of the I, Bgr and Mask in spatial and
    %      temporal axes. 
    %
    % Output:
    %   CleanObj: the struct containing the calculated Calcium release map
    %      and related metadata.
    %
    % Qinghai Tian
    % Institute for Molecular Cellbiology
    % Medical Facalty of University of
    % Saarland.
    % Homburg, Germany.
    % tian_qhcn@icloud.com

    %% Input parameters.
    In=inputParser;
    In.addParameter('CaThresholdOverwrite',[], @(x) isempty(x) || (numel(x)==1 && x>0));
    In.addParameter('CaCleanThreshold',3.5, @(x) numel(x)==1 && x>0);
    In.addParameter('PSFAmplitude',[],@(x) isempty(x) || (numel(x)==1 && x>0)); % in photon.
    In.addParameter('ApparentDiffusionK',60,@(x) numel(x)==1 && x>0); % um^2/s.
    In.addParameter('CleanDiffusionK',   30,@(x) numel(x)==1 && x>0); % um^2/s.
    In.addParameter('CaBlur',0.3,@(x) numel(x)==1 && x>0);
    In.addParameter('CleanCycle',2e6,@(x) numel(x)==1 && x>0);
    In.addParameter('ConsiderDiffusion',-1,@(x)isscalar(x) && numel(x)==1);
    
    parse(In, varargin{:});
    In=In.Results;
    
    %% Storage of parameters.
    In.CaCleanVersion='CaCleanMaxPixel_V_1_2_Simp';

    I=single(I);
    if ~isempty(Bgr)
        Bgr=single(Bgr);
    end
    
    CleanObj.InputParameters=In;
    CleanObj.Mask=Mask;
    CleanObj.xyt_dim=xyt_dim;
    if ~isempty(Bgr)
        CleanObj.DataBgr=mean(Bgr,3);
        I=bsxfun(@minus,single(I),single(CleanObj.DataBgr));
    end
    
    % if ~isempty(Ipre)
    %     Ipre=bsxfun(@minus,single(Ipre),single(CleanObj.DataBgr));
    % end
    CleanObj.Data=I;
    
    IsizRaw=size(I);
    
    %% PSF for Ca release processes by CLEAN algorithm.
    [PSF,CaDiffuse]=generatePSF(xyt_dim,In.ApparentDiffusionK);
    cleanPSF=generatePSF(xyt_dim,In.CleanDiffusionK);
    PSFsiz=ceil((size(cleanPSF,1)+1)/2);
  

    %% Ca release CaCleanThreshold.
    if isempty(In.CaThresholdOverwrite)
        Bgr=diff(Bgr,1,3);
        IapparantDiffSigma=std(Bgr(repmat(Mask,[1,1,size(Bgr,3)])));
        CleanObj.CaCleanThreshold=median(IapparantDiffSigma*In.CaCleanThreshold);
        clear('IapparantDiffSigma');
        CaCleanThreshold=CleanObj.CaCleanThreshold;
    else
        CleanObj.CaCleanThreshold=In.CaThresholdOverwrite;
        CaCleanThreshold=CleanObj.CaCleanThreshold;
   end
    clear('Bgr')
    % fprintf('   Threshold = %0.3f\n',CaCleanThreshold);
    %% Define the clean amplitude.
    if isempty(In.PSFAmplitude)
        I_photonNum=bsxfun(@times,diff(I,1,3),Mask);
        I_photonNum=sum(I_photonNum(I_photonNum>CaCleanThreshold));
        I_photonNum=I_photonNum/In.CleanCycle/size(I,3);
        CleanObj.PSFAmplitude=I_photonNum;
        clear('I_photonNum')
    else
        CleanObj.PSFAmplitude=In.PSFAmplitude;
    end
   
    cleanPSF=cleanPSF*CleanObj.PSFAmplitude;
    CleanObj.DiffusionPSF=PSF;
    CleanObj.CLEANPSF=cleanPSF;
    CleanObj.Diffuse=CaDiffuse;
    CleanObj.DiffuseK=In.ApparentDiffusionK; % um2/s
    CleanObj.DiffuseUnit='um^2/s';

    %% Main clean process here.
    % Clear the edge so that the clean kernel will not go out of boundary.
    I=padarray(I,[PSFsiz,PSFsiz,0],'both','replicate');
    Mask=padarray(Mask,[PSFsiz,PSFsiz],'both');
    
    % if isempty(Ipre)
    Ipre=cat(3,zeros(size(I,1),size(I,2)),I(:,:,1:end-1));
    % end
    
    
    CaReleaseCounting=zeros(size(I),'uint32');
    parfor k=1:size(I,3)
        I1=Ipre(:,:,k);
        I2=I(:,:,k);
        I1=imfilter(I1,PSF,'conv','replicate','same');
        I_diff=I2-I1.*Mask;
        I_diff=I_diff-CaCleanThreshold;
        currxy=CaCLEANXYKernel(double(I_diff),Mask,0,double(cleanPSF(:))); %#ok<PFBNS>
        CaReleaseCounting(:,:,k)=currxy;
    end

    CaReleaseCounting=CaReleaseCounting((PSFsiz+1):(PSFsiz+IsizRaw(1)),(PSFsiz+1):(PSFsiz+IsizRaw(2)),:);
    CaRelease_max=max(CaReleaseCounting(:));
    if CaRelease_max<=255
        CaReleaseCounting=uint8(CaReleaseCounting);
    elseif CaRelease_max<=65535
        CaReleaseCounting=uint16(CaReleaseCounting);
    end
    CleanObj.CaReleaseCounting=CaReleaseCounting;
    
    CaRelease2D=single(CaReleaseCounting)*CleanObj.PSFAmplitude;
    CleanObj.CaRelease2D=smoothn(sum(CaRelease2D,3),In.CaBlur);
    
    CaRelease=single(CaReleaseCounting)*CleanObj.PSFAmplitude;
    for k=1:size(CaRelease,3)
        CaRelease(:,:,k)=smoothn(CaRelease(:,:,k),In.CaBlur);
    end
    CleanObj.CaRelease=CaRelease;
    %% Rebuild image stack.
    if In.ConsiderDiffusion>=0
        CaReleaseRebuilt=CICRrebuild(CleanObj,'ConsiderDiffusion',In.ConsiderDiffusion);
        CaReleaseRebuilt=single(CaReleaseRebuilt*CleanObj.PSFAmplitude);
        CleanObj.CaReleaseRebuilt=bsxfun(@plus,CaReleaseRebuilt,CleanObj.DataBgr);
    end
end
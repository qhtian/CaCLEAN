function S=CICRrebuildSimp(S,varargin)
    %% CleanObj=CICRrebuildSimp(CleanObj,varargin);
    % CICRrebuildSimp calculates the upstroke of a calcium transient from the
    % calculated calcium release map that is derived with CaCLEAN algorithm.
    %
    % Inputs:
    %   CleanObj: the struct result from CaCLEAN (CICRcleanSimp) function.
    %   Name-Value Parameters:
    %   DecayFF0Rate: default 0.03. The program follows F = F - dF/F0*K*t.
    %   ConsiderDiffusion: default true.
    %
    % Output:
    %   CleanObj: the struct containing the calculated upstroke of a
    %   calcium transient.
    %
    % Qinghai Tian
    % Institute for Molecular Cellbiology
    % Medical Facalty of University of
    % Saarland.
    % Homburg, Germany.
    % tian_qhcn@icloud.com
    
    %%
    p=inputParser;
    p.addParameter('DecayFF0Rate',0.03,@(x)isscalar(x) && x>=0);              % In F = F - dF/F0 * K * t.
    p.addParameter('ConsiderDiffusion',1,@(x)isscalar(x) && x>=0);              % In F = F - dF/F0 * K * t.
    parse(p, varargin{:});
    p=p.Results;
    
    ConsiderDiffusion=p.ConsiderDiffusion;
    DecayFF0Rate=p.DecayFF0Rate;
    CPU_NumCores=feature('NumCores');

    %% Apply Ca diffusion.
    Snumel=numel(S);
    for j=1:Snumel
        fprintf('     %0.0f / %0.0f%-34s      0%%',j,Snumel,'    Calcium diffusing:');
        % PSFdiffusionHalf=generatePSF([S(j).xyt_dim(1:2),S(j).xyt_dim(3)/2],S(j).CaDiffuseK);
        PSFdiffusionHalf=S(j).CLEANPSF;
        PSFdiffusionHalf=PSFdiffusionHalf/sum(PSFdiffusionHalf(:));
        
        try
            PSFdiffusion=S(j).DiffusionPSF;
        catch
            PSFdiffusion=generatePSF(S.xyt_dim,S.CaDiffuseK);
        end
        
        CICRrebuiltStack=double(S(j).CaReleaseCounting);
    
        Isiz=size(CICRrebuiltStack);
        % Diffusion.
        for k=1:Isiz(3)
            if ConsiderDiffusion
                if k==1
                    CICRLastTimePieceDiffuse=ImDiffuseParallel(zeros(Isiz(1),Isiz(2)),S(j).Mask,PSFdiffusion,CPU_NumCores);
                else
                    CICRLastTimePieceDiffuse=ImDiffuseParallel(CICRrebuiltStack(:,:,k-1),S(j).Mask,PSFdiffusion,CPU_NumCores);
                end
            else
                if k==1
                    CICRLastTimePieceDiffuse=zeros(Isiz(1),Isiz(2));
                else
                    CICRLastTimePieceDiffuse=CICRrebuiltStack(:,:,k-1);
                end
            end
            
            CICRCurrTimePieceDiffuse=ImDiffuseParallel(CICRrebuiltStack(:,:,k),S(j).Mask,PSFdiffusionHalf,CPU_NumCores);
            
            CICRFF0Removal=(CICRLastTimePieceDiffuse+CICRCurrTimePieceDiffuse)*DecayFF0Rate*S(j).xyt_dim(3)/2;
            CICRrebuiltStack(:,:,k)=CICRLastTimePieceDiffuse+CICRCurrTimePieceDiffuse-CICRFF0Removal;
            % CICRrebuiltStack(:,:,k)=CICRrebuiltStack(:,:,k)*S(j).PSFAmplitude;
            fprintf('\b\b\b\b%3.0f%%',floor(k/size(CICRrebuiltStack,3)*100));
        end
        CICRrebuiltStack=CICRrebuiltStack*S(j).PSFAmplitude;
        clear('CICRLastTimePieceDiffuse','CICRCurrTimePieceDiffuse','CICRFF0Removal');
        fprintf('\n');
        
        S(j).CICRrebuilt=CICRrebuiltStack;
    end
end

% function [PSF,CaDiffuse]=generatePSF(xyt_dim,ApparentDiffusionK)
%     CaDiffuse=@(t,D,x,y) exp(-(x.^2+y.^2)/(4*D*t));
%     CaDiffuseSize=-30:xyt_dim(1):30;
%     if mod(numel(CaDiffuseSize),2)==0
%         CaDiffuseSize=numel(CaDiffuseSize)/2;
%         CaDiffuseSize=0:xyt_dim(1):xyt_dim(1)*CaDiffuseSize;
%         CaDiffuseSize=cat(2,-fliplr(CaDiffuseSize(2:end)),CaDiffuseSize);
%     end
%     [x,y]=ndgrid(CaDiffuseSize,CaDiffuseSize);
%     PSF=CaDiffuse(xyt_dim(3)/1000,ApparentDiffusionK,x,y);
%     PSF_bw=PSF>max(PSF(:))*0.001;
%     PSF_bw=sum(PSF_bw,2);
%     PSF_bw=ceil(sum(PSF_bw>0)/2);
%     PSF_siz=(size(PSF,1)+1)/2;
%     PSF=PSF((PSF_siz-PSF_bw):(PSF_siz+PSF_bw),(PSF_siz-PSF_bw):(PSF_siz+PSF_bw));
%     PSF=PSF/sum(PSF(:));
% end


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
    RecrdSize=size(I);
    Icon=zeros(RecrdSize);
    
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
            if x2>RecrdSize(1)
                xyshift=x2-RecrdSize(1);
                x2=RecrdSize(1);
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
            if y2>RecrdSize(2)
                xyshift=y2-RecrdSize(2);
                y2=RecrdSize(2);
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
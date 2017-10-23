function I=sCMOSdynamicOffsetCorrection(I,bw,nthdim)
    % Function sCMOSdynamicOffsetCorrection removes dynamic offset of sCMOS
    % camera.
    %
    %   I = sCMOSdynamicOffsetCorrection(I, bw, nthdim);
    %   I = sCMOSdynamicOffsetCorrection(I, sCMOSOffset);
    %
    %   Inputs:  1) I:      image stack need to correct;
    %            2) bw:     area that will be used to calculate the offset
    %                       of the current frame. Default = true(size(I));
    %            3) nthdim: dimension that the dynamic offset is constant.
    %                       This is optional, default value is 2, dim = 2;
    %            4) sCMOSOffset: [dim, x1, x2, x3, x4] for area x1:x2 and x3:x4;
    %
    % Note:
    %   It will also remove the fixed offset of the camera in the specified
    %   dimension.
    %
    % Developed by Qinghai Tian (tian_qhcn@icloud.com).
    % More info could be found www.lipplab.de.
    %
    
    %%
    if nargin==1
        disp('   Too few inputs. Skip Now and nothing done.');
        return;
    elseif nargin==2
        nthdim=2;
    elseif nargin==3
        if ~((nthdim==1) || (nthdim==2))
            disp('   Dim info is wrong. Skip Now and nothing done.');
            return;
        end
    elseif (nargin<1) || (nargin>3)
        disp('   Wrong number of inputs. Skip Now and nothing done.');
        return;
    end
    
    
    if numel(bw)==5
        sCMOSOffset=bw;
        nthdim=bw(1);
        if (nthdim==1) || (nthdim==2)
        else
            disp('   Dim info from 2nd input is wrong. Skip Now and nothing done.');
            return;
        end
        
        Isiz=size(I);
        bw=false(Isiz(1),Isiz(2));
        x1=round(sCMOSOffset(2)*Isiz(nthdim));
        if x1<1; x1=1; end;        if x1>Isiz(nthdim); x1=Isiz(nthdim); end
        x2=round(sCMOSOffset(3)*Isiz(nthdim));
        if x2<1; x2=1; end;        if x2>Isiz(nthdim); x2=Isiz(nthdim); end
        x3=round(sCMOSOffset(4)*Isiz(nthdim));
        if x3<1; x3=1; end;        if x3>Isiz(nthdim); x3=Isiz(nthdim); end
        x4=round(sCMOSOffset(5)*Isiz(nthdim));
        if x4<1; x4=1; end;        if x4>Isiz(nthdim); x4=Isiz(nthdim); end
        
        if nthdim==1
            bw(x1:x2,:)=true; bw(x3:x4,:)=true;
        elseif nthdim==2
            bw(:,x1:x2)=true; bw(:,x3:x4)=true;
        end
    end
    %%
    imClass=class(I);
    if ~strcmpi(imClass,'double')
        I=single(I);
    end
    
    Ibr=bsxfun(@times,I,bw);
    IbrDim2Raw=sum(Ibr,nthdim);
    clear('Ibr')
    
    bw_sum=sum(bw,nthdim);
    
    IbrDim2=bsxfun(@rdivide,IbrDim2Raw,bw_sum);
    IbrDim2(isinf(IbrDim2))=0;
    
    I=bsxfun(@minus,I,IbrDim2);
end
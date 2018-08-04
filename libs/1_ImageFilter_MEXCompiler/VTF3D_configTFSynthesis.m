function [TFSynthStruct,vxl_fin] = VTF3D_configTFSynthesis(ImgSize,maxitr,vxl_ini)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DFKsMO_KernelSize = [5 5 5]; % Predefined!!!
DFKsMO_pad = DFKsMO_KernelSize - 1;
DFKsFlat_KernelSize = [3 3 3]; % Predefined!!!
DFKsFlat_pad = DFKsFlat_KernelSize - 1;

PadImgSize = ImgSize + (2 * DFKsMO_pad);

IDXimgDFKsMO = struct('idxs',uint32(zeros(prod([5 5 5]),maxitr)));
IDXimgDFKsFlat = struct('idxs',uint32(zeros(prod([3 3 3]),maxitr)));

itr = 0;

for vxl = vxl_ini : 1 : vxl_ini + (maxitr - 1)
    
    if vxl < prod(ImgSize) + 1
        
        [Rvxl,Cvxl,Pvxl] = ind2sub(ImgSize,vxl);
        
        itr = itr + 1;
        
        [xDFKsMOBlk,yDFKsMOBlk,zDFKsMOBlk] = ndgrid(Rvxl : Rvxl+DFKsMO_pad(1) , Cvxl : Cvxl+DFKsMO_pad(2) , Pvxl : Pvxl+DFKsMO_pad(3));
        IDXsDFKsMO_temp = uint32( xDFKsMOBlk(:) + (yDFKsMOBlk(:)-1)*PadImgSize(1) + (zDFKsMOBlk(:)-1)*PadImgSize(1)*PadImgSize(2) );
        
        [xDFKsFlatBlk,yDFKsFlatBlk,zDFKsFlatBlk] = ndgrid(Rvxl : Rvxl+DFKsFlat_pad(1) , Cvxl : Cvxl+DFKsFlat_pad(2) , Pvxl : Pvxl+DFKsFlat_pad(3));
        IDXsDFKsFlat_temp = uint32( xDFKsFlatBlk(:) + (yDFKsFlatBlk(:)-1)*PadImgSize(1) + (zDFKsFlatBlk(:)-1)*PadImgSize(1)*PadImgSize(2) );
        
        
        pos = mod(vxl,maxitr); pos(pos == 0) = maxitr;
        
        IDXimgDFKsMO(1,1).idxs(:,pos) = IDXsDFKsMO_temp(:);
        IDXimgDFKsFlat(1,1).idxs(:,pos) = IDXsDFKsFlat_temp(:);
    
    else
        break;
    end
    
end

vxl_fin = vxl_ini + itr;

TFSynthStruct.IDX4DFKsMO = IDXimgDFKsMO;
TFSynthStruct.IDX4DFKsFlat = IDXimgDFKsFlat;
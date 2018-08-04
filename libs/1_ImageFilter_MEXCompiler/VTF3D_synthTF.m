function TFStruct = VTF3D_synthTF(TFStruct,TFSynthStruct,FiltResp,FRidxDFK,DFKsMO,DFKsFlat,vxl_ini,bbFlag)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(FRidxDFK) % IF it is NOT the BOUNDARIES & BACKGROUND Component

        % retrieve the last patch
        lastptc = find(TFSynthStruct.IDX4DFKsMO(1,1).idxs(1,:) == 0,2,'first');
        if ~isempty(lastptc)
            lastptc = lastptc(1)-1;
        else
            lastptc = size(TFSynthStruct.IDX4DFKsMO(1,1).idxs,2);
        end

        % Sliding Patch-Sweep
        for ptc = 1 : lastptc % for each voxel apply the kernel patch (= ptc)
            
            WAlin = FiltResp(ptc + vxl_ini) .* DFKsMO(FRidxDFK(ptc),1).k(:) .* TFStruct.DW(:);
            
            %%% Weights Accumulator
            TFStruct.WA(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.WA(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin;
            
            %%% Overlapping and Adding Adjacent Contributions for the
            %%% Synthetic Structure Tensor Field Reconstruction
            TFStruct.T11(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.T11(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).EIGLE.T11(:);
            
            TFStruct.T12(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.T12(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).EIGLE.T12(:);
            
            TFStruct.T13(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.T13(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).EIGLE.T13(:);
            
            TFStruct.T22(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.T22(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).EIGLE.T22(:);
            
            TFStruct.T23(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.T23(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).EIGLE.T23(:);
            
            TFStruct.T33(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.T33(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).EIGLE.T33(:);
            
            %%% Anisotropic Ratio
            TFStruct.ARelong(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.ARelong(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).FAblob(:);
        
            TFStruct.ARcross(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) = TFStruct.ARcross(TFSynthStruct.IDX4DFKsMO(1,1).idxs(:,ptc)) + WAlin .* DFKsMO(FRidxDFK(ptc),1).FAcros(:);
            
        end
    
    
else % IF it is the BOUNDARY and BACKGROUND Component
    
    % Integral WEIGHTs for BOUNDARY and BACKGROUND Components in the
    % Structure Tensor Field
    wBB = 0.5;
        
    % retrieve the last patch
    lastptc = find(TFSynthStruct.IDX4DFKsFlat(1,1).idxs(1,:) == 0,2,'first');
    if ~isempty(lastptc)
        lastptc = lastptc(1)-1;
    else
        lastptc = size(TFSynthStruct.IDX4DFKsFlat(1,1).idxs,2);
    end
    
    for ptc = 1 : lastptc % for each voxel apply the kernel patch (= ptc)
        
        TFStruct.WA(TFSynthStruct.IDX4DFKsFlat(1,1).idxs(:,ptc)) = TFStruct.WA(TFSynthStruct.IDX4DFKsFlat(1,1).idxs(:,ptc)) + ( wBB * FiltResp(ptc) .* DFKsFlat(bbFlag,1).k(:) );
        
    end
end
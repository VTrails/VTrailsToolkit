function pairs = VTR3D_getPairs2BeEvaluated(nElm)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pairs = zeros(((nElm.^2)-nElm)/2,2);
ini_pos = 1;

for iElm = 1 : nElm-1
    pairs(ini_pos: ini_pos + (nElm - iElm) - 1 , : ) = [ iElm * ones( (nElm - iElm) ,1) , (iElm + 1 : nElm)' ];
    ini_pos = ini_pos + (nElm - iElm);
end
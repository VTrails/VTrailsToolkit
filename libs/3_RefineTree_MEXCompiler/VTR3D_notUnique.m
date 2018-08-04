function [RepeatedValues,a_check] = VTR3D_notUnique(a,b)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a MUST be mx1 and b MUST be nx1

if ~isempty(a) && ~isempty(b)
    % Initialization
    a_eval = a(:);
    b_eval = b(:);
    
    a_check = false(size(a_eval));
    
    EvalVector = [ a_eval ; b_eval ];
    
    [EvalVectorSorted,~] = sort(EvalVector);
    
    Ridx = find( diff( EvalVectorSorted ) == 0);
    
    if ~isempty( Ridx )
        RepeatedValues = EvalVectorSorted( Ridx );
        for j = 1 : length( RepeatedValues )
            a_check( a_eval == RepeatedValues(j) ) = true;
        end
    else
        RepeatedValues = [];
    end
    
elseif ~isempty(a) && isempty(b)
    RepeatedValues = [];
    a_check = false(size(a(:)));
else
    RepeatedValues = [];
    a_check =logical([]);
end
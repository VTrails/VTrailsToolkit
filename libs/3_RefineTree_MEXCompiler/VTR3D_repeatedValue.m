function [RepeatedValues,a_check,b_check] = VTR3D_repeatedValue(a,b)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B.: Inputs 'a' and 'b' MUST be matrices mx3 and nx3.

RepeatedValues = NaN(size(a,1)*size(b,1),3);
a_check = false([size(a,1),1]);
b_check = false([size(b,1),1]);
entry = 0;
for ia = 1 : size(a,1)
    for ib = 1 : size(b,1)
        if sum( a(ia,:) == b(ib,:) ) == 3
           entry = entry + 1;
           RepeatedValues(entry,:) = a(ia,:);
           a_check(ia) = true;
           b_check(ib) = true;
        end
    end
end
RepeatedValues = RepeatedValues(~isnan(RepeatedValues(:,1)),:);
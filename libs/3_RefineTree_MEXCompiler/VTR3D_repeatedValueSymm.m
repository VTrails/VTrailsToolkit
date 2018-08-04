function [RepeatedValues,a_check,b_check] = VTR3D_repeatedValueSymm(a,b)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(size(a,2) == size(b,2),'size(a,2) NOT EQUAL TO size(b,2)');

if isempty(a) && isempty(b)

    RepeatedValues = [];
    a_check = logical([]);
    b_check = logical([]);
    
elseif isempty(a) && ~isempty(b)
    
    RepeatedValues = [];
    a_check = logical([]);
    b_check = false([size(b,1),1]);
    
elseif ~isempty(a) && isempty(b)
    
    RepeatedValues = [];
    a_check = false([size(a,1),1]);
    b_check = logical([]);
    
else % ~isempty(a) && ~isempty(b)
    
    RepeatedValues = NaN(size(a,1)*size(b,1),size(a,2));
    
    a_in = a;
    b_in = b;
    
    a_check_temp = false(size(a_in,1),1);
    b_check_temp = false(size(b_in,1),1);
    entry = 0;
    dim2 = size(a_in,2);

    for ai = 1 : size(a_in,1)
        for bi = 1 : size(b_in,1)
            if  sum(a_in(ai,:) == b_in(bi,:)) == dim2 || sum(a_in(ai,:) == fliplr(b_in(bi,:))) == dim2
                entry = entry + 1;
                RepeatedValues(entry,:) = a_in(ai,:);
                a_check_temp(ai) = true;
                b_check_temp(bi) = true;
            end
        end
    end

    a_check = a_check_temp;
    b_check = b_check_temp;
    RepeatedValues = RepeatedValues(~isnan(RepeatedValues(:,1)),:);
end
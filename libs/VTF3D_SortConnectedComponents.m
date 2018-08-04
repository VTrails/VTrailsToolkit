function [sortedOS,ConnectedSegments,BranchPoints] = VTF3D_SortConnectedComponents(OS,dispalyProgress)%
% Use this function to sort and remove redundancy in the Organised Seeds,
% so that they are presented as Sorted and Connected Components (Connected
% Paths).
%
% INPUTS:
%
% - OS: Logical Volume of Organised Seeds, aligned to Mid-Vessel Maxima.
% - dispalyProgress: Bool for verbose.
%
% OUTPUTS: 
%
% - sortedOS: Label MAP (that can be converted in logical Volume) of the
%             Sorted Organised Seeds
% - ConnectedSegments: Struct containing the Size of OS and the respective
%                      sorted list of indexed points for each connected
%                      component (read as Connected Path).
%
% Example:
%
% [sortedOS,ConnectedSegments] = VTF3D_SortConnectedComponents(OS);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, Sept-2017, stefano.moriconi.15@ucl.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dispalyProgress
    tic;
    disp('***********************************************');
    disp('   VTrails: Sort Seeds Connected Components    ');
end

OS = VTF3D_adjustSeedsForImageBoundary(OS);

CC = bwconncomp(OS, 26 );

sortedOS = zeros(size(OS));

ConnectedSegments.Size = CC.ImageSize;
ConnectedSegments.List = [];

BranchPoints = [];

entry = 0;

if dispalyProgress
    ProgressValue = 0;
    ProgressValue_step = 100.0/CC.NumObjects;
    ProgressMessage = 'Sorting Connected Components:';
    ProgressOldString = VTF3D_showProgressCM(ProgressMessage,'',ProgressValue,false);
end

for cc = 1 : CC.NumObjects

    SeedsList = CC.PixelIdxList{cc};
        
    if length(SeedsList) > 1
        [SortedPath,MissingPartsFlag,RemainingSeeds] = VTF3D_convertBinary3DConnCompToSortedPath(SeedsList,CC.ImageSize);
    else
        SortedPath = SeedsList;
        RemainingSeeds = [];
        MissingPartsFlag = false;
    end
    
    if MissingPartsFlag
        
        [SortedPaths,Furcations] = VTF3D_splitMultipleBranches(SortedPath,RemainingSeeds,CC.ImageSize);
        
        for b = 1 : length(SortedPaths)
            sortedOS(SortedPaths(b).Sequence) = entry + 1;
            entry = entry + 1;
            ConnectedSegments.List(entry,1).Sequence = SortedPaths(b).Sequence;
        end
        
        BranchPoints = unique(cat(1,BranchPoints,Furcations));
        
        %%%
        
        tempOS = false(size(OS));
        tempOS(RemainingSeeds) = true;
        
        [sortedtempOS,tempConnectedSegments,tempBranchPoints] = VTF3D_SortConnectedComponents(tempOS,false);
        
        sortedOS(sortedtempOS > 0) = entry + sortedtempOS(sortedtempOS>0);
        ConnectedSegments.List = cat(1,ConnectedSegments.List,tempConnectedSegments.List);
        entry = length(ConnectedSegments.List);
        
        BranchPoints = unique(cat(1,BranchPoints,tempBranchPoints));
        
    else
        
        sortedOS(SortedPath) = entry + 1;
        entry = entry + 1;
        ConnectedSegments.List(entry,1).Sequence = SortedPath;
        
    end
    
    if dispalyProgress
        % Progress Percentage
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTF3D_showProgressCM([],ProgressOldString,ProgressValue,cc == CC.NumObjects);
    end
    
end

sortedOS = single(sortedOS);

if dispalyProgress
    toc;
    disp('VTrails: Sort Connected Components - COMPLETE');
end

function [SortedPath,MissingPartsFlag,RemainingSeeds] = VTF3D_convertBinary3DConnCompToSortedPath(RandomList,Size)
%
% N.B.: The input RandomList MUST NOT have any furcation!
% RandomList can be a binary connected component which MUST represent a
% segment without furcations. Cycles are allowed.
BWimg = false(Size);
BWimg(RandomList) = true;
WeightMask = padarray(BWimg,ones(1,length(Size)),false);
ConnMap = zeros(Size + 2);

Conn6Kernel  = logical(cat(3,[0 0 0;0 1 0;0 0 0],[0 1 0;1 0 1;0 1 0],[0 0 0;0 1 0;0 0 0]));
Conn18Kernel = logical(cat(3,[0 1 0;1 1 1;0 1 0],[1 1 1;1 0 1;1 1 1],[0 1 0;1 1 1;0 1 0]));
Conn26Kernel = logical(cat(3,[1 1 1;1 1 1;1 1 1],[1 1 1;1 0 1;1 1 1],[1 1 1;1 1 1;1 1 1]));

idx = find(WeightMask);

for jj = 1 : length(idx)
    
    [row,col,pln] = ind2sub(Size + 2,idx(jj));
    ConnMap(row,col,pln) = sum( sum( sum( WeightMask(row-1:row+1,col-1:col+1,pln-1:pln+1) &  Conn6Kernel ) ) ) + ...
                           sum( sum( sum( WeightMask(row-1:row+1,col-1:col+1,pln-1:pln+1) & Conn18Kernel ) ) ) + ...
                           sum( sum( sum( WeightMask(row-1:row+1,col-1:col+1,pln-1:pln+1) & Conn26Kernel ) ) ) ;
    clear row col pln;
    
end

ConnMap = 1./ConnMap(2:end-1,2:end-1,2:end-1);

%%% Search 2-sided Minimal Path %%%

D06 = sqrt(3);
D18 = sqrt(2);
D26 = 1;
D00 = 0;

neigh27PrefEuclidDist = cat(3,[D26 D18 D26;D18 D06 D18;D26 D18 D26] , ...
                              [D18 D06 D18;D06 D00 D06;D18 D06 D18] , ...
                              [D26 D18 D26;D18 D06 D18;D26 D18 D26] );
                          
neigh27PrefEuclidDist = neigh27PrefEuclidDist(:);

[~,startIdx] = max(ConnMap(RandomList));
startPnt = RandomList(startIdx); clear startIdx;

minPath = cell(2,1);

for side = 1 : 2
    
    itr = 0;
    ContinueSearch = true;
    minPath{side} = startPnt;
    
    iniPnt = startPnt;
    
    while ContinueSearch
        
        itr = itr + 1;
        
        [n27n,n27nMSK] = VTF3D_get27Neighbours(iniPnt,Size,false);
         n27n = n27n(n27nMSK);
        
        if isempty( find( ConnMap( n27n ) ~= Inf, 1 ) )
            
            ConnMap(iniPnt) = Inf; % Set Isolated Pnt to Inf
            ContinueSearch = false;
            
        else
            [~,MinIdx] = min( ConnMap(n27n).*neigh27PrefEuclidDist(n27nMSK) );
            nxtPnt = n27n(MinIdx);
            minPath{side} = cat(1,minPath{side},nxtPnt);
            ConnMap = VTF3D_avoidGradientDescentLoopBack(iniPnt,nxtPnt,ConnMap);
            iniPnt = nxtPnt;
            
            ContinueSearch = true;
            clear MinIdx nxtPnt;
        end
        
        clear n26n*;
        
        if itr > 1000
            disp(['Too many iterations? Maybe While Loop? - Iterations: ',num2str(itr,'%d')]);
        end
        
    end
    
end

SortedPath = [ flipud(minPath{1}) ; minPath{2}(2:end) ];

MissingPartsFlag = false;

RemainingSeeds = [];

if ~isempty(find(~isinf(ConnMap), 1))
    % disp('WARNING: The Binary Connected Component has at least one furcation. Furcations will be ignored -> missing part in the Binary Connected Component!');
    MissingPartsFlag = true;
    
    RemainingSeeds = find(~isinf(ConnMap));
    
end

function [neigh27,mask] = VTF3D_get27Neighbours(pnt,Size,uniqueness)
% NB: if uniqueness == false, the output will always be a 27 element vector
% with possible NaNs. The ORDER of the Elements, in this case IS IMPORTANT!
% NB: if uniqueness == true, the output neigh27 may vary in length!
% Moreover it will be sorted in ascendig fashion. (ORDER no more important)

neigh27 = [];
mask = [];

for j = 1 : length(pnt)
    
    [PTNx,PTNy,PTNz] = ind2sub(Size,pnt(j));
    
    iniPTNx = PTNx-1;
    finPTNx = PTNx+1;
    
    iniPTNy = PTNy-1;
    finPTNy = PTNy+1;

    iniPTNz = PTNz-1;
    finPTNz = PTNz+1;
    
    [xG,yG,zG] = ndgrid(iniPTNx:finPTNx,iniPTNy:finPTNy,iniPTNz:finPTNz);
    
    ValidityChk = ~ (reshape( (xG<1 | xG>Size(1)) | (yG<1 | yG>Size(2)) | (zG<1 | zG>Size(3)) , [numel(xG),1]));
    
    % Fixing Boundaries
    xG(xG < 1) = 1;
    xG(xG > Size(1)) = Size(1);
    yG(yG < 1) = 1;
    yG(yG > Size(2)) = Size(2);
    zG(zG < 1) = 1;
    zG(zG > Size(3)) = Size(3);
    
    neigh27_temp = reshape(sub2ind(Size,xG,yG,zG),[numel(xG),1]);
    
    % Setting to NaN those invalid
    neigh27_temp(~ValidityChk) = NaN;
   
    % Identify the position of the Selected point
    SelPntChk = neigh27_temp == pnt(j);
    
    % Extracting only those different from the selected point
    neigh27_temp(SelPntChk) = NaN;
    mask_temp = ~isnan(neigh27_temp);
    
    neigh27 = cat( 1, neigh27 , neigh27_temp );
    mask = logical([mask ; mask_temp]);
    
    clear PTN* ini* fin* xG yG zG neigh26_temp mask_temp ValidityChk SelPntChk;
end

if uniqueness
    neigh27 = unique(neigh27(mask));
    mask = [];
end

function U = VTF3D_avoidGradientDescentLoopBack(iniPnt,endPnt,U)

Size = size(U);

[iniPntn26n,~] = VTF3D_get27Neighbours( iniPnt , Size , true );
[endPntn26n,~] = VTF3D_get27Neighbours( endPnt , Size , true );

[Intersection,~] = VTF3D_notUnique(endPntn26n,iniPntn26n);
Intersection = [Intersection;iniPnt];

U(Intersection) = max(U(:));

function [RepeatedValues,a_check] = VTF3D_notUnique(a,b)%
% Use this function to obtain a list of repeated values between two vectors
% a and b different or equal length.
% RepeatedValues will be a column vector listing the repeated values
% between 'a' and 'b'.
% a_check will be a column vector, same size of 'a', with logical true in
% correspondence of the repeated values.
% N.B.: Inputs 'a' and 'b' MUST be vectors.

assert( isempty(a) | isempty(b) | (isvector(a) && isvector(b)),'Inputs of IFM2D_notUnique MUST be vectors of double/integer' );

% Enforcing Column Vectors
a = a(:);
b = b(:);

% Initialization
RepeatedValues = [];
a_check = false(size(a));

EvalVector = [ a ; b ];

[EvalVectorSorted,~] = sort(EvalVector);

Ridx = find( diff( EvalVectorSorted ) == 0);
if ~isempty( Ridx )
    RepeatedValues = EvalVectorSorted( Ridx );
    
    for j = 1 : length( RepeatedValues )
        a_check( a == RepeatedValues(j) ) = true;
    end
end

function OS_out = VTF3D_adjustSeedsForImageBoundary(OS)

assert(~isequal(size(OS),[1 1 1]),'The Image Volume is INVALID!');

frame = VTF3D_getImageFrame3D(size(OS));

OSframe = OS & frame;

OS_out = OS & ~frame;

idx = find(OSframe, 1);

if ~isempty(idx)
   
    for jj = 1 : length(idx)
       
        [Xpos,Ypos,Zpos] = ind2sub(size(OS),idx(jj));
        
        % Correct Xpos
        if Xpos == 1
            Xpos = Xpos + 1;
        end
        
        if Xpos == size(OS,1)
            Xpos = size(OS,1) - 1;
        end
        
        % Correct Ypos
        if Ypos == 1
            Ypos = Ypos + 1;
        end
        
        if Ypos == size(OS,2)
            Ypos = size(OS,2) - 1;
        end
        
        % Correct Zpos
        if Zpos == 1
            Zpos = Zpos + 1;
        end
        
        if Zpos == size(OS,3)
            Zpos = size(OS,3) - 1;
        end
        
        OS_out(Xpos,Ypos,Zpos) = true;
        
    end
    
end


function frame = VTF3D_getImageFrame3D(Size)

frame = false(Size);

frame(1,:,:) = true;
frame(end,:,:) = true;
frame(:,1,:) = true;
frame(:,end,:) = true;
frame(:,:,1) = true;
frame(:,:,end) = true;

function OldValueStr = VTF3D_showProgressCM(ProgressMessage,OldValueStr,Value,Finished)

if ~isempty(ProgressMessage)
    fprintf(ProgressMessage);
end

if ~Finished
    NewValueStr = sprintf(' %3.2f', Value);
    fprintf([OldValueStr, NewValueStr]);
    OldValueStr = repmat(sprintf('\b'), 1, length(NewValueStr));
else
    NewValueStr = sprintf(' %3.2f', Value);
    fprintf([OldValueStr, NewValueStr,' -- DONE\n']);
    OldValueStr = '';
end

function [SortedPaths,Furcations] = VTF3D_splitMultipleBranches(SortedPath,RemainingSeeds,ImageSize)

[xSP,ySP,zSP] = ind2sub(ImageSize,SortedPath); SP = [xSP,ySP,zSP];
[xRS,yRS,zRS] = ind2sub(ImageSize,RemainingSeeds); RS = [xRS,yRS,zRS];

[ICP_SP,~,~,~] = VTF3D_computeICPfor3DPointsSequences(SP,RS);

closenessTHR = 1+sqrt(3);

FurcationPos = [ 1 ; find( diff( ICP_SP <= closenessTHR ) == 1 ) + 1 ; length(SortedPath)];
SortedPaths = struct([]);

sgmt = 0;
Furcations = [];

segmentLengthTHR = 3;

for ff = 1 : length(FurcationPos)-1
    
    Sequence = SortedPath(FurcationPos(ff):FurcationPos(ff+1)-1);
    if length(Sequence) > segmentLengthTHR
        sgmt = sgmt + 1;
        SortedPaths(sgmt).Sequence = Sequence;
        if FurcationPos(ff+1) ~= length(SortedPath)
            Furcations = cat(1,Furcations,SortedPath(FurcationPos(ff+1))); 
        end
    end
end

function [ICPa,ICPb,ICPaIDX,ICPbIDX] = VTF3D_computeICPfor3DPointsSequences(a,b)%
% N.B.: Inputs 'a' and 'b' MUST be matrices mx3 and nx3.

aBOX = repmat(reshape(a',[1,3,size(a,1)]),[size(b,1),1,1]);
bBOX = repmat(b,[1,1,size(a,1)]);

EuclidDist = squeeze( sqrt( sum( (aBOX - bBOX).^2 , 2) ) );
% EuclidDist must be a nxm matrix with the respective eucludean distances
% between b and a.

% Therefore the ICP for b is given by ICPb
[ICPb,ICPbIDX] = min(EuclidDist,[],2);

% And the ICP for a is given by ICPa
[ICPa,ICPaIDX] = min(EuclidDist,[],1);
ICPa = ICPa'; ICPaIDX = ICPaIDX';
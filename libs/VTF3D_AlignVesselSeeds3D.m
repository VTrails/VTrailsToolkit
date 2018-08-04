function OS = VTF3D_AlignVesselSeeds3D(CVM,US,quantileTHR,displayProgress)%
% Use This function to Align the Un-Organised Seeds to the Vessels Mid-Line
% by means of a constrained Gradient Descent.
%
% INPUTS:
%
% - CVM: Scalar Connected Vesselness Map obtained with SLoGS
% - US: Logical Volume of Un-Organised Seeds TO BE Aligned to the vessels
% - quantileTHR: Scalar value as Quantile Threshold within [0,1], being
%                0 -> all Seeds (no intensity Threshold Applied)
%                1 -> only the Seed(s) with MAX Intensity.
% - dispalyProgress: Bool for verbose.
%
% OUTPUTS:
%
% - OS: Logical Volume of Organised Seeds, aligned to Mid-Vessel Maxima.
%
% Example:
%
%  OS = VTF3D_AlignVesselSeeds3D( MaxMAP.CVMsmth, MaxMAP.US, 0.75 );
%
% NB. Here quantileTHR = 0.75, meaning that Seeds with Intensities in the
%     Upper Quartile of CVM and above will be eventually considered.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if displayProgress
    tic;
    disp('***********************************************');
    disp('         VTrails: Aligning Vessel Seeds        ');
end

% Compute Gradient and Element for Hessian (just ONCE) - TIME CONSUMING
VG = VTF3D_computeGradients2ndOrd_mex(CVM);

%msk1 = false(size(US)); %msk1(1:2:end,1:2:end,1:2:end) = true; %US1 = US & msk1;
msk1 = false(size(US));
msk1(1:2:end,1:2:end,1:2:end) = true;
msk1(2:2:end,2:2:end,1:2:end) = true;
msk1(:,:,2:2:end) = true;
msk1(1:2:end,1:2:end,2:2:end) = false;
msk1(2:2:end,2:2:end,2:2:end) = false;
US1 = US & msk1;

%msk2 = false(size(US)); %msk2(2:2:end,2:2:end,2:2:end) = true; %US2 = US & msk2;
msk2 = ~msk1;
US2 = US & msk2;

if displayProgress
    OS1 = AVS3D_ConstrainedGradientDescent3D(US1,VG,1);
    OS2 = AVS3D_ConstrainedGradientDescent3D(US2,VG,2);
else
    OS1 = AVS3D_ConstrainedGradientDescent3D(US1,VG,0);
    OS2 = AVS3D_ConstrainedGradientDescent3D(US2,VG,0);
end

OS = OS1 | OS2;

BinSeedsIdx = find(OS);
if ~isempty(quantileTHR)
    OS( BinSeedsIdx( CVM(BinSeedsIdx) < quantile(CVM(BinSeedsIdx),quantileTHR) ) ) = false;
else
    % Automatic Estimation of the THR as lower Whisker (as in boxplot to highlight outliers!)
    Q3 = quantile(CVM(BinSeedsIdx),0.75);
    Q1 = quantile(CVM(BinSeedsIdx),0.25);
    IQR = Q3-Q1;
    THR = Q1-1.5*IQR; THR(THR<0) = 0;
    
    OS( BinSeedsIdx( CVM(BinSeedsIdx) < THR ) ) = false;
end

if displayProgress
    toc;
    disp('VTrails: Aligning Vessel Seeds - COMPLETE');
end


function OS = AVS3D_ConstrainedGradientDescent3D(US,VG,phaseLabel)

% Creating the History Path (HP) Matrix:
% HP will be seeds x whileITER
% Each row will keep the history of the path of the specific seeds
% Positive indices mean an actual index of the seed position in the IMAGE
% Negative indices mean a loop, merge or cross had happened, so the abs(NegativeIndex) represent the linear index of HP where the paths join.
% Negative indices will not continue with iterations as a repetition is occurred. So those might be useful to trace back the path, or better
% identify the final seed to consider.
% In HP, NegativeIndices will be repeated as update to the matrix.
% Using this algoritm, HP is supposed to reduce iteratively the PositiveIndices every cycle, till convergence. 
% Convergence is obtained when all the elements of the last column of HP are NEGATIVE.

% Initialization
HP = find(US);
Flag = zeros(size(HP));

if phaseLabel > 0
    ProgressValue = 0;
    ProgressMessage = ['Aligning Vessel Seeds - Phase ',num2str(phaseLabel,'%d'),':'];
    ProgressOldString = AVS3D_showProgressCM(ProgressMessage,'',ProgressValue,false);
end

while ~isempty(find( HP(:,end) > 0, 1 ))
    
    IDX = HP( HP(:,end) > 0 , end);
    
    % Compute the EigenVectors of the Seeds to be CENTERED (must be aligned to a stable position)
    [~,VEvec2NS,VEvec3NS,~] = AVS3D_ComputeSelEVecs2(VG,IDX,true,false);
    %disp('Evects: DONE');
    
    % Estimate the Displacement for each Seed towards Stability
    DispAbsN = AVS3D_EstimateDisplacement(VG,VEvec2NS,VEvec3NS,IDX);
    %disp('Displacement: DONE');
    
    % Update Skeleton Seeds Position using the Estimated Displacement
    IDXnew = AVS3D_UpdateSkelSeeds(size(US),IDX,DispAbsN);
    %disp('Update Positions: DONE');
    
    % Updating the History Path Matrix
    [HP,Flag] = AVS3D_UpdateHistoryPath(HP,Flag,IDXnew);
    %disp('Update History: DONE');
    
    clear VEvec* IDX DispAbsN IDXnew;
    
    if phaseLabel > 0
        % Progress Percentage
        ProgressValue = 100*sum(HP(:,end) <= 0)/size(HP,1);
        ProgressOldString = AVS3D_showProgressCM([],ProgressOldString,ProgressValue, sum(HP(:,end) <= 0) == size(HP,1));
    end
    
end

% Need to prune the loops and cross!!!

Flag = sum(Flag,2);
% if ~isequal(unique(Flag),[1 3 5]')
%     warning('There MUST be a bug in the CODE relative to the Path of the Seeds!');
% end

OS = AVS3D_ExtractCenteredSeeds(HP,Flag,size(US));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  DispAbsN = AVS3D_EstimateDisplacement(VG,VEvec2NS,VEvec3NS,idx)%
% Project the Gradient to the EigenVector Bases
% Consider only the Orthogonal Bases to the main direction of the Vessel!

%                                                [.  .  . ]        [  .    .    .  ]
%                                              [.  .  . ]        [  .    .    .  ]
% Gradient Vector (GV): 2 x 3 x points       [Vx Vy Vz]   .*   [ e2x  e2y  e2z ]
% Bases Matrix (BM): 2 x 3 x points          [Vx Vy Vz]        [ e3x  e3y  e3z ]
% Projection Vector (PV): 2 x 1 x points                   

GV = repmat([reshape(VG.IxN(idx),[1,1,length(idx)]), reshape(VG.IyN(idx),[1,1,length(idx)]), reshape(VG.IzN(idx),[1,1,length(idx)])],[2,1,1]);
BM = [reshape(VEvec2NS',[1,3,length(idx)]) ; reshape(VEvec3NS',[1,3,length(idx)])];

PV = sum( GV .* BM , 2 );

% Create a unit displacement vector which belongs to the plane orthogonal
% to the main direction of the vessel (it has its own orientation)

DispRel = sum( repmat(PV,[1,3,1]) .* BM , 1);

DispRelN = DispRel ./ repmat((sum(DispRel.^2 , 2)).^(1/2) , [1 3 1] ); % Normalized by its norm

% Rearranging the Displacement vector as matrix: points x 3
DispRelN = squeeze(DispRelN)';

if isequal(size(DispRelN), [3,1])
    DispRelN = DispRelN';
end

% Map the resulting displacement (in relative reference system) to absolute
% rasterization reference system

RasterFactor = sqrt(3)/2; % Diagonal of Cube

% Estimate the displacement in terms of absolute voxel displacement
% (starting from the considered voxel) 

DispAbsN = sign( round( DispRelN ./ RasterFactor ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  IDXnew = AVS3D_UpdateSkelSeeds(Size,idx,DispAbsN)

[Sx,Sy,Sz] = ind2sub(Size,idx);

SxNEW = Sx + DispAbsN(:,2); %%% !!! PAY ATTENTION TO ROW & COLS (x and y) !!! %%%
SyNEW = Sy + DispAbsN(:,1); %%% !!! PAY ATTENTION TO ROW & COLS (x and y) !!! %%%
SzNEW = Sz + DispAbsN(:,3);

SxNEW(SxNEW < 1) = 1;
SyNEW(SyNEW < 1) = 1;
SzNEW(SzNEW < 1) = 1;

SxNEW(SxNEW > Size(1)) = Size(1);
SyNEW(SyNEW > Size(2)) = Size(2);
SzNEW(SzNEW > Size(3)) = Size(3);

IDXnew = sub2ind(Size,SxNEW,SyNEW,SzNEW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HistoryPath,Flag] = AVS3D_UpdateHistoryPath(HistoryPath,Flag,idxNEW)

%%% Branch1
% NB: It is Necessary to transfer the idxNEW to the idxNEWtemp (have a
% mixture of positive/negative indices) in order to keep accordingly the
% right number of the pointers.

idxNEWtemp = HistoryPath(:,end);
idxNEWtemp( HistoryPath(:,end) > 0 ) = idxNEW;

FlagNEWtemp = zeros(size(idxNEWtemp));

%%% Check Unique Values in the same step (MERGE) %%%
% within the values of idxNEW (that must be previously arranged in the idxNEWtemp)

[idxNEWsorted,arrangement] = sort(idxNEWtemp);
PositiveBlock = idxNEWsorted > 0;

% Work only on Positive Block
idxNEWsPB = idxNEWsorted(PositiveBlock);
arrangementPB = arrangement(PositiveBlock);

idxRepeated = diff([0;idxNEWsPB]) == 0;
if sum(idxRepeated) ~= 0 % if there is at least one Repetition...
    idxRepeatedFIN = find(diff([idxRepeated;0]) == -1);
    idxRepeatedINI = find(flipud(diff([flipud(idxRepeated);0])) == -1);
    idxRepeatedPOS = idxRepeatedINI - 1;
    
    for j = 1 : length(idxRepeatedINI)
        
        idxNEWtemp(arrangementPB(idxRepeatedINI(j):idxRepeatedFIN(j))) = - ( arrangementPB(idxRepeatedPOS(j)) + numel(HistoryPath) );
        FlagNEWtemp(arrangementPB(idxRepeatedINI(j):idxRepeatedFIN(j))) = 1; % Flag MERGE
        
    end
end

clear idxNEWsorted arrangement idxNEWsPB arrangementPB idxRepeated* j


%%% Check Unique Values in the Others/Own past paths (CROSS) %%%
% among the unique values of idxNEW, is there any value that had already
% shown up in previous steps of the algorithm? - checking the previous
% OTHER paths (flag3) and its OWN previous path (flag5)

idxNEWtempMatrix = repmat(idxNEWtemp,1,size(HistoryPath,2));
% Concatenate the data: [HistoryPath;idxNEWtempMatrix]
CMatrix = cat(1,HistoryPath,idxNEWtempMatrix);
[CMsorted,CMarrangement] = sort(CMatrix,1);
PositiveBlock = CMsorted > 0;

L = size(HistoryPath,1);

% Work only on Positive Block
for j = 1 : size(CMsorted,2) % for each step of the previous paths
    
    % Consider only the Posivive values of the paths @ step j
    CMsPB = CMsorted(PositiveBlock(:,j),j);
    CMarrPB = CMarrangement(PositiveBlock(:,j),j);
    
    idxRepeated = find(diff(CMsPB) == 0);
    
    if ~isempty(idxRepeated)
        
        arrPBMatrix = [CMarrPB(idxRepeated),CMarrPB(idxRepeated + 1)];
        
        arrMin = min(arrPBMatrix,[],2);
        arrMax = max(arrPBMatrix,[],2);
        
        % Computing Flags
        FlagCheck = arrMax - arrMin;
        F_temp = 3*ones(size(FlagCheck)); % Initialization to Flag = 3 : Cross
        F_temp(FlagCheck == L) = 5; % Correcting Flag for: Loop
        
        % Assignment of vals to idxNEWtemp and FlagNEWtemp
        idxNEWtemp(arrMax - L) = - ( (j-1)*L + arrMin );
        FlagNEWtemp(arrMax - L) = F_temp;

        clear arrPBMatrix arrM* FlagCheck F_temp
        
    end
    
    clear CMsPB CMarrPB idxRepeated
    
end

% Update
HistoryPath = [HistoryPath,idxNEWtemp];
Flag = [Flag,FlagNEWtemp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CenteredSkelSeeds = AVS3D_ExtractCenteredSeeds(HP,Flag,Size)

% Looking ONLY for LOOPS (the end of the Paths)
Loopflag = 5;

Sel = Flag == Loopflag;
SelHP = HP(Sel,:);
% Get the 'first'(repeated) negative value (pointer)
[c_val,c_fin] = min(SelHP,[],2);
[  ~  ,c_ini] = ind2sub(size(HP), - c_val);

Seeds = zeros(size(c_ini));

for j = 1 : length(c_ini)
    
    [LoopX,LoopY,LoopZ] = ind2sub(Size,SelHP(j,c_ini(j):c_fin(j) - 1));
    
    LoopX = round(mean(LoopX));
    LoopY = round(mean(LoopY));
    LoopZ = round(mean(LoopZ));
    
    if LoopX < 1
        LoopX = 1;
    end
    if LoopY < 1
        LoopY = 1;
    end
    if LoopZ < 1
        LoopZ = 1;
    end
    
    
    if LoopX > Size(1)
        LoopX = Size(1);
    end
    if LoopY > Size(2)
        LoopY = Size(2);
    end
    if LoopZ > Size(3)
        LoopZ = Size(3);
    end
    
    Seeds(j) = sub2ind(Size,LoopX,LoopY,LoopZ);
    
end

CenteredSkelSeeds = false(Size);
CenteredSkelSeeds(Seeds) = true;

function [VEvec1NS,VEvec2NS,VEvec3NS,idxReal] = AVS3D_ComputeSelEVecs2(VG,idx,SortByMagnitude,AbsoluteFlag)

block_len = 16^3;
num_block = ceil(size(idx,1)/block_len);

VEvec1NS = NaN(size(idx,1),3);
VEvec2NS = NaN(size(idx,1),3);
VEvec3NS = NaN(size(idx,1),3);

idxRealCheck = false(size(idx));

for j = 1 : num_block
   
    if j < num_block
       
       ini = ((j-1)*block_len)+1;
       fin = (j*block_len);
       
       idx_block = idx( ini : fin );
   else
       
       ini = ((j-1)*block_len)+1;
       fin = size(idx,1);
       
       idx_block = idx( ini : fin );
   end
   
   [~,~,~,Evc1,Evc2,Evc3,idx_discard] = VTF3D_eigenDecompBlock3D_mex(VG,idx_block,SortByMagnitude,AbsoluteFlag);
   
   VEvec1NS(ini:fin,:) = Evc1;
   VEvec2NS(ini:fin,:) = Evc2;
   VEvec3NS(ini:fin,:) = Evc3;
   
   idxRealCheck(ini:fin) = ~idx_discard; % list of indices where Eigenvectors are to be considered ( i.e. PURE REAL, so Complex and NaN values are excluded)
   
   clear ini fin idx_block Evc*;
   
end

VEvec1NS = VEvec1NS(idxRealCheck,:);
VEvec2NS = VEvec2NS(idxRealCheck,:);
VEvec3NS = VEvec3NS(idxRealCheck,:);
idxReal = idx(idxRealCheck);

function OldValueStr = AVS3D_showProgressCM(ProgressMessage,OldValueStr,Value,Finished)

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
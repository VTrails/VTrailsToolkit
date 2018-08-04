function [GeodesicMSTs,GeodesicMSTsMatrix] = VTrailsRefine3D_ConnectedGraph2VTrailsTree_MAIN(DATA,varargin)%
%
%VTrails: Extract and Refine Geodesic Minimum Spanning Tree from
%         Over-Connected Vascular Graphs
%
%%%%%%%%%
% Inputs:
%
% - DATA: structure containing Image-related features, and the
%         Over-connected Vascular Graph.
%         NB: DATA is the output of the function:
%         'VTrailsConnectingGeodesicGraph3D_MAIN.m' 
%        
%                    DATA *MUST* contain the following Fields:
%
%  * IMG: Header of the input nifty Images.
%  * GRID: Structure with Image Grid and Scales' infos for Resampling.
%  * GGM: Geodesically Weighted Connectivity Matrix of the Over-Connected
%         Vascular Graph.
%  * exploredminPath: struct containing the list of all the extracted and
%                     explored connecting geodesics (minimal paths).
%  * MST: Minimum Spanning Tree(s) Matrix obtained from GGM with
%         'Kruskal' method.
%  * OrigSeedsIDX: list of Original Seeds;
%  * ExtraSeedsIDX: list of automatically generated Seeds, during the
%                   self-arranging refinement; 
%
% CONFIGURATION INPUTS:
%
% - ['PruningLengthThresholdMM']: Length of leaves that should be pruned. -
%                                 Default: 5 [millimeters]
% - ['PruningIterativeSteps']: Number of steps for the iterative pruning. -
%                              Default: 1
% - ['VicinityTHR_Leaves']: Threshold for spatially close leaves, that
%                           should be merged together. -
%                           Default: 1.5 [normalised units] 
% - ['VicinityTHR_Nodes']: Threshold for spatially close nodes, that
%                          should be welded together. -
%                          Default: 3 [normalised units]
% - ['LeavesPruningROI']: ROI-based Mask for removing sputious leaves after
%                         masking the Over-Connected Geodesic Graph
%                         accodring to anatomical vascular territories. -
%                         Default: [] (empty) -- This will be removed in
%                         future releases!
% - ['TreeRootSeed']: Coordinate of the user-defined Root of the Tree(s).
%                     [Root_X Root_Y Root_Z] array. - Default: [0 0 0] 
% - ['SupervisedSeeds']: Bool flag for user-defined vs. automatic seeds. -
%                        Default: true 
% - ['GroupIDs2ConsiderRange']: max value of output trees (in case of a
%                               Forest of MSTs) - Default: 1 -- Set to
%                               'inf' for all MSTs.
% - ['GroupIDs2Consider']: when the identifier(s) of the MSTs is known
%                          beforehand, use this parameter to select the
%                          MSTs. - Default: [] (empty).
% - ['PadSize']: number of padding samples for smoothing the 3D coordinates
%                of the extracted MSTs edges. - Default: 5
% - ['ResampleRatio']: resampling ratio for smoothing/resampling the 3D
%                      coordinates of the extracted MSTs edges. -
%                      Default: 0.1
% - ['PathFolder']: exporting path folder. - Default: as in DATA.PathFolder
%
%%%%%%%%%%
% Outputs:
%   
% - GeodesicMSTs: structure containing the 3D embedding, geodesic
%                 information, and edge-connectivity of the extracted
%                 Minimum Spanning Tree(s). 
%                 NB. in case multiple MSTs are obtained (Forest of MSTs),
%                 the same structure lists all the elements with respecive
%                 group-identifiers.
%
%      GeodesicMSTs will contain the following Edges-related Fields:
%
%  * CGPathContinuous: 3D sequence of coordinates for the conneted path.
%  * CGPathID: Identifying label of the path (for the considered MST).
%  * GeodesicEnergy: Integral value of the Geodesic Energy along the
%                    3D connecting path.
%  * GroupID: Identifier of the MST (in case of a forest of MSTs)
%  * CGPathAdjacencyList: List of connecting edges (or paths) with respect
%                         to the considered one.
%  * EuclideanLength: Integral Euclidean length of the 3D connecting path.
%  * isLeaf: Bool flag (true if the connecting edge is a leaf).
%  * isRoot: Bool flag (true if the tree ROOT belongs to that edge).
%
% - GeodesicMSTsMatrix: structure containing the connectivity matrix
%                       (adjacencies) for vectorial visualisation and
%                       representation.
%
%        GeodesicMSTsMatrix will contain the following Fields:
%
%  * M: Sparse matrix encoding the adjacency of the connected edges.
%
%%%%%%%%%%%
% Example Call:
%
%    [GeodesicMSTs,...
%     GeodesicMSTsMatrix] = VTrailsRefine3D_ConnectedGraph2VTrailsTree_MAIN( DATA , 'PruningLengthThresholdMM' , 5 );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  author: Stefano Moriconi
%  e-mail: vtrailstoolkit@gmail.com 
%  release: August 2018
%
%  The following code is part of VTrails toolkit
%  https://vtrails.github.io/VTrailsToolkit/
%
%  When using this code, please cite the corresponding references:  
%    
% [1] @InProceedings{Moriconi2017VTrails,
%  	  title = {VTrails: Inferring Vessels with Geodesic Connectivity Trees},
%  	  author = {Stefano Moriconi and 
%               Maria A. Zuluaga and 
%               Hans Rolf J{\"a}ger and 
%               Parashkev Nachev and 
%               S{\'e}bastien Ourselin and
%               M. Jorge Cardoso},
%  	  booktitle = {IPMI},
%  	  year = {2017}
%	  }
%
% [2] @Article{Moriconi2018Inference,
%  	  title = {Inference of Cerebrovascular Topology with Geodesic Minimum Spanning Trees},
%  	  author = {Stefano Moriconi and 
%               Maria A. Zuluaga and 
%               Hans Rolf J{\"a}ger and 
%               Parashkev Nachev and 
%               S{\'e}bastien Ourselin and
%               M. Jorge Cardoso},
%  	  journal  = {IEEE Transactions on Medical Imaging},
%  	  year = {2018}
%  	  }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

includeVTrailsLibs;

SET = VTR3D_configureINPUTS(DATA,varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting the Connected Geodesic Graph Matrix (GGM) and the
% exploredminPaths into the Minimum Spanning Tree: VTRAILS 
disp('_______________________________________________________');
disp(' * Converting Connected Geodesic Graph ...');
disp('-------------------------------------------------------');
VTRAILSc = VTR3D_computeGeodesicConnectingPaths(SET.GGM,SET.exploredminPath,SET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Leaves that have as Free-EndPoint an ExtraSeed
VTRAILSc = VTR3D_removeFreeEndPointsWithExtraSeeds(VTRAILSc,SET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Possible Signle Points
VTRAILSc = VTR3D_removePossibleSinglePoint(VTRAILSc,SET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Possible Leaves according to SET.LeavesPruningROI (for masked Geodesic Graphs)
VTRAILSc = removeFuzzyLeavesInPruningROI(VTRAILSc,SET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welding Edges that run along each others
VTRAILSc = VTR3D_performIterativeEdgeWelding(VTRAILSc,SET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the Trees by their GroupIDs (in case of a FOREST of Minimum
% Spanning Trees) 
GroupIDs2consider = VTR3D_retireveGIDs2Consider(VTRAILSc,SET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pruning VTRAILS (as many times as specified in 'PruningIterativeSteps'
VTRAILSprun = VTR3D_processPruning4VTRAILS(VTRAILSc,SET,GroupIDs2consider);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct the Final Tree (Adjacency Matrix)
[GeodesicMSTsMatrix,GeodesicMSTs] = VTR3D_recoverVTrailsConnectingMatrix(VTRAILSprun,SET.TreeRootSeed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoothing VTRAILS and assigning the Physical-Space Coordinates
if SET.Paraview
    GeodesicMSTs = VTR3D_smoothConnectedPaths4Paraview(GeodesicMSTs,SET);
else
    GeodesicMSTs = VTR3D_smoothConnectedPaths(GeodesicMSTs,SET);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exporting
VTR3D_exportRefinedVTrails(GeodesicMSTs,GeodesicMSTsMatrix,SET);

disp('*********************************************************');
disp('  VTrails: Vascular Minimum Spanning Tree -- COMPLETE!');


function SET = VTR3D_configureINPUTS(DATA,optInputs)
% Function to Assert Correctness of INPUTS
disp('*********************************************************');
disp('         VTrails: Vascular Minimum Spanning Tree ');
disp('*********************************************************');

assert(isstruct(DATA) ,'Input DATA MUST be the Struct generated by: "VTrailsConnectingGeodesicGraph3D_MAIN.m"');

% Initializing SET with Default Values
SET.IMG = DATA.IMG;
SET.Spacing = DATA.GRID.physicalvoxelsize ./ DATA.GRID.Scale;
SET.PhysicalImageOrigin = [DATA.IMG.hdr.hist.qoffset_x,DATA.IMG.hdr.hist.qoffset_y,DATA.IMG.hdr.hist.qoffset_z];
SET.ImageOriginRef = [DATA.IMG.hdr.dime.dim(2)*DATA.GRID.physicalvoxelsize(1) - DATA.IMG.hdr.hist.qoffset_x , ...
                      DATA.IMG.hdr.dime.dim(3)*DATA.GRID.physicalvoxelsize(2) - DATA.IMG.hdr.hist.qoffset_y , ...
                      DATA.IMG.hdr.hist.qoffset_z ];
SET.GGM = DATA.GGM;
SET.exploredminPath = DATA.exploredminPath;

[ESPosX,ESPosY,ESPosZ] = ind2sub([length(DATA.GRID.Xr),length(DATA.GRID.Yr),length(DATA.GRID.Zr)],DATA.ExtraSeedsIDX);
SET.ExtraSeedsPOS = [ESPosY,ESPosX,ESPosZ];

% NewFields
SET.TreeRootSeed = DATA.IMG.hdr.hist.originator(1:3);
SET.PruningLengthThreshold_mm = 5;
SET.PruningIterativeSteps = 1;
SET.GroupIDs2ConsiderRange = 1;
SET.GroupIDs2Consider = [];
SET.PathFolder = DATA.PathFolder;
SET.vicinityTHR_leaves = 1.5; % normalised units
SET.vicinityTHR_nodes = 3;
SET.SupervisedSeeds = true;

SET.PadSize = 5;
SET.ResampleRatio = 0.1;
SET.Paraview = false;
SET.ExportOffsetMultiplier = [0 0 0];
SET.ExportCoordMultiplier = [1 1 1];
SET.LeavesPruningROI = [];

if ~isempty(optInputs)
    for j = 1 : 2 :length(optInputs)
        
        switch upper(optInputs{j})
            case 'PRUNINGLENGTHTHRESHOLDMM'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given PruningLengthThreshold (mm) flag: -- Default: 5 applied!');
                else
                    SET.PruningLengthThreshold_mm = abs(optInputs{j+1});
                end
            case 'TREEROOTSEED'
                if isempty(optInputs{j+1}) || ~isvector(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || sum(isfinite(optInputs{j+1}))~=numel(optInputs{j+1}) || length(optInputs{j+1}) ~= 3
                    disp('Invalid|Not Given TreeRootSeed [Px Py Pz] array: -- Default: [0 0 0] applied!');
                else
                    SET.TreeRootSeed = optInputs{j+1};
                end
            case 'PRUNINGITERATIVESTEPS'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given PruningIterativeSteps value: -- Default: 1 applied!');
                else
                    SET.PruningIterativeSteps = abs(optInputs{j+1});
                end
            case 'GROUPIDS2CONSIDERRANGE'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || optInputs{j+1} < 1
                    disp('Invalid|Not Given GroupIDs2ConsiderRange value: -- Default: 1 applied!');
                else
                    SET.GroupIDs2ConsiderRange = optInputs{j+1};
                end
            case 'GROUPIDS2CONSIDER'
                if isempty(optInputs{j+1}) || ~isvector(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || sum(~isnan(optInputs{j+1}))~=numel(optInputs{j+1})
                    disp('Invalid|Not Given GroupIDs2Consider array: -- Default: [] applied!');
                else
                    SET.GroupIDs2Consider = abs(optInputs{j+1});
                end
            case 'PATHFOLDER'
                if isempty(optInputs{j+1})
                    disp(['Invalid|Not Given Path Folder for VTRAILS Exporting: -- Default: ', SET.PathFolder ]);
                else
                    if exist(optInputs{j+1},'file') ~= 7
                        mkdir(optInputs{j+1});
                    end
                    SET.PathFolder = optInputs{j+1};
                end
            case 'SUPERVISEDSEEDS'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given SupervisedSeeds flag: -- Default: "true" applied!');
                else
                    SET.SupervisedSeeds = optInputs{j+1};
                end
            case 'PADSIZE'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given PadSize value: -- Default: 5 applied!');
                else
                    SET.PadSize = abs(optInputs{j+1});
                    SET.PadSize(SET.PadSize == 0) = 5;
                end
            case 'RESAMPLERATIO'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given ResampleRatio value: -- Default: 0.1 applied!');
                else
                    SET.ResampleRatio = abs(optInputs{j+1});
                    SET.ResampleRatio(SET.ResampleRatio == 0) = 0.1;
                end
            case 'PARAVIEW'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given Paraview flag: -- Default: "false" applied!');
                else
                    SET.Paraview = optInputs{j+1};
                end
            case 'EXPORTOFFSETMULTIPLIER'
                if isempty(optInputs{j+1}) || ~isvector(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || sum(isfinite(optInputs{j+1}))~=numel(optInputs{j+1}) || length(optInputs{j+1}) ~= 3
                    disp('Invalid|Not Given ExportOffsetMultiplier [Mx My Mz] array: -- Default: [0 0 0] applied!');
                else
                    SET.ExportOffsetMultiplier = optInputs{j+1};
                end
            case 'EXPORTCOORDINATEMULTIPLIER'
                if isempty(optInputs{j+1}) || ~isvector(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || sum(isfinite(optInputs{j+1}))~=numel(optInputs{j+1}) || length(optInputs{j+1}) ~= 3
                    disp('Invalid|Not Given ExportCoordinateMultiplier [Cx Cy Cz] array: -- Default: [1 1 1] applied!');
                else
                    SET.ExportCoordMultiplier = optInputs{j+1};
                end
            case 'LEAVESPRUNINGROI'
                if isempty(optInputs{j+1}) || ~isequal(size(optInputs{j+1}),DATA.IMG.hdr.dime.dim(2:4))
                    disp('Invalid|Not Given LeavesPruningROI: -- Default: []');
                else
                    SET.LeavesPruningROI = logical(optInputs{j+1});
                end
            case 'VICINITYTHR_LEAVES'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given VicinityTHR_leaves value: -- Default: 1.5 [normUnits] applied!');
                else
                    SET.vicinityTHR_leaves = abs(optInputs{j+1});
                    SET.vicinityTHR_leaves(SET.vicinityTHR_leaves == 0) = 1.5;
                end
            case 'VICINITYTHR_NODES'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isfinite(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given VicinityTHR_nodes value: -- Default: 3.0 [norm units] applied!');
                else
                    SET.vicinityTHR_nodes = abs(optInputs{j+1});
                    SET.vicinityTHR_nodes(SET.vicinityTHR_nodes == 0) = 3.0;
                end
            otherwise
                disp('Given Unknown parameter(s) for OPT. DEFAULT value(s) will be set.');
        end
        
    end
end

function VTRAILSc = VTR3D_computeGeodesicConnectingPaths(GGM,ExploredMinPaths,SET)

PATHSc = VTR3D_convertGGM2ConnectedGeodesicPaths(GGM,ExploredMinPaths);
BRANCHESc = VTR3D_convertCGPaths2Branches(PATHSc);
VTRAILSc = VTR3D_convertBranches2VTrails(BRANCHESc,SET.Spacing);

function PATHSc = VTR3D_convertGGM2ConnectedGeodesicPaths(GGM,ExploredMinPaths)

continueprocessing = true;

GroupID = 1;

PATHSc = struct([]);

entry = length(PATHSc);

while continueprocessing
    
    if ~isequal(GGM,true([length(ExploredMinPaths),1]))
        
        MST = graphminspantree(sparse(GGM));
        
        [a1,a2,~] = find(MST);
        
        if isempty(a1) || isempty(a2)
            MST = graphminspantree(sparse(GGM),'Method','Kruskal');
            [a1,a2,~] = find(MST);
        end
        
        validpairs = zeros(size(ExploredMinPaths,2),2);
        for jj = 1 : size(ExploredMinPaths,2)
            validpairs(jj,:) = ExploredMinPaths(jj).pairs;
        end
        
        [~,elms2considerCHK,~] = VTR3D_repeatedValueSymm_mex(validpairs,[a1,a2]);
        
        elms2consider = find(elms2considerCHK);
        
        % Processing
        ProgShowFlag = length(elms2consider) >= 1;
        
        if ProgShowFlag
            if GroupID == 1
                ProgressOldString = '';
            end
            ProgressMessage = ['a) Computing Geodesic Paths -- ',num2str(GroupID,'%d'),'/? :'];
            ProgressValue = 0;
            ProgressValue_step = 100.0/size(elms2consider,1);
            ProgressOldString = VTR3D_showProgressCM(ProgressMessage,ProgressOldString,ProgressValue,false);
        end
        
        for elm = 1 : size(elms2consider,1)
            
            entry = entry + 1;
            
            PATHSc(entry).CGPathID = elm;
            PATHSc(entry).GroupID = GroupID;
            
            PATHSc(entry).CGPathContinuous = ExploredMinPaths(elms2consider(elm)).CGPathContinuous;
            PATHSc(entry).GeodesicEnergy = ExploredMinPaths(elms2consider(elm)).geodesicRAW;
            
            PATHSc(entry).CGPathAdjacencyList = VTR3D_getAdjacencyList(ExploredMinPaths,elms2consider,elm);
           
            if ProgShowFlag
                % Progress Percentage
                ProgressValue = ProgressValue + ProgressValue_step;
                ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, elm == size(elms2consider,1) );
            end
            
        end
        
        % Update GGM and GroupID
        GGM(a1,:) = 0; GGM(:,a1) = 0;
        GGM(a2,:) = 0; GGM(:,a2) = 0;
        GroupID = GroupID + 1;
        
        % Update while condition
        if isequal( graphminspantree(sparse(GGM)) , sparse(zeros(size(GGM))) )
            continueprocessing = false;
        end
        
        if ProgShowFlag
            ProgressOldString = repmat(sprintf('\b'), 1, length(ProgressMessage) + length(ProgressOldString));
        end
        
    else
        
        TotExploredMinPaths = ExploredMinPaths;
        
        GroupIDs = unique([TotExploredMinPaths(:).GroupID]);
        
        for gid = 1 : length(GroupIDs)
            
            ExploredMinPaths = TotExploredMinPaths([TotExploredMinPaths(:).GroupID] == GroupIDs(gid));
            elms2consider = (1:length(ExploredMinPaths))';
            
            % Processing
            ProgShowFlag = length(elms2consider) >= 1;
        
            if ProgShowFlag
                if gid == 1
                   ProgressOldString = '';
                end
                ProgressMessage = ['a) Computing Geodesic Paths -- ',num2str(gid,'%d'),'/',num2str(length(GroupIDs),'%d'),' :'];
                ProgressValue = 0;
                ProgressValue_step = 100.0/size(elms2consider,1);
                ProgressOldString = VTR3D_showProgressCM(ProgressMessage,ProgressOldString,ProgressValue,false);
            end
            
            for elm = 1 : size(elms2consider,1)
                
                entry = entry + 1;
                
                PATHSc(entry).CGPathID = elm;
                PATHSc(entry).GroupID = ExploredMinPaths(elms2consider(elm)).GroupID;
                
                PATHSc(entry).CGPathContinuous = ExploredMinPaths(elms2consider(elm)).CGPathContinuous;
                PATHSc(entry).GeodesicEnergy = ExploredMinPaths(elms2consider(elm)).GeodesicEnergy;
                
                PATHSc(entry).CGPathAdjacencyList = VTR3D_getAdjacencyList(ExploredMinPaths,elms2consider,elm);
                
                if ProgShowFlag
                    % Progress Percentage
                    ProgressValue = ProgressValue + ProgressValue_step;
                    ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, elm == size(elms2consider,1) );
                end
                
            end
            
            if ProgShowFlag
                ProgressOldString = repmat(sprintf('\b'), 1, length(ProgressMessage) + length(ProgressOldString));
            end
            
        end
        
        continueprocessing = false;
        
    end

end

function AdjLst = VTR3D_getAdjacencyList(ExploredMinPaths,elms2consider,elm)

MinPaths2consider = ExploredMinPaths(elms2consider);

refMinPath = ExploredMinPaths(elms2consider(elm));

points2eval = refMinPath.CGPathContinuous ; 

AdjLst = [];

for jj = 1 : size(elms2consider)
   
    if jj ~= elm % if it is NOT the same element (comparing points2eval against istelf!)
        
        otherpoints2eval = MinPaths2consider(jj).CGPathContinuous;
        
        if ~isempty(VTR3D_repeatedValue_mex(points2eval,otherpoints2eval)) % if there is a repetition in the endpoints of the connecting paths, then the Adjacency List MUST be updated!
            AdjLst = cat(1 , AdjLst , jj);
        end
        
    end
    
end

function BRANCHES_Final = VTR3D_convertCGPaths2Branches(TotPATHSc)

BRANCHES_Final = struct([]);

GroupIDs = unique([TotPATHSc(:).GroupID]);

for gid = 1 : length(GroupIDs)
    
    PATHSc = TotPATHSc([TotPATHSc(:).GroupID] == GroupIDs(gid));
    
    BRANCHES = struct([]);
    
    newEntry = 0;
    
    % Processing
    ProgShowFlag = size(PATHSc,2) >= 1;
        
    if ProgShowFlag
        if gid == 1
            ProgressOldString = '';
        end
        ProgressMessage = ['b) Computing Geodesic Paths -- ',num2str(gid,'%d'),'/',num2str(length(GroupIDs),'%d'),' :'];
        ProgressValue = 0;
        ProgressValue_step = 100.0/size(PATHSc,2);
        ProgressOldString = VTR3D_showProgressCM(ProgressMessage,ProgressOldString,ProgressValue,false);
    end
    
    for jj = 1 : size(PATHSc,2)
        
        if ~isempty(PATHSc(jj).CGPathContinuous)
            
            newEntry = newEntry + 1;
            
            mergedSeq = PATHSc(jj).CGPathContinuous;
            mergedGeodesic = PATHSc(jj).GeodesicEnergy;
            
            PATHSc(jj).CGPathContinuous = [];
            
            finalID = newEntry;
            mergedID = PATHSc(jj).CGPathID;
            unmergedID = [];
            
            Adjacent2check = PATHSc(jj).CGPathAdjacencyList;
            [~,A2C] = VTR3D_notUnique_mex(Adjacent2check,mergedID);
            Adjacent2check = unique(Adjacent2check(~A2C)); % Excluding those already visited
            
            while ~isempty(Adjacent2check)
                
                otherSeq = PATHSc(Adjacent2check(1)).CGPathContinuous;
                otherGeodesic = PATHSc(Adjacent2check(1)).GeodesicEnergy;
                
                
                if ~isempty(otherSeq)
                    
                    otherID = PATHSc(Adjacent2check(1)).CGPathID;
                    
                    [~,otherSeqChk,mergedSeqChk] = VTR3D_repeatedValue_mex(otherSeq,mergedSeq);
                    
                    posA = find(mergedSeqChk);
                    posB = find(otherSeqChk);
                    
                    if isscalar(posA) && isscalar(posB)
                    
                        if ((posB == 1) || (posB == size(otherSeqChk,1))) && ((posA == 1)|| (posA == size(mergedSeq,1)))% To be MERGED Together
                            
                            %%% MERGING THE SEQUENCES %%%
                            if     posB == 1 && posA == 1
                                
                                mergedSeq = [flipud(otherSeq);mergedSeq(2:end,:)];
                                mergedGeodesic = [flipud(otherGeodesic);mergedGeodesic(2:end,:)];
                                
                            elseif posB == 1 && posA == size(mergedSeq,1)                   
                                
                                % Compact
                                mergedSeq = cat(1,mergedSeq,otherSeq(2:end,:));
                                mergedGeodesic = cat(1,mergedGeodesic,otherGeodesic(2:end,:));
                                
                            elseif posB == size(otherSeq,1) && posA == 1
                                
                                mergedSeq = [otherSeq;mergedSeq(2:end,:)];
                                mergedGeodesic = [otherGeodesic;mergedGeodesic(2:end,:)];
                                
                            elseif posB == size(otherSeq,1) && posA == size(mergedSeq,1)
                                
                                mergedSeq = [mergedSeq(1:end-1,:);flipud(otherSeq)];
                                mergedGeodesic = [mergedGeodesic(1:end-1,:);flipud(otherGeodesic)];
                                
                            else
                                error('There MUST be an Error!');
                            end
                            
                            %%% Resetting the otherSeq %%%
                            
                            mergedID = sort([mergedID;otherID]);
                            
                            % Removing the otherSeq
                            PATHSc(Adjacent2check(1)).CGPathContinuous = [];
                            
                            % Updating the AdjacemncyList
                            Adjacent2check = unique([Adjacent2check(2:end);PATHSc(Adjacent2check(1)).CGPathAdjacencyList]);
                            [~,A2C] = VTR3D_notUnique_mex(Adjacent2check,mergedID);
                            Adjacent2check = unique(Adjacent2check(~A2C)); % Excluding those already visited
                            
                        else % To be LEFT As It Is
                            
                            unmergedID = unique([unmergedID;Adjacent2check(1)]);
                            Adjacent2check = Adjacent2check(2:end);
                            
                        end
                    else
                        warning('Entry with more than one point of overlap in the sequence!');
                        unmergedID = unique([unmergedID;Adjacent2check(1)]);
                        % Forget about the buggy entry!
                        Adjacent2check = Adjacent2check(2:end);
                    end
                else
                    unmergedID = unique([unmergedID;PATHSc(Adjacent2check(1)).LinkingNode]);
                    Adjacent2check = Adjacent2check(2:end);
                end
                
            end
            
            for kk = 1 : length(mergedID)
                PATHSc(mergedID(kk)).LinkingNode = finalID;
            end
            
            BRANCHES(finalID).Sequence = mergedSeq;
            BRANCHES(finalID).ID = finalID;
            BRANCHES(finalID).mergedID = mergedID;
            BRANCHES(finalID).AdjacencyList = unmergedID;
            BRANCHES(finalID).GeodesicEnergy = mergedGeodesic;
            
            clear mergedSeq mergedID remainingAdjacents mergedGeodesic;
            
        end
        
        if ProgShowFlag
            % Progress Percentage
            ProgressValue = ProgressValue + ProgressValue_step;
            ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, jj == size(PATHSc,2) );
        end
        
    end
    
    BRANCHES_Output = struct([]);
    
    for ii = 1 : size(BRANCHES,2)
        
        BRANCHES_Output(ii).CGPathContinuous = BRANCHES(ii).Sequence;
        BRANCHES_Output(ii).CGPathID = BRANCHES(ii).ID;
        BRANCHES_Output(ii).GeodesicEnergy = BRANCHES(ii).GeodesicEnergy;
        
        BRANCHES_Output(ii).GroupID = GroupIDs(gid);
        
        UpdtAdjLst = [];
        
        for jj = 1 : size(BRANCHES,2)
            
            if ii~=jj
                
                rep = VTR3D_repeatedValue_mex(BRANCHES(ii).Sequence,BRANCHES(jj).Sequence);
                
                if ~isempty(rep)
                    
                    assert(size(rep,1) == 1,'More than one Point In Common: There might be a Loop in the Graph(?)');
                    
                    UpdtAdjLst = cat(1,UpdtAdjLst,jj);
                    
                end
                
            end
            
        end
        
        BRANCHES_Output(ii).CGPathAdjacencyList = UpdtAdjLst;
        
    end
    
    if ProgShowFlag
        ProgressOldString = repmat(sprintf('\b'), 1, length(ProgressMessage) + length(ProgressOldString));
    end
    
    BRANCHES_Final = cat(2,BRANCHES_Final,BRANCHES_Output);

end

function VTRAILS = VTR3D_convertBranches2VTrails(TotBRANCHES,Spacing)

% Initialization
VTRAILS = struct([]);

GroupIDs = unique([TotBRANCHES(:).GroupID]);

ProgShowFlag = length(GroupIDs) >= 1;

for gid = 1 : length(GroupIDs)
    
    tempVTRAILS = TotBRANCHES([TotBRANCHES(:).GroupID] == GroupIDs(gid));

    % Processing
    if ProgShowFlag
        if gid == 1
            ProgressOldString = '';
        end
        ProgressMessage = ['c) Computing Geodesic Paths -- ',num2str(gid,'%d'),'/',num2str(length(GroupIDs),'%d'),' :'];
        ProgressValue = 0;
        ProgressValue_step = 100.0/6;
        ProgressOldString = VTR3D_showProgressCM(ProgressMessage,ProgressOldString,ProgressValue,false);
    end
    
    newVTRAILS = VTR3D_makeExtraCuts(tempVTRAILS);
    if ProgShowFlag
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, false );
    end
    
    newVTRAILS = VTR3D_checkConsistencyFullOverlap(newVTRAILS);
    if ProgShowFlag
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, false );
    end
    
    newVTRAILS = VTR3D_checkConsistencyPartOverlap(newVTRAILS);
    if ProgShowFlag
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, false );
    end
    
    newVTRAILS = VTR3D_updateAdjacencyList(newVTRAILS,Spacing);
    if ProgShowFlag
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, false );
    end
    
    newVTRAILS = VTR3D_removeUnconnectedIsles(newVTRAILS,Spacing);
    if ProgShowFlag
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, false );
    end
    
    newVTRAILS = VTR3D_assignIsLeaf(newVTRAILS);
    if ProgShowFlag
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTR3D_showProgressCM([],ProgressOldString,ProgressValue, true );
    end
    
    VTRAILS = cat(2, VTRAILS, newVTRAILS);
    
    if ProgShowFlag
        ProgressOldString = repmat(sprintf('\b'), 1, length(ProgressMessage) + length(ProgressOldString));
    end
    
end

function EuclideanLength = VTR3D_computeEuclideanLength(Sequence,Spacing)%
% Use this function to compute the Euclidean Length of a path constituted
% by the Ordered 'Sequence' of the indices in an Image of size equal to
% 'Size' and of spacing equal to 'spacing'

switch size(Sequence,2)
    case 2 % case 2D
        
        CGPtempCOORD = [ Spacing(1) .* Sequence(:,1) , Spacing(2) .* Sequence(:,2) ]; % Inverted XY in Matlab!
        
        EuclideanLength = sum( sqrt( sum( diff( CGPtempCOORD , 1 ) .^2 , 2) ) ) ;

    case 3 % case 3D
       
        CGPtempCOORD = [ Spacing(1) .* Sequence(:,1) , Spacing(2) .* Sequence(:,2) , Spacing(3) .* Sequence(:,3)]; 
        
        EuclideanLength = sum( sqrt( sum( diff( CGPtempCOORD , 1 ).^2 , 2 ) ) );
    otherwise
        disp('computeEuclideanLength supports only 2D or 3D indexed matrices.');
        EuclideanLength = [];
end

function newVTRAILS = VTR3D_checkConsistencyFullOverlap(tempVTRAILS)

% Creating Pairs
pairs = VTR3D_getPairs2BeEvaluated_mex(length(tempVTRAILS));

removelist = [];

for pp = 1 : size(pairs,1)
    
    [Rvals,Rvalspp1CHK,Rvalspp2CHK] = VTR3D_repeatedValue_mex(tempVTRAILS(pairs(pp,1)).CGPathContinuous,tempVTRAILS(pairs(pp,2)).CGPathContinuous);
    
    if ~isempty(Rvals) && ~isequal(size(Rvals,1),1)
        
        if isequal(Rvalspp1CHK,true(size(Rvalspp1CHK))) && isequal(Rvalspp2CHK,true(size(Rvalspp2CHK))) % they are both EQUAL!
            % Remove one (the second)!
            removelist = cat(1,removelist,pairs(pp,2));
            
        elseif isequal(Rvalspp1CHK,true(size(Rvalspp1CHK))) && ~isequal(Rvalspp2CHK,true(size(Rvalspp2CHK))) % the first is fully included in the second
            % Remove the first one!
            removelist = cat(1,removelist,pairs(pp,1));
            
        elseif isequal(Rvalspp2CHK,true(size(Rvalspp2CHK))) && ~isequal(Rvalspp1CHK,true(size(Rvalspp1CHK))) % the second is fully included in the first
            % Remove the second one!
            removelist = cat(1,removelist,pairs(pp,2));
            
        end
        
    end

end

removelist = unique(removelist);

if ~isempty(removelist)
    disp('Full-Overlap Found! -- Please Check your data!');
end

[~,entryCHK] = VTR3D_notUnique_mex(1:length(tempVTRAILS),removelist);
newVTRAILS = tempVTRAILS(~entryCHK);

function newVTRAILS = VTR3D_checkConsistencyPartOverlap(tempVTRAILS)

appendVTRAILS = struct([]);

entry = length(tempVTRAILS);

% Creating Pairs
pairs = VTR3D_getPairs2BeEvaluated_mex(length(tempVTRAILS));

removelist = [];

for pp = 1 : size(pairs,1)
    
    [Rvals,Rvalspp1CHK,Rvalspp2CHK] = VTR3D_repeatedValue_mex(tempVTRAILS(pairs(pp,1)).CGPathContinuous,tempVTRAILS(pairs(pp,2)).CGPathContinuous);
    
    if ~isempty(Rvals) && ~isequal(size(Rvals,1),1)
        
        if ~isequal(Rvalspp1CHK,true(size(Rvalspp1CHK))) && ~isequal(Rvalspp2CHK,true(size(Rvalspp2CHK))) % Partial Overlap...
            
            removelist = cat(1,removelist,pairs(pp,1),pairs(pp,2));
            
            tmpVTRLS1 = tempVTRAILS(pairs(pp,1));
            tmpVTRLS1.CGPathContinuous = tmpVTRLS1.CGPathContinuous( ~Rvalspp1CHK  , :);
            tmpVTRLS1.GeodesicEnergy = tmpVTRLS1.GeodesicEnergy( ~Rvalspp1CHK  , :);
            
            tmpVTRLS2 = tempVTRAILS(pairs(pp,2));
            tmpVTRLS2.CGPathContinuous = tmpVTRLS2.CGPathContinuous( ~Rvalspp2CHK  , :);
            tmpVTRLS2.GeodesicEnergy = tmpVTRLS2.GeodesicEnergy( ~Rvalspp2CHK  , :);
            
            tmpVTRLS3 = tempVTRAILS(pairs(pp,1));
            tmpVTRLS3.CGPathContinuous = tmpVTRLS3.CGPathContinuous( imdilate ( Rvalspp1CHK ,logical([1 1 1]') ), :);
            tmpVTRLS3.GeodesicEnergy = tmpVTRLS3.GeodesicEnergy( imdilate ( Rvalspp1CHK ,logical([1 1 1]') ), :);
            
            tmpVTRLS3.isLeaf = false;
            entry = entry + 1;
            tmpVTRLS3.CGPathID = entry;
            tmpVTRLS3.Parent = NaN;
            
            appendVTRAILS = cat(2,appendVTRAILS,tmpVTRLS1,tmpVTRLS2,tmpVTRLS3);
            
        end
        
    end
    
end
    

removelist = unique(removelist);

if ~isempty(removelist)
    disp('Partial-Overlap Found! -- Please Check your data!');
end

[~,entryCHK] = VTR3D_notUnique_mex(1:length(tempVTRAILS),removelist);
newVTRAILS = tempVTRAILS(~entryCHK);
newVTRAILS = cat(2,newVTRAILS,appendVTRAILS);

function newVTRAILS = VTR3D_makeExtraCuts(tempVTRAILS)

removelist = [];

appendVTRAILS = struct([]);

entry = length(tempVTRAILS);

for elm = 1 : size(tempVTRAILS,2)
    
    CutPositions = [];
    
    for adj = 1 : size(tempVTRAILS(elm).CGPathAdjacencyList,1)
        
        [Rvals,Rvalspp1CHK,~] = VTR3D_repeatedValue_mex( tempVTRAILS(elm).CGPathContinuous , tempVTRAILS(tempVTRAILS(elm).CGPathAdjacencyList(adj)).CGPathContinuous );
        
        if ~isempty(Rvals) && isequal( size(Rvals,1) , 1 )
            
            endpts1 = (find(Rvalspp1CHK) == 1 || find(Rvalspp1CHK) == length(Rvalspp1CHK));
            
            if ~endpts1

                removelist = unique( cat(1, removelist, elm ) );
                
                % first Path Split
                CutPositions = sort( cat(1, CutPositions, find(Rvalspp1CHK) ) );
                
            end
                
        end
            
    end
    
    if ~isempty( CutPositions )
        
        tmpVTRLS = struct([]);
        
        CutPositions = sort( cat(1, CutPositions, 1 , length(Rvalspp1CHK) ) );
        
        % Determine Here the New Segments According to CutPositions
        
        for cut = 1 : length(CutPositions) - 1
            
            entry = entry + 1;
            tmpVTRLS(cut).CGPathID = entry;
            tmpVTRLS(cut).CGPathContinuous = tempVTRAILS(elm).CGPathContinuous( CutPositions(cut) : CutPositions(cut+1) , : );
            tmpVTRLS(cut).GeodesicEnergy = tempVTRAILS(elm).GeodesicEnergy( CutPositions(cut) : CutPositions(cut+1) , : );
            tmpVTRLS(cut).CGPathAdjacencyList = NaN;
            tmpVTRLS(cut).GroupID = tempVTRAILS(elm).GroupID;
            
        end
        
        appendVTRAILS = cat(2,appendVTRAILS,tmpVTRLS);
        
        clear tmpVTRLS;
    end

end

removelist = unique(removelist);

[~,entryCHK] = VTR3D_notUnique_mex(1:length(tempVTRAILS),removelist);
newVTRAILS = tempVTRAILS(~entryCHK);
newVTRAILS = cat(2,newVTRAILS,appendVTRAILS);

function newVTRAILS = VTR3D_updateAdjacencyList(newVTRAILS,Spacing)

% Resetting Values Before Update
for elm = 1 : size(newVTRAILS,2)
   
    newVTRAILS(elm).CGPathID = elm;
    newVTRAILS(elm).CGPathAdjacencyList = [];
    newVTRAILS(elm).EuclideanLength = VTR3D_computeEuclideanLength(newVTRAILS(elm).CGPathContinuous,Spacing);
    
end

% Creating Pairs
pairs = VTR3D_getPairs2BeEvaluated_mex(length(newVTRAILS));

for pp = 1 : size(pairs,1)
    
    [Rvals,~,~] = VTR3D_repeatedValue_mex(newVTRAILS(pairs(pp,1)).CGPathContinuous , newVTRAILS(pairs(pp,2)).CGPathContinuous );
    
    if ~isempty(Rvals)
        newVTRAILS(pairs(pp,1)).CGPathAdjacencyList = cat(1, newVTRAILS(pairs(pp,1)).CGPathAdjacencyList , pairs(pp,2) );
        newVTRAILS(pairs(pp,2)).CGPathAdjacencyList = cat(1, newVTRAILS(pairs(pp,2)).CGPathAdjacencyList , pairs(pp,1) );
    end
    
end

function newVTRAILS = VTR3D_removeUnconnectedIsles(newVTRAILS,Spacing)

if length(newVTRAILS) > 1 % Evaluate only if the Tree has more than one connected Branch
    
    Paths2BeKept = true(1,length(newVTRAILS));
    
    for elm = 1 : length(newVTRAILS)
        
        if isempty(newVTRAILS(elm).CGPathAdjacencyList)
            
            Paths2BeKept(elm) = false;
            
        end
        
    end
    
    newVTRAILS = newVTRAILS(Paths2BeKept);
    
    if ~isequal(sum(Paths2BeKept),length(Paths2BeKept)) % Update only if some Unconnected Islands have been removed!
        newVTRAILS = VTR3D_updateAdjacencyList(newVTRAILS,Spacing);
    end
    
end

function newVTRAILS = VTR3D_assignIsLeaf(newVTRAILS)

% Creating Pairs
pairs = VTR3D_getPairs2BeEvaluated_mex(length(newVTRAILS));

EndpointsConnectionCounter = struct('ep1c',[],'ep2c',[]);

for pp = 1 : size(pairs,1)
    
    EndpointsConnectionCounter(pairs(pp,1)).ep1 = newVTRAILS(pairs(pp,1)).CGPathContinuous( 1  , : );
    EndpointsConnectionCounter(pairs(pp,1)).ep2 = newVTRAILS(pairs(pp,1)).CGPathContinuous(end , : );
    
    points1_toEval = [newVTRAILS(pairs(pp,1)).CGPathContinuous( 1  , : );newVTRAILS(pairs(pp,1)).CGPathContinuous(end , : )];
    
    EndpointsConnectionCounter(pairs(pp,2)).ep1 = newVTRAILS(pairs(pp,2)).CGPathContinuous( 1  , : );
    EndpointsConnectionCounter(pairs(pp,2)).ep2 = newVTRAILS(pairs(pp,2)).CGPathContinuous(end , : );
    
    points2_toEval = [newVTRAILS(pairs(pp,2)).CGPathContinuous( 1  , : );newVTRAILS(pairs(pp,2)).CGPathContinuous(end , : )];
    
    [~,p1tE,p2tE] = VTR3D_repeatedValue_mex(points1_toEval,points2_toEval);
    
    EndpointsConnectionCounter(pairs(pp,1)).ep1c = cat(1,EndpointsConnectionCounter(pairs(pp,1)).ep1c,p1tE(1));
    EndpointsConnectionCounter(pairs(pp,1)).ep2c = cat(1,EndpointsConnectionCounter(pairs(pp,1)).ep2c,p1tE(2));
    
    EndpointsConnectionCounter(pairs(pp,2)).ep1c = cat(1,EndpointsConnectionCounter(pairs(pp,2)).ep1c,p2tE(1));
    EndpointsConnectionCounter(pairs(pp,2)).ep2c = cat(1,EndpointsConnectionCounter(pairs(pp,2)).ep2c,p2tE(2));
    
end

for elm = 1 : size(EndpointsConnectionCounter,2)
   
    SUMep1c = sum(EndpointsConnectionCounter(elm).ep1c);
    SUMep2c = sum(EndpointsConnectionCounter(elm).ep2c);
    
    if SUMep1c < 1 || SUMep2c < 1
        newVTRAILS(elm).isLeaf = true;
    else
        newVTRAILS(elm).isLeaf = false;
    end
    
end

function prunedVTRAILS = VTR3D_pruneVTrailsLeaves(VTRAILS,GroupID,SET)

euclidTHR = SET.PruningLengthThreshold_mm;

prunedVTRAILS = struct([]);

for gID = 1 : length(GroupID)
   
    selVTRAILS = VTRAILS([VTRAILS(:).GroupID] == GroupID(gID));
    
    ELs = [selVTRAILS(:).EuclideanLength];
    Leaves = [selVTRAILS(:).isLeaf];
    Root = [selVTRAILS(:).isRoot];
    
    if euclidTHR == 0
        seleuclidTHR = quantile(ELs(Leaves),0.5);
    else
        seleuclidTHR = euclidTHR;
    end
    
    if ~isempty(SET.LeavesPruningROI)
        isInPruningROI = getLeavesInPruningROI(selVTRAILS,SET.LeavesPruningROI,0.25);
    else
        isInPruningROI = false(size(Leaves));
    end
    
    isTooShort = ELs < seleuclidTHR;
    
    ExclusionCriterion = ( Leaves(:) & isTooShort(:) & ~Root(:) ) | ( Leaves(:) & ~Root(:) & isInPruningROI(:) );

    selVTRAILS = selVTRAILS(~ExclusionCriterion);

    % Concatenation
    prunedVTRAILS = cat(2,prunedVTRAILS,selVTRAILS);
    
end

function hasFreeEndPointExtraSeed = VTR3D_checkExtraSeedsAtEndPoint(gidVTRAILSc,ExtraSeedsPOS)

hasFreeEndPointExtraSeed = false(length(gidVTRAILSc),1);

for el = 1 : length(gidVTRAILSc)
    
    if gidVTRAILSc(el).isLeaf
        
        [repES,~,~] = VTR3D_repeatedValueSymm_mex(gidVTRAILSc(el).CGPathContinuous,ExtraSeedsPOS);

        if ~isempty(repES)
            
            % the ExtraSeed is located in the Sequence:
            % Check now whether it is connected to the rest of the Tree, or
            % if it a Free-EndPoint (no other Connections with the
            % Sequences of the AdjacencyList)
            
            repESchkCAT = [];
            
            for cc = 1 : length(gidVTRAILSc(el).CGPathAdjacencyList)
                
                [~,repESchk,~] = VTR3D_repeatedValue_mex(repES,gidVTRAILSc( gidVTRAILSc(el).CGPathAdjacencyList(cc) ).CGPathContinuous);
                
                repESchkCAT = cat(2,repESchkCAT,repESchk);
                
            end
            
            % If all ExtraSeeds are Connected to the Tree, then the
            % repESchkCAT = true(size(repESchkCAT)
            if sum(repESchkCAT(:)) ~= numel(repESchkCAT)
                hasFreeEndPointExtraSeed(el) = true;
            end
            
        end
        
    end
   
end

function VTRAILSwoFEPES = VTR3D_removeFreeEndPointsWithExtraSeeds(VTRAILS,SET)

if ~(SET.SupervisedSeeds)
    disp('_______________________________________________________');
    disp(' * Removing Leaves with ExtraSeed as Free-EndPoint ...');
    
    VTRAILS = VTR3D_assignTreeRoot(VTRAILS,SET.TreeRootSeed);
    
    VTRAILSwoFEPES = struct([]);
    
    ExtraSeedsPOS = SET.ExtraSeedsPOS;
    
    GroupIDs = unique([VTRAILS(:).GroupID]);
    
    for gID = 1 : length(GroupIDs)
        
        selVTRAILS = VTRAILS([VTRAILS(:).GroupID] == GroupIDs(gID));
        
        Leaves = [selVTRAILS(:).isLeaf];
        Root = [selVTRAILS(:).isRoot];
        
        hasFreeEndPointExtraSeed = VTR3D_checkExtraSeedsAtEndPoint(selVTRAILS,ExtraSeedsPOS);
        
        ExclusionCriterion = ( Leaves(:) & hasFreeEndPointExtraSeed(:) & ~Root(:) );
        
        % KEEP those that MUST NOT be EXCLUDED
        selVTRAILS = selVTRAILS(~ExclusionCriterion);
        
        % Concatenation
        VTRAILSwoFEPES = cat(2,VTRAILSwoFEPES,selVTRAILS);
        
    end
else
    VTRAILSwoFEPES = VTRAILS;
end

function [VTRAILStree,VTRAILScat] = VTR3D_recoverVTrailsConnectingMatrix(VTRAILS,TreeRootSeed)

VTRAILStree = struct([]);
VTRAILScat = struct([]);

GroupIDs = unique([VTRAILS(:).GroupID]);

for gid = 1 : length(GroupIDs)
    
    selVTRAILS = VTRAILS([VTRAILS(:).GroupID] == GroupIDs(gid));
    
    d = zeros(length(selVTRAILS),1);
    
    % new version for MarbalR2015b up
    Madj = zeros([length(selVTRAILS),length(selVTRAILS)]);
    
    for el = 1 : length(selVTRAILS)
        
        d(el) = min( sqrt( sum( ( selVTRAILS(el).CGPathContinuous - repmat(TreeRootSeed, [ size(selVTRAILS(el).CGPathContinuous,1), 1] ) ).^2 , 2) ) );
        
        selVTRAILS(el).isRoot = false;
        
        selVTRAILS(el).GeodesicEnergy = log( 1 + sum(selVTRAILS(el).GeodesicEnergy) );
        
        Madj(selVTRAILS(el).CGPathID,selVTRAILS(el).CGPathAdjacencyList) = 1;
        Madj(selVTRAILS(el).CGPathAdjacencyList,selVTRAILS(el).CGPathID) = 1;
        
    end
    
    Leaves = [selVTRAILS(:).isLeaf];
    d(~Leaves) = inf;
    
    [~,RootID] = min(d);
    
    selVTRAILS(RootID).isRoot = true;
    
    VTRAILStree(gid).M = adjacency( shortestpathtree( graph( Madj ) , RootID ) );
    
    VTRAILScat = cat(2,VTRAILScat,selVTRAILS);
    
end

function [VTRAILStree,visitedList] = VTR3D_getChildrenNodes(VTRAILS,fatherNode,VTRAILStree,visitedList)

for ff = 1 : length(fatherNode)
    
    continueSearch = true;
    
    while continueSearch
        
        visitedList = unique(cat(1,visitedList,fatherNode(ff)));
        
        children = VTRAILS(fatherNode(ff)).CGPathAdjacencyList;
        [~,childrenCHK] = VTR3D_notUnique_mex(children,visitedList);
        
        children = children(~childrenCHK);
        
        if ~isempty(children)
            VTRAILStree(fatherNode(ff),children) = 1;
            visitedList = unique(cat(1,visitedList,children(:)));
            [VTRAILStree,visitedList] = VTR3D_getChildrenNodes(VTRAILS,children,VTRAILStree,visitedList);
        else
            continueSearch = false;
        end
        
    end
    
end

function gdp_Uniform = VTR3D_smoothResampleSequence(gdp,padsize,UniformSpacing)

gdp = padarray(gdp,[padsize,0],'replicate');
gWin = gausswin((2*padsize)-1)./sum(gausswin((2*padsize)-1));
gdp_filttemp = [ conv(gdp(:,1),gWin,'same') , conv(gdp(:,2),gWin,'same') , conv(gdp(:,3),gWin,'same') ];
gdp = [gdp(1,:) ; gdp_filttemp(padsize+1:end-padsize,:) ; gdp(end,:)];

if size(gdp,1) > 5
    
    tt = 1:size(gdp,1);
    ratio = UniformSpacing/mean([mean(abs(diff(gdp(:,1)))),...
                                 mean(abs(diff(gdp(:,2)))),...
                                 mean(abs(diff(gdp(:,3))))]);
    
    tt_new = linspace(tt(1),tt(end),round((tt(end)-tt(1))./ratio)+1);
    
    gdp_Uniform = [linterp_matlab(tt,gdp(:,1),tt_new)',...
                   linterp_matlab(tt,gdp(:,2),tt_new)',...
                   linterp_matlab(tt,gdp(:,3),tt_new)'];
    
else
    gdp_Uniform = gdp;
end

function VTRAILSroot = VTR3D_assignTreeRoot(VTRAILS,TreeRootSeed)

VTRAILSroot = struct([]);

GroupIDs = unique([VTRAILS(:).GroupID]);

for gid = 1 : length(GroupIDs)
    
    selVTRAILS = VTRAILS([VTRAILS(:).GroupID] == GroupIDs(gid));
    
    d = zeros(length(selVTRAILS),1);
    
    for el = 1 : length(selVTRAILS)
        
        d(el) = min( sqrt( sum( ( selVTRAILS(el).CGPathContinuous - repmat(TreeRootSeed, [ size(selVTRAILS(el).CGPathContinuous,1), 1] ) ).^2 , 2) ) );
        
        selVTRAILS(el).isRoot = false;
        
    end
    
    Leaves = [selVTRAILS(:).isLeaf];
    d(~Leaves) = inf;
    
    [~,RootID] = min(d);
    
    selVTRAILS(RootID).isRoot = true;

    VTRAILSroot = cat(2,VTRAILSroot,selVTRAILS);
    
end

function VTR3D_exportRefinedVTrails(GeodesicMSTs,GeodesicMSTsMatrix,SET)

IMGhdr = SET.IMG;

PathDirName = SET.PathFolder;

DATAFileName = strcat([PathDirName,SET.IMG.ImgName,'_GeodesicMSTs.mat']);

disp([' * Exporting VTrails - Vascular Minimum Spanning Tree to: ',DATAFileName]);

save(DATAFileName,'GeodesicMSTs','GeodesicMSTsMatrix','IMGhdr','-v7.3');

function GroupIDs2consider = VTR3D_retireveGIDs2Consider(VTRAILSc,SET)

if isempty(SET.GroupIDs2Consider)
    
    GroupIDs = unique([VTRAILSc(:).GroupID]);
    CumLengths = zeros(size(GroupIDs));
    for gid = 1 : length(GroupIDs)
        CumLengths(gid) = sum([VTRAILSc([VTRAILSc(:).GroupID] == GroupIDs(gid)).EuclideanLength]);
    end
    [~,CLidx] = sort(CumLengths,'descend');
    if length(CLidx)>1
        GroupIDs2consider = GroupIDs(CLidx(1:SET.GroupIDs2ConsiderRange));%%% !!! %%% GroupID Selection!
    else
        GroupIDs2consider = GroupIDs(CLidx(1));
    end
    
else
    if isinf(SET.GroupIDs2Consider)
        GroupIDs2consider = unique([VTRAILSc(:).GroupID]);
    else
        GroupIDs2consider = SET.GroupIDs2Consider;
    end
end

function VTRAILSprun = VTR3D_processPruning4VTRAILS(VTRAILSc,SET,GroupIDs2consider)
disp('_______________________________________________________');
disp(' * Pruning and Refining Minimum Spanning Tree ...');

if SET.PruningIterativeSteps > 0
    
    tempVTRAILSc = VTRAILSc;
    
    for pp = 1 : SET.PruningIterativeSteps
        % Assign The Tree Root (ot Forest Roots) so that it cannot be
        % pruned!
        tempVTRAILSc = VTR3D_assignTreeRoot(tempVTRAILSc,SET.TreeRootSeed);
        
        % Pruning
        VTRAILScPruned = VTR3D_pruneVTrailsLeaves( tempVTRAILSc , GroupIDs2consider , SET );
        
        SET.PruningLengthThreshold_mm = SET.PruningLengthThreshold_mm * (1 + (pp/SET.PruningIterativeSteps)/10);
        
        % Recovening back the fully connected TREE
        disp('-------------------------------------------------------');
        tempVTRAILSc = VTR3D_computeGeodesicConnectingPaths(true([length(VTRAILScPruned),1]),VTRAILScPruned,SET);

        % Remove Degenerate Leaves
        tempVTRAILSc = VTR3D_removeDegenerateEntries(tempVTRAILSc);
        
    end
    
    VTRAILSprun = tempVTRAILSc;
    clear tempVTRAILSc VTRAILScPruned tempPATHSc tempBRANCHESc;
    
else
    VTRAILScPruned = [];
    
    for gid = 1 : length(GroupIDs2consider)
        
        VTRAILStemp = VTRAILSc([VTRAILSc(:).GroupID] == GroupIDs2consider(gid));
        
        VTRAILScPruned = cat(2, VTRAILScPruned , VTRAILStemp );
    end
    
    % Recovening back the fully connected TREE
    disp('-------------------------------------------------------');
    VTRAILScPruned = VTR3D_computeGeodesicConnectingPaths(true([length(VTRAILScPruned),1]),VTRAILScPruned,SET);

    % Remove Degenerate Leaves
    VTRAILScPruned = VTR3D_removeDegenerateEntries(VTRAILScPruned);
    
    VTRAILSprun = VTRAILScPruned;
end

function VTRAILSprun = VTR3D_smoothConnectedPaths(VTRAILSprun,SET)

disp('_______________________________________________________');
disp(' * Smoothing Final VTrails ...');

padsize = SET.PadSize;
resampleRatio = SET.ResampleRatio;%0.1;

ref = SET.IMG.hdr.dime.pixdim(2:4) .* SET.IMG.hdr.dime.dim(2:4);
EOM = SET.ExportOffsetMultiplier;
ECM = SET.ExportCoordMultiplier;

for elm = 1 : size(VTRAILSprun,2)
   
    VTRAILSprun(elm).CGPathContinuous = VTR3D_smoothResampleSequence(VTRAILSprun(elm).CGPathContinuous, padsize , resampleRatio ); 
    
%    %%% Uncomment here for the physical space coordinates! %%%
%     VTRAILSprun(elm).CGPathContinuous = [ SET.ImageOriginRef(2) - ((VTRAILSprun(elm).CGPathContinuous(:,2)) * SET.Spacing(2)) , ...
%                                          SET.ImageOriginRef(1) - ((VTRAILSprun(elm).CGPathContinuous(:,1)) * SET.Spacing(1)) , ...
%                                          SET.ImageOriginRef(3) + ((VTRAILSprun(elm).CGPathContinuous(:,3)) * SET.Spacing(3)) ];
%     
%     VTRAILSprun(elm).CGPathContinuous = [ (EOM(1) * ref(1)) + ECM(1) * VTRAILSprun(elm).CGPathContinuous(:,1),...
%                                          (EOM(2) * ref(2)) + ECM(2) * VTRAILSprun(elm).CGPathContinuous(:,2),...
%                                          (EOM(3) * ref(3)) + ECM(3) * VTRAILSprun(elm).CGPathContinuous(:,3)];
    
    VTRAILSprun(elm).EuclideanLength = sum( sqrt( sum( diff( VTRAILSprun(elm).CGPathContinuous ).^2 , 2) ) );
    
end

function VTRAILSprun = VTR3D_smoothConnectedPaths4Paraview(VTRAILSprun,SET)

disp('_______________________________________________________');
disp(' * Smoothing Final VTrails ...');

padsize = SET.PadSize;
resampleRatio = SET.ResampleRatio;
EOM = SET.ExportOffsetMultiplier;

for elm = 1 : size(VTRAILSprun,2)
   
    VTRAILSprun(elm).CGPathContinuous = VTR3D_smoothResampleSequence(VTRAILSprun(elm).CGPathContinuous, padsize , resampleRatio ); 
    
    VTRAILSprun(elm).CGPathContinuous = [ (EOM(1) * SET.IMG.hdr.hist.qoffset_x) + ((VTRAILSprun(elm).CGPathContinuous(:,2)) * SET.Spacing(2)) , ...
                                          (EOM(2) * SET.IMG.hdr.hist.qoffset_y) + ((VTRAILSprun(elm).CGPathContinuous(:,1)) * SET.Spacing(1)) , ...
                                          (EOM(3) * SET.IMG.hdr.hist.qoffset_z) + ((VTRAILSprun(elm).CGPathContinuous(:,3)) * SET.Spacing(3)) ];
    
    VTRAILSprun(elm).EuclideanLength = sum( sqrt( sum( diff( VTRAILSprun(elm).CGPathContinuous ).^2 , 2) ) );
    
end

%%% Section: Edge Welding

function VTRAILSconverged = VTR3D_performIterativeEdgeWelding(VTRAILSc,SET)

disp('_______________________________________________________');
disp(' * Welding Edges running along each others ...');

ContinueWelding = true;
weld_itr = 0;
maxitr = 15;

VTRAILSconverged = struct([]);

while ContinueWelding && weld_itr <= maxitr
    
    weld_itr = weld_itr + 1;
    
    disp('-------------------------------------------------------');
    VTRAILSc = VTR3D_computeGeodesicConnectingPaths( true([length(VTRAILSc),1]) , VTRAILSc , SET);
    
    gIDs = unique([VTRAILSc(:).GroupID]);
    tempLLs = zeros(1,length(gIDs));
    for gid = 1 : length(gIDs)
        sIDsCHK = [VTRAILSc(:).GroupID] == gIDs(gid);
        tempLLs(1,gid) = sum(sIDsCHK);
    end
        
    if weld_itr > 1
        
        LLs = cat(1,LLs,tempLLs);
        
        % Check for change state of ContinueWelding
        Converged = LLs(2,:) >= LLs(1,:);
        
        VTRAILStoProcess = struct([]);
        
        for gid = 1 : length(gIDs)
            if Converged(gid)
                tempVTRAILS = VTRAILSc([VTRAILSc(:).GroupID] == gIDs(gid));
                VTRAILSconverged = cat(2,VTRAILSconverged,tempVTRAILS);
            else
                tempVTRAILS = VTRAILSc([VTRAILSc(:).GroupID] == gIDs(gid));
                VTRAILStoProcess = cat(2,VTRAILStoProcess,tempVTRAILS);
            end
        end
        
        LLs = LLs(2,~Converged);
        
    else
        LLs = tempLLs;
        VTRAILStoProcess = VTRAILSc;
    end
    
    if ~isempty(VTRAILStoProcess) && weld_itr < maxitr
        % Edge Welding
        [VTRAILSc,~] = VTR3D_refineByWeldingCloseEdges(VTRAILStoProcess,SET.vicinityTHR_leaves,SET.Spacing);
        % Restore Structure
        disp('-------------------------------------------------------');
        VTRAILSc = VTR3D_computeGeodesicConnectingPaths( true([length(VTRAILSc),1]) , VTRAILSc , SET);
        % Node Welding
        [VTRAILSc,~] = VTR3D_refineByWeldingCloseNodes(VTRAILSc,SET.vicinityTHR_nodes,SET.Spacing);
    else
        ContinueWelding = false;
    end
    
end

if weld_itr >= maxitr
    warning('Iterative Edge Welding Process: maxitr REACHED!');
end

function [VTRAILS,ChangeOccurred] = VTR3D_refineByWeldingCloseEdges(VTRAILSc,vicinityTHRleaves,Spacing)

VTRAILS = struct([]);

GroupIDs = unique([VTRAILSc(:).GroupID]);

ChangeOccurred = false(1,length(GroupIDs));
VerifyStruct = false(1,length(GroupIDs));

% Processing Each Tree (in case of a Forest) according to the GroupID
for gid = 1 : length(GroupIDs)
    
    tempVTRAILSc = VTRAILSc([VTRAILSc(:).GroupID] == GroupIDs(gid));
    
    [procVTRAILSc,ChangeOccurred(gid)] = VTR3D_refineEdgesEachTree(tempVTRAILSc,vicinityTHRleaves,Spacing);
    
    VerifyStruct(gid) = ~isequal(length(tempVTRAILSc),length(procVTRAILSc)); % A Little Rough, but...
    
    VTRAILS = cat(2, VTRAILS, procVTRAILSc);
    
end

ChangeOccurred = logical(sum(ChangeOccurred & VerifyStruct));

function [VTRAILS,ChangeOccurred] = VTR3D_refineByWeldingCloseNodes(VTRAILSc,vicinityTHRnodes,Spacing)

VTRAILS = struct([]);

GroupIDs = unique([VTRAILSc(:).GroupID]);

ChangeOccurred = false(1,length(GroupIDs));
VerifyStruct = false(1,length(GroupIDs));

% Processing Each Tree (in case of a Forest) according to the GroupID
for gid = 1 : length(GroupIDs)
    
    tempVTRAILSc = VTRAILSc([VTRAILSc(:).GroupID] == GroupIDs(gid));
    
    [procVTRAILSc,ChangeOccurred(gid)] = VTR3D_refineNodesEachTree(tempVTRAILSc,vicinityTHRnodes,Spacing);
    
    VerifyStruct(gid) = ~isequal(length(tempVTRAILSc),length(procVTRAILSc)); % A Little Rough, but...
    
    VTRAILS = cat(2, VTRAILS, procVTRAILSc);
    
end

ChangeOccurred = logical(sum(ChangeOccurred & VerifyStruct));

function [VTRAILSTree,ChangeOccurred] = VTR3D_refineEdgesEachTree(VTRAILSc,vicinityTHRleaves,Spacing)

VTRAILSTree = struct([]);

for jj = 1 : length(VTRAILSc)
   VTRAILSc(jj).hasChanged = false; 
end

% Remove Single-Point Sequences
[VTRAILSc,RemovalFlag] = VTR3D_removeSingleOrIdenticalPointSequences(VTRAILSc);

WeldFlag = false(1,length(VTRAILSc));

% Process for Edge Welding
for el = 1 : length(VTRAILSc)
   
    if ~isempty(VTRAILSc(el).CGPathContinuous)
        [VTRAILSc,WeldFlag(el)] = VTR3D_refineByEdgeWelding(VTRAILSc, el , vicinityTHRleaves, Spacing );
    end
    
end

% Removing Empty Elements
for el = 1 : length(VTRAILSc)
    if ~isempty(VTRAILSc(el).CGPathContinuous)
        VTRAILSTree = cat(2,VTRAILSTree,VTRAILSc(el));
    end
end

% Determine if ChangeOccurred
ChangeOccurred = RemovalFlag | logical(sum(WeldFlag));

function [VTRAILSTree,ChangeOccurred] = VTR3D_refineNodesEachTree(VTRAILSc,vicinityTHRnodes,Spacing)

VTRAILSTree = struct([]);

for jj = 1 : length(VTRAILSc)
   VTRAILSc(jj).hasChanged = false; 
end

% Remove Single-Point Sequences
[VTRAILSc,RemovalFlag] = VTR3D_removeSingleOrIdenticalPointSequences(VTRAILSc);

%%%This part has been introduced on April 18th 2018
if RemovalFlag
    % Remove Empty entries
    VTRAILSc = VTRAILSc([VTRAILSc(:).EuclideanLength] > 0);
    % Recompute Brute-Force All AdjacencyList
    VTRAILSc = VTR3D_UpdateBruteForceTreeAdjacencyList(VTRAILSc);
end
%%%

WeldFlag = false(1,length(VTRAILSc));

% Process for Node welding
for el = 1 : length(VTRAILSc)
   
    if (~isempty(VTRAILSc(el).CGPathContinuous)) && (~VTRAILSc(el).isLeaf)
        [VTRAILSc,WeldFlag(el)] = VTR3D_refineByNodeWelding(VTRAILSc, el , vicinityTHRnodes, Spacing);
    end
    
end

% Removing Empty Elements
for el = 1 : length(VTRAILSc)
    if ~isempty(VTRAILSc(el).CGPathContinuous)
        VTRAILSTree = cat(2,VTRAILSTree,VTRAILSc(el));
    end
end

% Determine if ChangeOccurred
ChangeOccurred = RemovalFlag | logical(sum(WeldFlag));

function VTRAILSc = VTR3D_UpdateBruteForceTreeAdjacencyList(VTRAILSc)

for ii = 1 : size(VTRAILSc,2)

    AdjLst = [];
    
    for jj = 1 : size(VTRAILSc,2)
        
        if ii~=jj
            
            rep = VTR3D_repeatedValue_mex(VTRAILSc(ii).CGPathContinuous,VTRAILSc(jj).CGPathContinuous);
            
            if ~isempty(rep)
                
                assert(size(rep,1) == 1,'More than one Point In Common: There might be a Loop in the Graph(?)');
                
                AdjLst = cat(1,AdjLst,jj);
                
            end
            
        end
        
    end
    
    if isempty(AdjLst) && ~isempty(VTRAILSc(ii).CGPathAdjacencyList)
        warning('Initial Adjacency List is NOT empty -- will be replaced by VOID Adjacency List! -- Disconnected Element! (Bug?)');
    end 
    
    VTRAILSc(ii).CGPathAdjacencyList = AdjLst;
    
end

function [VTRAILSc,WeldFlag] = VTR3D_refineByEdgeWelding(VTRAILSc, idx, vicinityTHR, Spacing)

WeldFlag = false;

% INITIALISATION
SeqA = VTRAILSc(idx);
Seqs2Evaluate = SeqA.CGPathAdjacencyList;
SeqsVisited = [];
ChangedWith = [];

while ( ~isempty( Seqs2Evaluate ) && ~isequal( sort(Seqs2Evaluate) , sort(SeqsVisited) ) )
    
    [~,chk] = VTR3D_notUnique_mex(Seqs2Evaluate,SeqsVisited);
    elm2consider = Seqs2Evaluate( find(~chk,1,'first') );
    
    SeqB = VTRAILSc(elm2consider);
    if ~isempty(SeqB.CGPathContinuous)
        % Aligning Sequences
        [SeqA,SeqB,skip] = VTR3D_alignSequences(SeqA,SeqB);
        
        if ~skip
            % Compute the ICP between the SeqA and SeqB AND get the EdgeWeldFlag
            [ICP,EdgeWeldFlag] = VTR3D_evaulateICPwithvicinityTHR(SeqA,SeqB,vicinityTHR);
            
            newChange = isempty(VTR3D_notUnique_mex(elm2consider,ChangedWith));
            
            if EdgeWeldFlag && newChange % There is an Edge to Weld!!!
                
                % WeldFlag
                WeldFlag = true;
                
                % Process SeqA and SeqB and Update Data in VTRAILSc
                VTRAILSc = VTR3D_processEdgeWelding(VTRAILSc,SeqA,SeqB,ICP,Spacing);
                
                % Store the element that produced the CHANGE
                ChangedWith = cat(1,ChangedWith,elm2consider);
                % NB: an entry cannot change more than once with the same
                % considered element!
                
                % Reset: Seqs2Evaluate   AND   SeqsVisited
                SeqA = VTRAILSc(idx);
                Seqs2Evaluate = SeqA.CGPathAdjacencyList;
                SeqsVisited = [];
                
            else % No edges to Weld!
                %Updating
                SeqsVisited = cat(1,SeqsVisited,elm2consider);
            end
        else
            SeqsVisited = cat(1,SeqsVisited,elm2consider);
        end
        
    else
        SeqsVisited = cat(1,SeqsVisited,elm2consider);
        warning('Welding Paths: a sequence entry is empty! -- IGNORED');
    end
 
end

function [SeqA,SeqB,skip] = VTR3D_alignSequences(SeqA,SeqB)

skip = false;

% Alignment so that Both SeqA and SeqB start with the same shared and
% connecting point 3D
[~,Acheck,Bcheck] = VTR3D_repeatedValue_mex(SeqA.CGPathContinuous,SeqB.CGPathContinuous);

posA = find(Acheck);
posB = find(Bcheck);

if isempty(posA) || isempty(posB)
    skip = true;
else
    
    if posA == 1 && posB == 1
        % Do Nothing
    elseif posA == 1 && posB == length(Bcheck)
        SeqB.CGPathContinuous = flipud(SeqB.CGPathContinuous);
        SeqB.GeodesicEnergy = flipud(SeqB.GeodesicEnergy);
    elseif posA == length(Acheck) && posB == 1
        SeqA.CGPathContinuous = flipud(SeqA.CGPathContinuous);
        SeqA.GeodesicEnergy = flipud(SeqA.GeodesicEnergy);
    elseif posA == length(Acheck) && posB == length(Bcheck)
        SeqA.CGPathContinuous = flipud(SeqA.CGPathContinuous);
        SeqA.GeodesicEnergy = flipud(SeqA.GeodesicEnergy);
        
        SeqB.CGPathContinuous = flipud(SeqB.CGPathContinuous);
        SeqB.GeodesicEnergy = flipud(SeqB.GeodesicEnergy);
    else
        assert(sum(Acheck) == 1 && sum(Bcheck) == 1,'Bug: Multiple identical values within the 3D point sequence!');
        assert(SeqA.hasChanged || SeqB.hasChanged,'Bug: Connecting Point in the middle of a Sequence, without having the Sequence been CHANGED!');
        skip = true;
    end
    
end

function [ICP,EdgeWeldFlag] = VTR3D_evaulateICPwithvicinityTHR(SeqA,SeqB,vicinityTHR)

% Compute the ICP between the SeqA and SeqB
[ICP.A2B,ICP.B2A,...
 ICP.A2Bidx,ICP.B2Aidx] = VTR3D_computeICPfor3DPointsSequences(SeqA.CGPathContinuous,SeqB.CGPathContinuous);

% Evaluate the ICPs with the threshold
% IMPORTANT: the sequence will be considered all suprathreshold from
% the first element to the last suprathreshold element!!!
% This means that sequences like:
%                       A2Bthr = [ 1 1 1 0 0 1 1 0 0 0 0 0 ]
% WILL BE CONVERTED TO:                  | |
%                       A2Bthr = [ 1 1 1 1 1 1 1 0 0 0 0 0 ]
%
% Note: the final sequence has NO GAP(s)!
A2Bthr = ICP.A2B <= vicinityTHR ;
ICP.A2Bthr = false(size(A2Bthr));
ICP.A2Bthr(1:find(A2Bthr,1,'last')) = true;

B2Athr = ICP.B2A <= vicinityTHR ;
ICP.B2Athr = false(size(B2Athr));
ICP.B2Athr(1:find(B2Athr,1,'last')) = true;

if sum(ICP.A2Bthr) > 1 || sum(ICP.B2Athr) > 1
    % There are portion of sequences that can be fused!
    EdgeWeldFlag = true;
else
    % No fusion between sequences (except for the only common connecting point)
    EdgeWeldFlag = false;
end

function [ICPa,ICPb,ICPaIDX,ICPbIDX] = VTR3D_computeICPfor3DPointsSequences(a,b)%
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

function VTRAILSc = VTR3D_processEdgeWelding(VTRAILSc,SeqA,SeqB,ICP,Spacing)

% Identification of the portion to fuse (which to be removed and which to be kept)
% Identification based on First: on isLeaf flag (preserve inner branches). 
% Then: on the Geodesic Energy!

Leaves = [SeqA.isLeaf,SeqB.isLeaf];

if ~xor(SeqA.isLeaf,SeqB.isLeaf) % if both leaves or both inner branches
    GEs = [sum(SeqA.GeodesicEnergy(ICP.A2Bthr)),sum(SeqB.GeodesicEnergy(ICP.B2Athr))];
    [~,toChange] = max(GEs);
else
    toChange = find(Leaves);
end

isLeaf = Leaves(toChange);

if toChange == 1
    WeldedSeq = SeqA;
    UntouchSeq = SeqB;
    
    ICPsel.W2U = ICP.A2B;
    ICPsel.W2Uidx = ICP.A2Bidx;
    ICPsel.W2Uthr = ICP.A2Bthr;
else
    WeldedSeq = SeqB;
    UntouchSeq = SeqA;
    
    ICPsel.W2U = ICP.B2A;
    ICPsel.W2Uidx = ICP.B2Aidx;
    ICPsel.W2Uthr = ICP.B2Athr;
end

switch isLeaf
    
    case true
        % Processing for Welding Leaves
        % Change CGPathContinuous Sequence
        WeldedSeq.CGPathContinuous = cat( 1 ,...
                                          UntouchSeq.CGPathContinuous( ICPsel.W2Uidx( find( ICPsel.W2Uthr , 1 , 'last' ) ) , : ), ...
                                          WeldedSeq.CGPathContinuous( ~(ICPsel.W2Uthr) , : ) );

        % Check HERE the consistency of the WeldedSequence:
        % if the WeldedSeq.CGPathContinuous has only one 3D Point,
        % then it means that it has been WELDED COMPLETELY with the
        % UntouchSeq: therefore WeldedSeq has to be 'removed' by RESET.
        % AKA: Replaced by empty, as well as every field.
        if size(WeldedSeq.CGPathContinuous,1) > 1
            
            % Change GeodesicEnergy Sequence
            WeldedSeq.GeodesicEnergy   = cat( 1 ,...
                                              UntouchSeq.GeodesicEnergy( ICPsel.W2Uidx( find( ICPsel.W2Uthr , 1 , 'last' ) ) , : ), ...
                                              WeldedSeq.GeodesicEnergy( ~(ICPsel.W2Uthr) , : ) );
            % Update the EuclideanLength
            WeldedSeq.EuclideanLength = VTR3D_computeEuclideanLength(WeldedSeq.CGPathContinuous,Spacing);
        
            % Update the OWN CGPathAdjacencyList
            OldAdjacencyList = WeldedSeq.CGPathAdjacencyList;
            WeldedSeq.CGPathAdjacencyList = VTR3D_updateAdjacencyList4EdgeWelding(VTRAILSc,WeldedSeq.CGPathContinuous,OldAdjacencyList);
            
        else
            
            % Reset to Void Entry!
            WeldedSeq.CGPathContinuous = [];
            WeldedSeq.GeodesicEnergy = [];
            WeldedSeq.EuclideanLength = 0;
            % Update the OWN CGPathAdjacencyList
            OldAdjacencyList = WeldedSeq.CGPathAdjacencyList;
            WeldedSeq.CGPathAdjacencyList = [];
            
        end
        
        % Flag: hasCHANGED
        WeldedSeq.hasChanged = true;
        
        % Replace newSeqA into VTRAILSc HERE:
        VTRAILSc(WeldedSeq.CGPathID) = WeldedSeq;
        
        % Update the OTHERS CGPathAdjacencyList
        VTRAILSc = VTR3D_updateOthersAdjacencyList(VTRAILSc,OldAdjacencyList);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % End Processing For Leaves %
        
    case false
        % Processing for Welding Connecting Branches
        
        WeldedSeqBKP = WeldedSeq;
        
        % Change CGPathContinuous Sequence
        WeldedSeq.CGPathContinuous = cat( 1 ,...
                                          UntouchSeq.CGPathContinuous( ICPsel.W2Uidx( find( ICPsel.W2Uthr , 1 , 'last' ) ) , : ), ...
                                          WeldedSeq.CGPathContinuous( ~(ICPsel.W2Uthr) , : ) );
                                      
        % Check HERE the consistency of the WeldedSequence:
        % if the WeldedSeq.CGPathContinuous has only one 3D Point,
        % then it means that it has been WELDED COMPLETELY with the
        % UntouchSeq: therefore WeldedSeq has to be 'removed' by RESET.
        % AKA: Replaced by empty, as well as every field.
        if size(WeldedSeq.CGPathContinuous,1) > 1
            
            % Change GeodesicEnergy Sequence
            WeldedSeq.GeodesicEnergy   = cat( 1 ,...
                                              UntouchSeq.GeodesicEnergy( ICPsel.W2Uidx( find( ICPsel.W2Uthr , 1 , 'last' ) ) , : ), ...
                                              WeldedSeq.GeodesicEnergy( ~(ICPsel.W2Uthr) , : ) );
            % Update the EuclideanLength
            WeldedSeq.EuclideanLength = VTR3D_computeEuclideanLength(WeldedSeq.CGPathContinuous,Spacing);
        
            % Update the OWN CGPathAdjacencyList
            OldAdjacencyList = WeldedSeq.CGPathAdjacencyList;
            WeldedSeq.CGPathAdjacencyList = VTR3D_updateAdjacencyList4EdgeWelding(VTRAILSc,WeldedSeq.CGPathContinuous,OldAdjacencyList);
            
        else
            
            % Reset to Void Entry!
            WeldedSeq.CGPathContinuous = [];
            WeldedSeq.GeodesicEnergy = [];
            WeldedSeq.EuclideanLength = 0;
            % Update the OWN CGPathAdjacencyList
            OldAdjacencyList = WeldedSeq.CGPathAdjacencyList;
            WeldedSeq.CGPathAdjacencyList = [];
            
        end
        
        % Working on the Welded Portion
        VTRAILSc = VTR3D_updateBranchesWithinWeldedPortion(VTRAILSc,WeldedSeqBKP,UntouchSeq,ICPsel,Spacing);
        
        % Flag: hasCHANGED
        WeldedSeq.hasChanged = true;
        
        % UPDATE:
        % Replace WeldedSeq into VTRAILSc HERE:
        VTRAILSc(WeldedSeq.CGPathID) = WeldedSeq;
        
        % Update the OTHERS CGPathAdjacencyList
        VTRAILSc = VTR3D_updateOthersAdjacencyList(VTRAILSc,OldAdjacencyList);
        
    otherwise
        error('isLeaf Flag not Defined!');
end

function VTRAILSc = VTR3D_updateBranchesWithinWeldedPortion(VTRAILSc,WeldedSeqBKP,UntouchSeq,ICPsel,Spacing)

WPortion = WeldedSeqBKP.CGPathContinuous( (ICPsel.W2Uthr) , : ); %portion that is removed because of welding
UPortion = UntouchSeq.CGPathContinuous( 1 : ICPsel.W2Uidx( find( ICPsel.W2Uthr , 1 , 'last' ) ) , : );
geoUPortion = UntouchSeq.GeodesicEnergy( 1 : ICPsel.W2Uidx( find( ICPsel.W2Uthr , 1 , 'last' ) ) , : );

OldAdjacencyList = WeldedSeqBKP.CGPathAdjacencyList;

% Scan all the Connected Branches that HAVE shared a Connecting
% Point with the Welded Portion
[ListOfSeqs2Consider,ConnectingPoints] = VTR3D_scanConnectedSeqsWithinWeldedPortion(VTRAILSc,WPortion,OldAdjacencyList);

% Scan which Branches in the ListOfSeqs2Consider DO NOT HAVE a shared
% point with the Untouched Portion: UPortion
% ONLY these MUST be processed!
[ListOfSeqs2Update,IniPoint2StartConcatenating,ClosestPoints2Concatenate,Geodesics2Concatenate] = VTR3D_scanConsideredSeqs2BeUpdated(ListOfSeqs2Consider,UPortion,ConnectingPoints,geoUPortion);

% Update and Substitute Values in VTRAILSc so that the Connected Branches
% WITHIN the Welded Portion are now Connected with the Untouched Sequence
VTRAILSc = VTR3D_updateConnectedBranchesInStruct(VTRAILSc,ListOfSeqs2Update,IniPoint2StartConcatenating,ClosestPoints2Concatenate,Geodesics2Concatenate,Spacing,UntouchSeq);

function [ListOfSeqs2Consider,ConnectingPoints] = VTR3D_scanConnectedSeqsWithinWeldedPortion(VTRAILSc,WPortion,OldAdjacencyList)

ListOfSeqs2Consider = [];
ConnectingPoints = [];

for ss = 1 : length(OldAdjacencyList)
    
    seq2chk = VTRAILSc(OldAdjacencyList(ss)).CGPathContinuous;
    
    if ~isempty(WPortion) && ~isempty(seq2chk)
        
        [repVals,~,chk] = VTR3D_repeatedValue_mex(WPortion,seq2chk);
        
        if ~isempty(repVals)
            
            if size(repVals,1) ~= 1
                warning('More than one Point shared between Sequences: LOOP!');
            end
            
            if ~ ( find(chk) == 1 || find(chk) == length(chk) )
                warning('Shared Point Found in the Middle of a Sequence!');
            end

            ListOfSeqs2Consider = cat(1,ListOfSeqs2Consider,OldAdjacencyList(ss));
            ConnectingPoints = cat(1,ConnectingPoints,repVals);
            
        end
        
    end
    
end

function [ListOfSeqs2Update,IniPoint2StartConcatenating,ClosestPoints2Concatenate,Geodesics2Concatenate] = VTR3D_scanConsideredSeqs2BeUpdated(ListOfSeqs2Consider,UPortion,ConnectingPoints,geoUPortion)
ListOfSeqs2Update = [];
IniPoint2StartConcatenating = [];
ClosestPoints2Concatenate = [];
Geodesics2Concatenate = [];

for ss = 1 : length(ListOfSeqs2Consider)
    
    [repVals,~] = VTR3D_repeatedValue_mex(UPortion,ConnectingPoints(ss,:));
    
    if isempty(repVals)
        
        ListOfSeqs2Update = cat(1,ListOfSeqs2Update,ListOfSeqs2Consider(ss));
        IniPoint2StartConcatenating = cat( 1, IniPoint2StartConcatenating, ConnectingPoints(ss,:) );
        
        [~,~,UCPidx,~] = VTR3D_computeICPfor3DPointsSequences(ConnectingPoints(ss,:), UPortion);
        UClosestConnectingPoint = UPortion(UCPidx,:);
        
        ClosestPoints2Concatenate = cat( 1, ClosestPoints2Concatenate, UClosestConnectingPoint );
        Geodesics2Concatenate = cat( 1, Geodesics2Concatenate, geoUPortion(UCPidx) );
    end
    
end

function VTRAILSc = VTR3D_updateConnectedBranchesInStruct(VTRAILSc,ListOfSeqs2Update,IniPoint2StartConcatenating,ClosestPoints2Concatenate,Geodesics2Concatenate,Spacing,UntouchSeq)
% Updating the Connected Branches:
if ~isempty(ListOfSeqs2Update)
    for elm = 1 : length(ListOfSeqs2Update)
        
        tempSeq = VTRAILSc(ListOfSeqs2Update(elm)).CGPathContinuous;
        tempGeo = VTRAILSc(ListOfSeqs2Update(elm)).GeodesicEnergy;
        
        [~,chk,~] = VTR3D_repeatedValue_mex(tempSeq,IniPoint2StartConcatenating(elm,:));
        if find(chk) == 1
            tempSeq = cat( 1, ClosestPoints2Concatenate(elm,:) , tempSeq(~chk,:) );
            tempGeo = cat( 1, Geodesics2Concatenate(elm), tempGeo(~chk) );
        elseif find(chk) == length(chk)
            tempSeq = cat( 1, tempSeq(~chk,:) , ClosestPoints2Concatenate(elm,:) );
            tempGeo = cat( 1, tempGeo(~chk) , Geodesics2Concatenate(elm));
        else
            error('In between Connecting Point! -- BUG');
        end
        
        % Updating Values: Sequence and Geodesics
        VTRAILSc(ListOfSeqs2Update(elm)).CGPathContinuous = tempSeq;
        VTRAILSc(ListOfSeqs2Update(elm)).GeodesicEnergy = tempGeo;
        % Updating the Euclidean Length
        VTRAILSc(ListOfSeqs2Update(elm)).EuclideanLength = VTR3D_computeEuclideanLength(VTRAILSc(ListOfSeqs2Update(elm)).CGPathContinuous,Spacing);
        % Updating the AdjacencyList
        VTRAILSc(ListOfSeqs2Update(elm)).CGPathAdjacencyList = unique( cat( 1, VTRAILSc(ListOfSeqs2Update(elm)).CGPathAdjacencyList, UntouchSeq.CGPathID ) );
        
        VTRAILSc(ListOfSeqs2Update(elm)).hasChanged = true;
        
    end
end

function NewAdjacencyList = VTR3D_updateAdjacencyList4EdgeWelding(VTRAILSc,Sequence1,OldAdjacencyList)

NewAdjacencyList = [];

for elm = 1 : length(OldAdjacencyList)
   
    Sequence2 = VTRAILSc(OldAdjacencyList(elm)).CGPathContinuous;
    
    if ~ isempty(Sequence2)
        
        [repVals,~,~] = VTR3D_repeatedValue_mex(Sequence1,Sequence2);
        
        if ~isempty(repVals)
            
            assert(size(repVals,1) == 1,'More than one Point in Common: Possible LOOP!');
            
            NewAdjacencyList = cat(1,NewAdjacencyList,OldAdjacencyList(elm));
            
        end
    end
    
end

function VTRAILSc = VTR3D_updateOthersAdjacencyList(VTRAILSc,OldAdjacencyList)

for jj = 1 : length(OldAdjacencyList)
    
    Local_OldAdjacencyList = VTRAILSc(OldAdjacencyList(jj)).CGPathAdjacencyList;
    Local_Sequence = VTRAILSc(OldAdjacencyList(jj)).CGPathContinuous;
    
    Local_NewAdjacencyList = VTR3D_updateAdjacencyList4EdgeWelding(VTRAILSc,Local_Sequence,Local_OldAdjacencyList);
    VTRAILSc(OldAdjacencyList(jj)).CGPathAdjacencyList = Local_NewAdjacencyList;
    
    
    clear Local_OldAdjacencyList Local_Sequence Local_NewAdjacencyList;
end

function [VTRAILSc,RemovalFlag] = VTR3D_removeSingleOrIdenticalPointSequences(VTRAILSc)

RemovalFlag = false;

for jj = 1 : length(VTRAILSc)
   
    if ~isempty(VTRAILSc(jj).CGPathContinuous)
        
        if size(VTRAILSc(jj).CGPathContinuous,1) < 2 || sum(sqrt(sum(diff(VTRAILSc(jj).CGPathContinuous,1).^2,2))) < 1e-6
        
            RemovalFlag = true;
            
            % Reset to Void Entry!
            VTRAILSc(jj).CGPathContinuous = [];
            VTRAILSc(jj).GeodesicEnergy = [];
            VTRAILSc(jj).EuclideanLength = 0;
            
            % Update the OWN CGPathAdjacencyList
            OldAdjacencyList = VTRAILSc(jj).CGPathAdjacencyList;
            VTRAILSc(jj).CGPathAdjacencyList = [];
            
            % Flag: hasCHANGED
            VTRAILSc(jj).hasChanged = true;
            
            % Update the OTHERS CGPathAdjacencyList
            VTRAILSc = VTR3D_updateOthersAdjacencyList(VTRAILSc,OldAdjacencyList);
            
        end
    end
    
end

function VTRAILSc_out = VTR3D_removeDegenerateEntries(VTRAILSc)

VTRAILSc_out = struct([]);

GIDs = unique([VTRAILSc(:).GroupID]);

for gid = 1 : length(GIDs)
    
    VTRAILStemp = VTRAILSc([VTRAILSc(:).GroupID] == GIDs(gid));
    
    for jj = 1 : length(VTRAILStemp)
        
        if ~isempty(VTRAILStemp(jj).CGPathContinuous)
            
            if size(VTRAILStemp(jj).CGPathContinuous,1) < 2 || sum(sqrt(sum(diff(VTRAILStemp(jj).CGPathContinuous,1).^2,2))) < 1e-6
                
                % Reset to Void Entry!
                VTRAILStemp(jj).CGPathContinuous = [];
                VTRAILStemp(jj).GeodesicEnergy = [];
                VTRAILStemp(jj).EuclideanLength = 0;
                
            end
        end
        
    end

    % Remove Void Entries (from the above step)
    VTRAILStemp = VTRAILStemp([VTRAILStemp(:).EuclideanLength] > 0);
    % Re-Enstablish GCPathID
    for jj = 1 : length(VTRAILStemp)
        VTRAILStemp(jj).CGPathID = jj;
    end
    % Recompute Brute-Force All AdjacencyList
    VTRAILStemp = VTR3D_UpdateBruteForceTreeAdjacencyList(VTRAILStemp);
    
    % Concatenate Results
    VTRAILSc_out = cat(2, VTRAILSc_out , VTRAILStemp );

end

function OldValueStr = VTR3D_showProgressCM(ProgressMessage,OldValueStr,Value,Finished)

if isempty(ProgressMessage)
    ProgressMessage = '';
end

if ~Finished
    NewValueStr = sprintf(' %3.2f', Value);
    fprintf([OldValueStr,ProgressMessage,NewValueStr]);
    OldValueStr = repmat(sprintf('\b'), 1, length(NewValueStr));
else
    NewValueStr = sprintf(' %3.2f', Value);
    fprintf([OldValueStr,ProgressMessage,NewValueStr,' -- DONE\n']);
    if strcmp(ProgressMessage,'')
        OldValueStr = repmat(sprintf('\b'), 1, length([ProgressMessage,NewValueStr,' -- DONE ']));
    else
        OldValueStr = repmat(sprintf('\b'), 1, length([OldValueStr,ProgressMessage,NewValueStr,' -- DONE ']));
    end
    pause(0.1);
end

function [VTRAILSc,WeldFlag] = VTR3D_refineByNodeWelding(VTRAILSc, idx , vicinityTHRnodes, Spacing)

WeldFlag = false;

% INITIALISATION
SeqA = VTRAILSc(idx);

if SeqA.EuclideanLength <= vicinityTHRnodes
    
    Seqs2Evaluate = SeqA.CGPathAdjacencyList;
    
    CGPathContinuous2Replace = mean(SeqA.CGPathContinuous,1);
    GeodesicEnergy2Replace = mean(SeqA.GeodesicEnergy);
    
    for adj = 1 : length(Seqs2Evaluate)
        
        SeqB = VTRAILSc(Seqs2Evaluate(adj));
        
        d_ini = sqrt( sum( ( CGPathContinuous2Replace - SeqB.CGPathContinuous( 1 ,:)  ).^2 ) );
        d_fin = sqrt( sum( ( CGPathContinuous2Replace - SeqB.CGPathContinuous(end,:)  ).^2 ) );
        
        if d_ini < d_fin
            % Add the CGPathContinuous2Replace at the BEGINNING of the Adjacent Sequence
            SeqB.CGPathContinuous = cat(1,CGPathContinuous2Replace,SeqB.CGPathContinuous(2:end,:));
            SeqB.GeodesicEnergy = cat(1,GeodesicEnergy2Replace,SeqB.GeodesicEnergy(2:end,:));
            SeqB.EuclideanLength = VTR3D_computeEuclideanLength(SeqB.CGPathContinuous,Spacing);
            SeqB.CGPathAdjacencyList = unique( cat(1,Seqs2Evaluate,SeqB.CGPathAdjacencyList) );
            SeqB.CGPathAdjacencyList = SeqB.CGPathAdjacencyList(SeqB.CGPathAdjacencyList ~= idx);
            SeqB.CGPathAdjacencyList = SeqB.CGPathAdjacencyList(SeqB.CGPathAdjacencyList ~= Seqs2Evaluate(adj));
            SeqB.hasChanged = true;
            
            % Replace after Update
            VTRAILSc(Seqs2Evaluate(adj)) = SeqB;
            
            WeldFlag(Seqs2Evaluate(adj)) = true;
            
        else % d_fin < d_ini
            % Add the CGPathContinuous2Replace at the END of the Adjacent Sequence
            SeqB.CGPathContinuous = cat(1,SeqB.CGPathContinuous(1:end-1,:),CGPathContinuous2Replace);
            SeqB.GeodesicEnergy = cat(1,SeqB.GeodesicEnergy(1:end-1,:),GeodesicEnergy2Replace);
            SeqB.EuclideanLength = VTR3D_computeEuclideanLength(SeqB.CGPathContinuous,Spacing);
            SeqB.CGPathAdjacencyList = unique( cat(1,Seqs2Evaluate,SeqB.CGPathAdjacencyList) );
            SeqB.CGPathAdjacencyList = SeqB.CGPathAdjacencyList(SeqB.CGPathAdjacencyList ~= idx);
            SeqB.CGPathAdjacencyList = SeqB.CGPathAdjacencyList(SeqB.CGPathAdjacencyList ~= Seqs2Evaluate(adj));
            SeqB.hasChanged = true;
            
            % Replace after Update
            VTRAILSc(Seqs2Evaluate(adj)) = SeqB;
            
        end
        
        clear seqB d_ini d_fin;
        
    end
    
    % Reset to void the considered ENTRY idx!
    SeqA.CGPathContinuous = [];
    SeqA.GeodesicEnergy = [];
    SeqA.EuclideanLength = 0;
    SeqA.CGPathAdjacencyList = [];
    SeqA.hasChanged = true;
    
    VTRAILSc(idx) = SeqA;
    
    WeldFlag = true;
    
end


%%% ROI-based iterative Pruning

function isInPruningROI = getLeavesInPruningROI(selVTRAILS,LeavesPruningROI,thrhld)

isInPruningROI = false(size([selVTRAILS(:).isLeaf]));

LeavesPruningROIperm = double(permute( LeavesPruningROI , [3 1 2] ));

xM = 1:size(LeavesPruningROI,2);
yM = 1:size(LeavesPruningROI,1);
zM = 1:size(LeavesPruningROI,3);

for jj = 1 : size(selVTRAILS,2)
    if selVTRAILS(jj).isLeaf
        vals = linterp_matlab( xM, yM, zM , LeavesPruningROIperm , selVTRAILS(jj).CGPathContinuous(:,1) , selVTRAILS(jj).CGPathContinuous(:,2) , selVTRAILS(jj).CGPathContinuous(:,3) );
        dd = sqrt(sum(diff(selVTRAILS(jj).CGPathContinuous,1,1).^2,2));
        M = sum(vals(2:end)'.*dd)/sum(dd);
        if M > thrhld
            isInPruningROI(jj) = true;
        end
    end
end

function [VTRAILSc,RemovalFlag] = VTR3D_removePossibleSinglePoint(VTRAILSc,SET)

RemovalFlag = false;

if ~isempty(SET.LeavesPruningROI)
    
    for jj = 1 : length(VTRAILSc)
        
        if ~isempty(VTRAILSc(jj).CGPathContinuous)
            
            if size(VTRAILSc(jj).CGPathContinuous,1) < 2 || sum(sqrt(sum(diff(VTRAILSc(jj).CGPathContinuous,1).^2,2))) < 1e-6 || isempty(VTRAILSc(jj).CGPathAdjacencyList)
                
                RemovalFlag = true;
                
                % Reset to Void Entry!
                VTRAILSc(jj).CGPathContinuous = [];
                VTRAILSc(jj).GeodesicEnergy = [];
                VTRAILSc(jj).EuclideanLength = 0;

            end
        end
        
    end
    
    if RemovalFlag
        VTRAILSc = VTRAILSc([VTRAILSc(:).EuclideanLength] > 0);
    end
    
end

function VTRAILSprun = removeFuzzyLeavesInPruningROI(VTRAILSc,SET)

if ~isempty(SET.LeavesPruningROI)
    
    GroupIDs2consider = VTR3D_retireveGIDs2Consider(VTRAILSc,SET);
    
    disp('_______________________________________________________');
    disp(' * Removing Extra Leaves from GeodesicGraph Masking ...');
    
    if SET.PruningIterativeSteps > 0
        
        tempVTRAILSc = VTRAILSc;
        
        converged = false;
        maxitr = SET.PruningIterativeSteps;
        prunitr = 0;
        
        while ~converged
            
            % Assign The Tree Root (ot Forest Roots) so that it cannot be
            % pruned!
            tempVTRAILSc = VTR3D_assignTreeRoot(tempVTRAILSc,SET.TreeRootSeed);
            
            % Pruning
            VTRAILScPruned = VTR3D_pruneVTrailsLeavesOnROI( tempVTRAILSc , GroupIDs2consider , SET );
            
            % Recovening back the fully connected TREE
            disp('-------------------------------------------------------');
            tempVTRAILScPruned = VTR3D_computeGeodesicConnectingPaths(true([length(VTRAILScPruned),1]),VTRAILScPruned,SET);
            
            prunitr = prunitr + 1;
            
            if (prunitr > maxitr) || (length(tempVTRAILSc) == length(tempVTRAILScPruned))
                VTRAILSprun = tempVTRAILScPruned;
                converged = true;
            else
                tempVTRAILSc = tempVTRAILScPruned;
                converged = false;
            end
            
        end
        
        clear tempVTRAILSc VTRAILScPruned tempPATHSc tempBRANCHESc;
        
    else
        VTRAILScPruned = [];
        
        for gid = 1 : length(GroupIDs2consider)
            
            VTRAILStemp = VTRAILSc([VTRAILSc(:).GroupID] == GroupIDs2consider(gid));
            
            VTRAILScPruned = cat(2, VTRAILScPruned , VTRAILStemp );
        end
        
        % Recovening back the fully connected TREE
        disp('-------------------------------------------------------');
        VTRAILScPruned = VTR3D_computeGeodesicConnectingPaths(true([length(VTRAILScPruned),1]),VTRAILScPruned,SET);
        
        VTRAILSprun = VTRAILScPruned;
    end
else
    VTRAILSprun = VTRAILSc;
end

function prunedVTRAILS = VTR3D_pruneVTrailsLeavesOnROI(VTRAILS,GroupID,SET)

prunedVTRAILS = struct([]);

for gID = 1 : length(GroupID)
   
    selVTRAILS = VTRAILS([VTRAILS(:).GroupID] == GroupID(gID));
    
    Leaves = [selVTRAILS(:).isLeaf];
    Root = [selVTRAILS(:).isRoot];
    
    isInPruningROI = getLeavesInPruningROI(selVTRAILS,SET.LeavesPruningROI,0.1);

    ExclusionCriterion = ( Leaves(:) & ~Root(:) & isInPruningROI(:) );

    selVTRAILS = selVTRAILS(~ExclusionCriterion);

    % Concatenation
    prunedVTRAILS = cat(2,prunedVTRAILS,selVTRAILS);
    
end

function includeVTrailsLibs

addpath(strcat(VTRootDir,'libs/'));
addpath(strcat(VTRootDir,'libs/0_NIfTI_IO/'));
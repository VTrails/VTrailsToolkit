function DATA = VTrailsConnectingGeodesicGraph3D_MAIN(MAP,varargin)%
%
%VTrails: Geodesic Connectivity Paradigm for Over-Connected Vascular Graphs
%
% Use this function to infer the over-connected vascular topology (graph)
% by running an Anisotropic Level-Set (Fast Marching) over the Riemannian
% Vesselness potential (CVM and TF) and by enforcing a self-organising
% exhaustive connectivity paradigm.
%
%%%%%%%%%
% Inputs:
%
% - MAP: structure containing Image-related features, the Riemannian
%        Vesselnes potential and the initial Seeds for the Level-Set.
%        NB: MAP can be either a single-scale output from the function:
%        'VTrailsFilter3D_MAIN.m' or it can be a multi-scale integrated
%        filter-response (MaxMAP) as from the function:
%        'VTF3D_IntegrateMaxMultiScaleResponses.m'.
%        
%                    MAP *MUST* contain the following Fields:
%
%  * IMG: Header of the input nifty Images.
%  * GRID: Structure with Image Grid and Scales' infos for Resampling.
%  * CVM: Integral Scalar Connected Vesselness Map obtained with SLoGS
%  * TF: Integral Vascular Tensors Field synthesized with SLoGS in the
%        Euclidean Domain -- can be [empty] for automatically FITTING
%        the Tensor Field over the CVM, or it can be a scalar equal to 1,
%        for the ISOTROPIC Front Propagation.
%  * ConnectedSegmentsOS: struct containing the list of SORTED Connected
%                         Components from the Organised Seeds.
%
%                            Other [optional] Fields:
%
%  - [CVMsmth]: As CVM, but slightly smoothed for local maxima.
%  - [BDM]: Integral Scalar Vessel Boundaries Map obtained with \deltaSLoGS
%  - [BGM]: Integral Scalar Vessel Background Map obtained with \nuSLoGS
%  - [TFLE]: Integral Vascular Tensors Field synthesized with SLoGS in
%            the Log-Euclidean Domain 
%  - [TFValidityMSK]: Binary Mask for non-isotropic Tensors.
%  - [ACC]: Internal Variable for Across-Scales Integration Regularisation
%  - [US]: Scalar Volume of Un-Organised Seeds TO BE Aligned to the
%          vessels with 'VTF3D_AlignVesselSeeds3D.m'
%  - [OS]: Scalar Volume of Organised Seeds, aligned to the
%          vessels with 'VTF3D_AlignVesselSeeds3D.m'
%  - [SortedOS]: Scalar Volume of Organised and SORTED Seeds, each
%                consisting in a labelled (numeric) path, obtained with
%                'VTF3D_SortConnectedComponents'.
%  - [BranchPoints]: list of early junctions detected with
%                    'VTF3D_SortConnectedComponents'.
%
% CONFIGURATION INPUTS:
%
% - ['Mask']: Validity Mask for Domain exploration (i.e. ~Barriers). -
%             Default: true(size(MAP.CVM))
% - ['SPDPowerFactor']: Power of the Vesselness Speed Potential (SPD). -
%                       Default: 1
% - ['VNRadiusMM']: Radius of Visible Seeds Neighborhood - limit the
%                   exploration of the domain to a restricted spherical
%                   neighborhood of the given Seed. - Default: 20 [mm]
% - ['NodeDistTHR']: Node Distance Threshold for Adaptive Graph Refinement.
%                    - Default: 3 [mm]
% - ['PathFolder']: Path Folder for Data Exporting - Default: [pwd,'/']
% - ['QExplore']: Quantile of Geodesic Exploration - to restrict the
%                 exploration of the domain accordingly with a geodesic
%                 metric. - Default: 0.5
% - ['FullDump']: Binary Scalar Flag for enabling FULL DUMP of all
%                 variables. - Default: false
%
%%%%%%%%%%
% Outputs:
%   
% - DATA: structure containing the Over-connected Vascular Graph.
%
%               DATA will contain the following Fields:
%
%  * GGM: Geodesically Weighted Connectivity Matrix of the Over-Connected
%         Vascular Graph.
%  * exploredminPath: struct containing the list of all the extracted and
%                     explored connecting geodesics (minimal paths).
%  * MST: Minimum Spanning Tree(s) Matrix obtained from GGM with
%         'Kruskal' method.
%  * OrigSeedsIDX: list of Original Seeds;
%  * ExtraSeedsIDX: list of automatically generated Seeds, during the
%                   self-arranging refinement; 
%  * IMG: as in MAP (input)
%  * GRID: as in MAP (input)
%
%               Other functional fields (for traceability):
%
%  * SPDPowerFactor: value of 'SPDPOWERFACTOR'.
%  * VisibleNeighborhoodRadius_mm: value of 'VNRADIUSMM'.
%  * NodeDistanceTHR_mm: value of 'NODEDISTTHR'.
%  * GradientDescentStopTHR: automatically determined  value.
%  * PathFolder: string with the Path Folder for Data Exporting.
%  * QExplore: value of 'QEXPLORE'.
%  * VisibleNeighborhoodRadius_vx: as in 'VNRADIUSMM', but VOXELs.
%  * VisibilityBall: Logical Volume created by 'VNRADIUSMM'.
%
%          [Other Fields are present when 'FULLDUMP' is enabled]
%
%%%%%%%%%%%
% Example Call:
%
%    DATA = VTrailsConnectingGeodesicGraph3D_MAIN( MaxMAP , 'PathFolder' , JobDumpDirPath , 'SPDPowerFactor' , 2);
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

DATA = VTC3D_configureINPUTS(MAP,varargin);

while DATA.NewNodesFLAG
    
    tic;
    
    DATA = computeAnisotropicGeodesics(DATA);
    
    DATA = determineGraphPairs(DATA);
    
    DATA = determineGeodesicGraph(DATA);
    
    toc;

    disp('---------------------------------------------------------');
end

DATA.GGM(isinf(DATA.GGM) | isnan(DATA.GGM)) = 0;

[~,CHKseeds] = VTC3D_notUnique(DATA.idx,DATA.OrigSeedsIDX);
DATA.ExtraSeedsIDX = DATA.idx(~CHKseeds);

[~] = VTC3D_exportGraphConnectivity(DATA);

disp('VTrails: Connecting Geodesic Graph -- COMPLETE!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DATA = VTC3D_configureINPUTS(MAP,optInputs)

disp('*********************************************************');
disp('  VTrails: Connecting Geodesic Graph  --  Vascular Tree  ');
disp('*********************************************************');

disp(strcat( datestr(datetime),' @ ', computer ));
disp('Data Initialisation: ...');

% Copying Fields
DATA.IMG = MAP.IMG;
DATA.GRID = MAP.GRID;
DATA.ConnectedSegmentsOS = MAP.ConnectedSegmentsOS;

% Managing Varargin
% Initializing OPT with Default Values
DATA.MSK = true(size(MAP.CVM));
DATA.SPDPowerFactor = 1;
DATA.VisibleNeighborhoodRadius_mm = 20;
DATA.NodeDistanceTHR_mm = 3;
DATA.GradientDescentStopTHR = 1.5;
DATA.PathFolder = [pwd,'/'];
DATA.QExplore = 0.5;
DATA.FullDump = false;

if ~isempty(optInputs)
    for j = 1 : 2 :length(optInputs)
        
        switch upper(optInputs{j})
            case 'MASK'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isequal(size(optInputs{j+1}),size(MAP.CVM))
                    disp('Invalid|Not Given Validity Mask: -- Default applied!');
                else
                    DATA.MSK = optInputs{j+1};
                end
            case 'SPDPOWERFACTOR'
                if isempty(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || ~isnumeric(optInputs{j+1})
                    disp('Invalid|Not Given SPDPowerFactor Value: -- Default: 1 applied!');
                else
                    DATA.SPDPowerFactor = abs(optInputs{j+1});
                end
            case 'VNRADIUSMM'
                if isempty(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || ~isnumeric(optInputs{j+1})
                    disp('Invalid|Not Given Visible Neighborhood Radius (mm): -- Default: 20 applied!');
                else
                    DATA.VisibleNeighborhoodRadius_mm = abs(optInputs{j+1});
                end
            case 'NODEDISTTHR'
                if isempty(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || ~isnumeric(optInputs{j+1})
                    disp('Invalid|Not Given Node Distance Threshold for Adaptive Graph Refinement (mm): -- Default: 3 applied!');
                else
                    DATA.NodeDistanceTHR_mm = abs(optInputs{j+1});
                end
            case 'PATHFOLDER'
                if isempty(optInputs{j+1})
                    disp(['Invalid|Not Given Path Folder for Data Exporting: -- Default: ', DATA.PathFolder ]);
                else
                    DATA.PathFolder = optInputs{j+1};
                end
            case 'QEXPLORE'
                if isempty(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || ~isnumeric(optInputs{j+1})
                    disp('Invalid|Not Given Quantile of Geodesic Exploration: -- Default: 0.5 applied!');
                else
                    DATA.QExplore = abs(optInputs{j+1});
                    DATA.QExplore(DATA.QExplore > 1) = 1;
                    DATA.QExplore(DATA.QExplore < 0) = 0.5;
                end
            case 'FULLDUMP'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given FullDump : -- Default: "false" applied!');
                else
                    DATA.FullDump = optInputs{j+1};
                end
            otherwise
                disp('Given Unknown parameter(s). DEFAULT value(s) will be set.');
        end
        
    end
end

%%%
% Defining the Scalar Speed Potential
%minMapVal = 1e-3; % Minimal tollerance value for the Anisotropic Fast Marching
minMapVal = 1e-2; % Minimal tollerance value for the Anisotropic Fast Marching
DATA.SPD = abs(double(MAP.CVM).^DATA.SPDPowerFactor) + minMapVal; % Speed Potential

if isempty(MAP.TF)
    wrnstr = warning('Tensor Field is Empty! -- Automatic Fitting from Input Data');
    
    % Fitting the Tensor on Data (If you are not using SLoGS).
    TensorValidTHRSHLD = quantile(MAP.CVM(MAP.CVM>0),0.75);
    TensorValidMSK = MAP.CVM > TensorValidTHRSHLD;
    % NB: all the other tensors associated with a Maximal Vesselness Map lower
    % than VmaxThreshold will be ISOTROPIC!
    
    % Fitting Process
    [TF,~] = VTC3D_get3DImageVoxelTensorLE( double(MAP.CVM) , DATA.MSK & TensorValidMSK , 1e-3 );
    
    % Remove the temporary Scale-Space Analysis Framework
    VTC3D_deleteSSAFW3D;

elseif ~isstruct(MAP.TF)
    wrnstr = warning('Unknown Tensor Field! -- Dummy Isotropic Tensor Field Generated');
    TF.El1 = ones(size(MAP.CVM));
    TF.El2 = ones(size(MAP.CVM));
    TF.El3 = ones(size(MAP.CVM));
    TF.Ev1 = cat(4,ones(size(MAP.CVM)),zeros(size(MAP.CVM)),zeros(size(MAP.CVM)));
    TF.Ev2 = cat(4,zeros(size(MAP.CVM)),ones(size(MAP.CVM)),zeros(size(MAP.CVM)));
    TF.Ev3 = cat(4,zeros(size(MAP.CVM)),zeros(size(MAP.CVM)),ones(size(MAP.CVM)));
else
    wrnstr = [];
    TF.El1 = double(MAP.TF.El1);
    TF.El2 = double(MAP.TF.El2);
    TF.El3 = double(MAP.TF.El3);
    TF.Ev1 = double(MAP.TF.Ev1);
    TF.Ev2 = double(MAP.TF.Ev2);
    TF.Ev3 = double(MAP.TF.Ev3);
end

DATA.VisibleNeighborhoodRadius_vx = round(DATA.VisibleNeighborhoodRadius_mm./(DATA.GRID.physicalvoxelsize ./ DATA.GRID.Scale));

[iVBx,iVBy,iVBz] = ndgrid( 1 : 2 * DATA.VisibleNeighborhoodRadius_vx(1) + 1 ,...
                           1 : 2 * DATA.VisibleNeighborhoodRadius_vx(2) + 1 ,...
                           1 : 2 * DATA.VisibleNeighborhoodRadius_vx(3) + 1);

DATA.VisibilityBall = false( (2 * DATA.VisibleNeighborhoodRadius_vx .* [1 1 1]) + 1);

DATA.VisibilityBall( ( ((iVBx - DATA.VisibleNeighborhoodRadius_vx(1) - 1).^2)/(DATA.VisibleNeighborhoodRadius_vx(1)^2) + ...
                       ((iVBy - DATA.VisibleNeighborhoodRadius_vx(2) - 1).^2)/(DATA.VisibleNeighborhoodRadius_vx(2)^2) + ...
                       ((iVBz - DATA.VisibleNeighborhoodRadius_vx(3) - 1).^2)/(DATA.VisibleNeighborhoodRadius_vx(3)^2) ) <= 1 ) = true;

if size(DATA.VisibilityBall,1) >= size(DATA.SPD,1) || size(DATA.VisibilityBall,2) >= size(DATA.SPD,2) || size(DATA.VisibilityBall,3) >= size(DATA.SPD,3) || DATA.QExplore >= 1
   DATA.VisibilityBall = true(size(DATA.SPD));                
end

% Converting the Tensor Field for the AFM
DATA.TF = VTC3D_TLICFromTEIG_3D(TF);

NV = isnan(DATA.TF(:,:,:,1))|isinf(DATA.TF(:,:,:,1)); % NotValid Mask
DATA.MSK = DATA.MSK & ~NV;

% Initialization of Processing Variables
DATA.GGM = []; % Geodesic Graph Matrix (Weighted Adjacency Matrix)
DATA.GM = struct([]); % Geodesic Map Struct computed from the AFM (based on a minimal bounding box of the ROI)
DATA.BB = struct([]); % Barriers Logical volume (should be similar to DATA.GM)
DATA.idx = []; % List of Seeds for the AFM Front Propagation
DATA.idxnew = [];
DATA.pairs = [];
DATA.exploredpairs = []; % Nx2 Matrix with the explored Pairs of Nodes (in order to avoid unnecessary computations)
DATA.minPath = [];
DATA.exploredminPath = struct([]); % Struct of already-visited Minimal Connecting Paths
DATA.territorypairs = [];

DATA.MST = []; % Minimum Spanning Tree

if ~isempty(DATA.ConnectedSegmentsOS)
    DATA.NewNodesFLAG = true;
else
    DATA.NewNodesFLAG = false;
end

if isempty(wrnstr)
    fprintf('\b\b\b\bDONE\n');
else
    disp('Data Initialisation: DONE');
end
disp(['Initial Connected Elements: ',num2str(length(DATA.ConnectedSegmentsOS.List),'%d')]);
disp('---------------------------------------------------------');

function DATA = determineGraphPairs(DATA)

BB = DATA.BB;
idx = DATA.idx;
exploredpairs = DATA.exploredpairs;

pairs = [];

LBL = zeros(BB(1).size);
LBL(idx) = (1:length(idx))';

ProgressValue = 0;
ProgressValue_step = 100.0/size(BB,2);
ProgressMessage = 'Determining Graph Pairs:';
ProgressOldString = VTC3D_showProgressCM(ProgressMessage,'',ProgressValue, false);

for jj = 1 : size(BB,2)
   
    BBtemp = VTC3D_retrieveOriginalMapFromCrop(BB,jj);
    
    list = nonzeros( LBL .* (~BBtemp) );
    
    pairs_new = [ jj*ones([length(list),1]) , list(:)];
    
    diagelm = pairs_new(:,1) == pairs_new(:,2);
    
    [~,pairs_newCHK,~] = VTC3D_repeatedValueSymm( pairs_new , pairs );
    
    pairs2cat = pairs_new( ~pairs_newCHK & ~diagelm , : );
    
    pairs = cat(1,pairs,pairs2cat);
    
    % Progress Percentage
    ProgressValue = ProgressValue + ProgressValue_step;
    ProgressOldString = VTC3D_showProgressCM([],ProgressOldString,ProgressValue, jj == size(BB,2));
    
end

[~,repeatedpairsCHK,~] = VTC3D_repeatedValueSymm(pairs,exploredpairs);
pairs = pairs(~repeatedpairsCHK,:);

DATA.pairs = pairs;

function DATA = computeAnisotropicGeodesics(DATA)

if isempty(DATA.GM)
    
    ProgressValue = 0;
    ProgressValue_step = 100.0/length(DATA.ConnectedSegmentsOS.List);
    ProgressMessage = 'Computing Geodesics:';
    ProgressOldString = VTC3D_showProgressCM(ProgressMessage,'',ProgressValue, false);
    
    for jj = 1 : length(DATA.ConnectedSegmentsOS.List)
        
        iniMSK = false(size(DATA.SPD));
        iniMSK(DATA.ConnectedSegmentsOS.List(jj).Sequence) = true;
        iniMSK = DATA.MSK & imdilate(iniMSK,DATA.VisibilityBall);
        
        if DATA.QExplore >= 1
           iniMSK = DATA.MSK; 
        end
        
        
        DATA.idxnew = unique([DATA.ConnectedSegmentsOS.List(jj).Sequence(1); DATA.ConnectedSegmentsOS.List(jj).Sequence(end)]);
        
        [GMnew,BBnew] = VTC3D_computeGeodesicMaps( DATA , iniMSK, Inf , false );
        
        DATA.GM = cat(2,DATA.GM,GMnew);
        DATA.BB = cat(2,DATA.BB,BBnew);
        DATA.idx = cat(1, DATA.idx, DATA.idxnew);
        
        % Progress Percentage
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTC3D_showProgressCM([],ProgressOldString,ProgressValue, jj == length(DATA.ConnectedSegmentsOS.List));
        
    end
    
    DATA.OrigSeedsIDX = DATA.idx;
    
else
    
    [GMnew,BBnew] = VTC3D_computeGeodesicMaps( DATA, [], Inf, true );
    
    DATA.GM = cat(2,DATA.GM,GMnew);
    DATA.BB = cat(2,DATA.BB,BBnew);
    DATA.idx = cat(1, DATA.idx, DATA.idxnew);
    
end

function DATA = determineGeodesicGraph(DATA)

if ~isempty(DATA.pairs) % The new distance metric MUST be computed for those pairs of Nodes!
    
    DATA.exploredpairs = cat(1,DATA.exploredpairs,DATA.pairs);
    
    [DATA.GGM,DATA.minPath] = VTC3D_computeGeodesicGraphMatrix(DATA);
    
    DATA.exploredminPath = cat(2,DATA.exploredminPath,DATA.minPath);
    
    % Extract Minimum Spanning Tree
    [DATA.MST,~] = graphminspantree(sparse(DATA.GGM),'Method','Kruskal');
    
    % Refine and Update Graph's Nodes (actually only on the Minimum
    % Spanning Tree!)
    [DATA.idxnew,DATA.NewNodesFLAG,DATA.territorypairs] = VTC3D_updateMSTNodes(DATA);
    
else
    
    % Extract Minimum Spanning Tree
    [DATA.MST,~] = graphminspantree(sparse(DATA.GGM),'Method','Kruskal');
    
    DATA.NewNodesFLAG = false;
    
    % Resetting Internal Variables
    assert(isempty(DATA.idxnew),'Something Went Wrong: Algorithm is stopping with new nodes to evaluate!');
    DATA.pairs = [];
    DATA.minPath = [];
    DATA.territorypairs = [];
    
end

function DATA = VTC3D_exportGraphConnectivity(DATA)

PathDirName = DATA.PathFolder;

warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir(PathDirName);

DATAFileName = strcat([PathDirName,DATA.IMG.ImgName,'_ConnGeoGraph.mat']);

disp(['Exporting Connecting Geodesic Graph to: ',DATAFileName]);

if DATA.FullDump
    save(DATAFileName,'DATA','-v7.3');
else
    % Export Only Useful Infos
    DATAexport.IMG = DATA.IMG;
    DATAexport.GRID = DATA.GRID;
    DATAexport.SPDPowerFactor = DATA.SPDPowerFactor;
    DATAexport.VisibleNeighborhoodRadius_mm = DATA.VisibleNeighborhoodRadius_mm;
    DATAexport.NodeDistanceTHR_mm = DATA.NodeDistanceTHR_mm;
    DATAexport.GradientDescentStopTHR = DATA.GradientDescentStopTHR;
    DATAexport.PathFolder = DATA.PathFolder;
    DATAexport.QExplore = DATA.QExplore;
    DATAexport.VisibleNeighborhoodRadius_vx = DATA.VisibleNeighborhoodRadius_vx;
    DATAexport.VisibilityBall = DATA.VisibilityBall;
    DATAexport.GGM = DATA.GGM;
    DATAexport.exploredminPath = DATA.exploredminPath;
    DATAexport.MST = DATA.MST;
    DATAexport.OrigSeedsIDX = DATA.OrigSeedsIDX;
    DATAexport.ExtraSeedsIDX = DATA.ExtraSeedsIDX;
    
    % Overwrite
    DATA = DATAexport;
    
    % Exporting
    save(DATAFileName,'DATA','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gdp_Uniform = VTC3D_smoothResamplePath(gdp,padsize,UniformSpacing)

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

function [GM,BB] = VTC3D_computeGeodesicMaps( DATA , MSK , Gmax , showProgressFlag )

% Initialisation
BB = struct([]);
GM = struct([]);

SPD = DATA.SPD;
T = DATA.TF;

frame = VTC3D_getImageFrame3D(size(SPD));

idx = DATA.idxnew;
BB_old = DATA.BB;
territorypairs = DATA.territorypairs;
% Scaled Normalised Voxel Size
ScaledNormalisedVoxelSize = (DATA.GRID.physicalvoxelsize/max(DATA.GRID.physicalvoxelsize)) * DATA.GRID.Scale;

if isempty(MSK)
    MSK = DATA.MSK;
    idxEvalTHRflag = false;
else
    if ~isscalar(idx)
        idxEvalTHRflag = true;
    else
        idxEvalTHRflag = false;
    end
end

if showProgressFlag
    ProgressValue = 0;
    ProgressValue_step = 100.0/length(idx);
    ProgressMessage = 'Computing Geodesics:';
    ProgressOldString = VTC3D_showProgressCM(ProgressMessage,'',ProgressValue, false);
end

for jj = 1 : length(idx)
    
    if showProgressFlag
        % Progress Percentage
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTC3D_showProgressCM([],ProgressOldString,ProgressValue, jj == length(idx));
    end
    
    if ~isempty(BB_old) && ~isempty(territorypairs)
        
        B1 = VTC3D_retrieveOriginalMapFromCrop(BB_old,territorypairs(jj,1));
        B2 = VTC3D_retrieveOriginalMapFromCrop(BB_old,territorypairs(jj,2));
        
        Barriers = ~( ~B1 | ~B2 );
        Barriers = Barriers | ~MSK | frame ;
        
    else
        
        Barriers = false(size(SPD));
        Barriers = Barriers | ~MSK | frame ;
        
    end
    
    SS = false(size(SPD));
    SS(idx(jj)) = true;
    
    % Cropping
    BoundingBoxCrop = VTC3D_computeBoundingBox(~Barriers);
    offset = BoundingBoxCrop(1,:) - 1;
    SPDbbc = SPD( BoundingBoxCrop(1,1):BoundingBoxCrop(2,1) , BoundingBoxCrop(1,2):BoundingBoxCrop(2,2) , BoundingBoxCrop(1,3):BoundingBoxCrop(2,3) );
    Tbbc   = T  ( BoundingBoxCrop(1,1):BoundingBoxCrop(2,1) , BoundingBoxCrop(1,2):BoundingBoxCrop(2,2) , BoundingBoxCrop(1,3):BoundingBoxCrop(2,3) , : );
    SSbbc  = SS ( BoundingBoxCrop(1,1):BoundingBoxCrop(2,1) , BoundingBoxCrop(1,2):BoundingBoxCrop(2,2) , BoundingBoxCrop(1,3):BoundingBoxCrop(2,3) );
    Barriersbbc = Barriers( BoundingBoxCrop(1,1):BoundingBoxCrop(2,1) , BoundingBoxCrop(1,2):BoundingBoxCrop(2,2) , BoundingBoxCrop(1,3):BoundingBoxCrop(2,3) );
    
    if idxEvalTHRflag
        OthrSS = false(size(SPD));
        idxEvalTHR = idx( 1:length(idx) ~= jj);
        OthrSS(idxEvalTHR) = true;
        OthrSSbbc = OthrSS ( BoundingBoxCrop(1,1):BoundingBoxCrop(2,1) , BoundingBoxCrop(1,2):BoundingBoxCrop(2,2) , BoundingBoxCrop(1,3):BoundingBoxCrop(2,3) );
    end
    
    % Computing Geodesics
    GeodesicMap = mxAnisoDistanceTransform4Graph( SPDbbc , Tbbc, double(SSbbc), double(Barriersbbc), double(ScaledNormalisedVoxelSize), Gmax);
    
    GeodesicMap(GeodesicMap < 0) = Inf;

    if idxEvalTHRflag
        %if DATA.QExplore < 1
            realVals = GeodesicMap(~(isinf(GeodesicMap)|isnan(GeodesicMap)));
            bb = GeodesicMap > ( GeodesicMap(OthrSSbbc) + 2*std( realVals(:) ) );
        %else
        %    bb = GeodesicMap > quantile(GeodesicMap(~(isinf(GeodesicMap)|isnan(GeodesicMap))), DATA.QExplore );
        %end
    else
        bb = GeodesicMap > quantile(GeodesicMap(~(isinf(GeodesicMap)|isnan(GeodesicMap))), DATA.QExplore );
    end
    
    GeodesicMap(bb) = inf;
    
    boundingbox = VTC3D_computeBoundingBox(~bb);
    
    GM(jj).map = GeodesicMap( boundingbox(1,1):boundingbox(2,1) , boundingbox(1,2):boundingbox(2,2) , boundingbox(1,3):boundingbox(2,3) );
    GM(jj).boundingbox = boundingbox + [offset ; offset] ;

    GM(jj).size = size(SPD);
    
    BB(jj).map = bb( boundingbox(1,1):boundingbox(2,1) , boundingbox(1,2):boundingbox(2,2) , boundingbox(1,3):boundingbox(2,3) );
    BB(jj).boundingbox = boundingbox + [offset ; offset];

    BB(jj).size = size(SPD);

end

function [RepeatedValues,a_check,b_check] = VTC3D_repeatedValueSymm(a,b)%

if isempty(a) && isempty(b)

    RepeatedValues = [];
    a_check = [];
    b_check = [];
    
elseif isempty(a) && ~isempty(b)
    
    RepeatedValues = [];
    a_check = [];
    b_check = false([size(b,1),1]);
    
elseif ~isempty(a) && isempty(b)
    
    RepeatedValues = [];
    a_check = false([size(a,1),1]);
    b_check = [];
    
else % ~isempty(a) && ~isempty(b)
    
    % This is to optimize speed!
    if size(a,1)*size(b,1) < 10e7
        
        aBOX = repmat(reshape(a',[1,size(a,2),size(a,1)]),[size(b,1),1,1]);
        bBOX = repmat(b,[1,1,size(a,1)]);
        cBOX = repmat(fliplr(b),[1,1,size(a,1)]);
        
        D = sum(abs(aBOX-bBOX),2);
        E = sum(abs(aBOX-cBOX),2);
        
        idx = find(D == 0 | E == 0);
        
        [bpos,~,apos] = ind2sub(size(D),idx);
        
        RepeatedValues = a(apos,:);
        
        a_check = false(size(a,1),1);
        a_check(apos) = true;
        
        b_check = false(size(b,1),1);
        b_check(bpos) = true;
        
    else
        
        RepeatedValues = [];
        a_check = false(size(a,1),1);
        b_check = false(size(b,1),1);
        
        swapCHK = false;
        
        if size(b,1) < size(a,1)
            swapCHK = true;
            temp = a;
            a = b;
            b = temp;
            
            a_check = false(size(a,1),1);
            b_check = false(size(b,1),1);
            
            clear temp;
        end
        
        for jj = 1 : size(a,1)
            
            arep = repmat(a(jj,:),[size(b,1),1]);
            
            D = sum(abs(arep - b),2);
            E = sum(abs(arep - fliplr(b)),2);
            
            idx = find(D == 0 | E == 0);
            
            if ~isempty(idx)
                RepeatedValues = cat(1,RepeatedValues,a(jj,:));
                a_check(jj) = true;
                b_check(idx) = true;
            end
            
        end
        
        if swapCHK
            temp = a_check;
            a_check = b_check;
            b_check = temp;
            clear temp;
        end
        
    end
    
end

function [GGM,minPath] = VTC3D_computeGeodesicGraphMatrix(DATA)

GM = DATA.GM;
BB = DATA.BB;
idx = DATA.idx;
pairs = DATA.pairs;
GGMold = DATA.GGM;
SPD = DATA.SPD;

minPath = struct([]);

SPDperm = permute( SPD , [3 1 2] );

xG = 1:GM(1).size(2);
yG = 1:GM(1).size(1);
zG = 1:GM(1).size(3);

padsize = 9;
FineUniformSpacing = 0.1;
CoarseUniformSpacing = 1.0;

GGM = double(~diag(true(size(idx))));
GGM(logical(GGM)) = Inf;
GGM(1:size(GGMold,1),1:size(GGMold,2)) = GGMold;

if ~isempty(GGMold)
    GGMold(GGMold == Inf) = 0;
    
    [ST,~] = graphminspantree(sparse(GGMold),'Method','Kruskal');
    [a1,a2,gvals] = find(ST);
    
    idxold_maxGeodesicSTVals = inf * ones(1,size(GGMold,1));
    
    for kk = 1 : size(GGMold,1)
        idxold_nodemsk = a1 == kk | a2 == kk;
        if ~isempty(find(idxold_nodemsk, 1))
            idxold_maxGeodesicSTVals(kk) = max(gvals(idxold_nodemsk));
        else
            idxold_maxGeodesicSTVals(kk) = inf;
        end
    end
    
    IGVref = inf * ones(size(idx));
    IGVref( 1 : length(idxold_maxGeodesicSTVals)) = idxold_maxGeodesicSTVals(:);
    
else
    IGVref = inf * ones(size(idx));
end

entry = 0;

% Progress
ProgressValue = 0;
ProgressValue_step = 100.0/size(pairs,1);
ProgressMessage = 'Computing Geodesic Graph Adjacency Matrix:';
ProgressOldString = VTC3D_showProgressCM(ProgressMessage,'',ProgressValue, false);

for el = 1 : size(pairs,1)

    % Progress Percentage
    ProgressValue = ProgressValue + ProgressValue_step;
    ProgressOldString = VTC3D_showProgressCM([],ProgressOldString,ProgressValue, el == size(pairs,1) );
        
    BB1 = VTC3D_retrieveOriginalMapFromCrop(BB,pairs(el,1));
    BB2 = VTC3D_retrieveOriginalMapFromCrop(BB,pairs(el,2));
    
    boundingbox = VTC3D_computeBoundingBox(~BB1|~BB2);

    clear BB1 BB2;
    
    GM1 = VTC3D_retrieveOriginalMapFromCrop(GM,pairs(el,1));
    GM2 = VTC3D_retrieveOriginalMapFromCrop(GM,pairs(el,2));
    Enrgy = log( (GM1 + GM2) + abs( GM1 - GM2 ) );
    EEnrgyPerm = permute( (GM1 + GM2) + abs( GM1 - GM2 ) , [3 1 2] );
    
    [ptA1,ptA2,ptA3] = ind2sub(size(Enrgy),idx(pairs(el,1))); % pointA
    ptA = [ptA2,ptA1,ptA3];
    
    [ptB1,ptB2,ptB3] = ind2sub(size(Enrgy),idx(pairs(el,2))); % pointB
    ptB = [ptB2,ptB1,ptB3];
    
    idxC = find( Enrgy == min(Enrgy(:))); % mutual GLOBAL Minimum: pointC
    [ptC1,ptC2,ptC3] = ind2sub(size(Enrgy),idxC);
    
    ptC = [ptC2(1),ptC1(1),ptC3(1)];
    
    % Crop to Relative BoundingBox
    EnrgyCrp = Enrgy( boundingbox(1,1):boundingbox(2,1) , boundingbox(1,2):boundingbox(2,2) , boundingbox(1,3):boundingbox(2,3) );
    offset = (boundingbox(1,:) - 1);
    
    ptAcrp = ptA - [offset(2),offset(1),offset(3)];
    ptBcrp = ptB - [offset(2),offset(1),offset(3)];
    ptCcrp = ptC - [offset(2),offset(1),offset(3)];
    
    [minPathHalves,ValidFlag] = VTC3D_GDPathSearch3D(EnrgyCrp,[ptAcrp;ptBcrp],ptCcrp,IGVref(pairs(el,1)),offset,DATA.GradientDescentStopTHR);
    
    if ValidFlag % PATH EXTRACTION CONVERGED
        
        entry = entry + 1;
        
        gdp = [minPathHalves(1).sequence ; ptC ; flipud(minPathHalves(2).sequence)];
        
        gdp_Uniform = VTC3D_smoothResamplePath(gdp,padsize,FineUniformSpacing);
        
        minPath(entry).sequence = gdp_Uniform;
        minPath(entry).CGPathContinuous = VTC3D_smoothResamplePath(gdp,padsize,CoarseUniformSpacing);
        
        minPath(entry).geodesic = [ 0 , abs( diff( linterp_matlab( xG, yG, zG , EEnrgyPerm , minPath(entry).sequence(:,1) , minPath(entry).sequence(:,2) , minPath(entry).sequence(:,3) ) ) ) ]';
        minPath(entry).geodesicRAW = [ 0 , abs( diff( linterp_matlab( xG, yG, zG , EEnrgyPerm , minPath(entry).CGPathContinuous(:,1) , minPath(entry).CGPathContinuous(:,2) , minPath(entry).CGPathContinuous(:,3) ) ) ) ]';
        intensityRAW = linterp_matlab( xG, yG, zG , SPDperm , minPath(entry).CGPathContinuous(:,1) , minPath(entry).CGPathContinuous(:,2) , minPath(entry).CGPathContinuous(:,3) )';
        minPath(entry).intensity = linterp_matlab( xG, yG, zG , SPDperm , minPath(entry).sequence(:,1) , minPath(entry).sequence(:,2) , minPath(entry).sequence(:,3) )';
        
        minPath(entry).geodesic(isnan(minPath(entry).geodesic)) = 0;
        minPath(entry).geodesicRAW(isnan(minPath(entry).geodesicRAW)) = 0;
        minPath(entry).geodesicRAW = minPath(entry).geodesicRAW + (1./(intensityRAW).^2);
        minPath(entry).pairs = pairs(el,:);
        
        GGM(pairs(el,1),pairs(el,2)) = sum( minPath(entry).geodesic ) + max(1./(minPath(entry).intensity.^2)) + (sum(sqrt(sum(diff(minPath(entry).sequence,1).^2,2)),1)).^2; % Orig 
        GGM(pairs(el,2),pairs(el,1)) = sum( minPath(entry).geodesic ) + max(1./(minPath(entry).intensity.^2)) + (sum(sqrt(sum(diff(minPath(entry).sequence,1).^2,2)),1)).^2; % Orig 
        
    end 
    
end

function [idxnew,NewNodesFLAG,territorypairs] = VTC3D_updateMSTNodes(DATA)

%Initialisation
MST = DATA.MST;
minPath = DATA.exploredminPath;
Size = size(DATA.SPD);
idx = DATA.idx;
physicalVoxelSize_mm = DATA.GRID.physicalvoxelsize ./ DATA.GRID.Scale;
NodeClosenessThreshold_mm = DATA.NodeDistanceTHR_mm;

idxnew = [];
territorypairs = [];
coords = [];

[a1,a2,~] = find(MST);
validpairs = zeros(size(minPath,2),2);
for jj = 1 : size(minPath,2)
   validpairs(jj,:) = minPath(jj).pairs; 
end

[~,elms2considerCHK,~] = VTC3D_repeatedValueSymm(validpairs,[a1,a2]);

elms2consider = find(elms2considerCHK);

for jj = 1 : length(elms2consider)
    
    [~,pos] = min( abs( cumsum(sqrt(sum(diff(minPath(elms2consider(jj)).sequence,1,1).^2,2))) - (sum(sqrt(sum(diff(minPath(elms2consider(jj)).sequence,1,1).^2,2)))/2) ) );% find the MIDPOINT!
    
    % along-path distance bewteen the mid-point and the respective
    % connecting nodes (that are the endpoints of the minPath). 
    edm1 = sum( sqrt( sum( diff( minPath(elms2consider(jj)).sequence(  1  : pos, : ) .* repmat(physicalVoxelSize_mm,[ size(minPath(elms2consider(jj)).sequence(  1  : pos, : ),1 ) ,1] ) , 1 ).^2 , 2 ) ) );
    edm2 = sum( sqrt( sum( diff( minPath(elms2consider(jj)).sequence( pos : end, : ) .* repmat(physicalVoxelSize_mm,[ size(minPath(elms2consider(jj)).sequence( pos : end, : ),1 ) ,1] ) , 1 ).^2 , 2 ) ) );
    
    % if the mid-point is far enough from the endpoints of the minPath:
    % Generate a new node @ mid-point
    if edm1 >= NodeClosenessThreshold_mm && edm2 >= NodeClosenessThreshold_mm
        coords = cat(1,coords,round(minPath(elms2consider(jj)).sequence(pos,:)));
        territorypairs = cat(1,territorypairs,minPath(elms2consider(jj)).pairs);
    end
    
end

if ~isempty(coords)
    [idxnew,idxnewCHK1,~] = unique( sub2ind(Size,coords(:,2),coords(:,1),coords(:,3)) );
    territorypairs = territorypairs(idxnewCHK1,:);
    
    SS = false(Size);
    SS(idx) = true;
    SS = imdilate(SS,true([3 3 3]));
    idxDilated = find(SS);
    
    [~,idxnewCHK2] = VTC3D_notUnique(idxnew,idxDilated);
    idxnew = idxnew(~idxnewCHK2);
    territorypairs = territorypairs(~idxnewCHK2,:);
    
end

if isempty(idxnew)
    NewNodesFLAG = false;
else
    NewNodesFLAG = true;
end

function [RepeatedValues,a_check] = VTC3D_notUnique(a,b)%
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

function boundingbox = VTC3D_computeBoundingBox(LogicalMap)

idx = find(LogicalMap);
[rr,cc,pp] = ind2sub(size(LogicalMap),idx);
boundingbox = [min(rr) - 1 , min(cc) - 1 , min(pp) - 1;...
               max(rr) + 1 , max(cc) + 1 , max(pp) + 1];
           
function Map = VTC3D_retrieveOriginalMapFromCrop(Strct,idx)

if islogical(Strct(idx).map)
    Map = true(Strct(idx).size);
    Map(Strct(idx).boundingbox(1,1):Strct(idx).boundingbox(2,1),...
        Strct(idx).boundingbox(1,2):Strct(idx).boundingbox(2,2),...
        Strct(idx).boundingbox(1,3):Strct(idx).boundingbox(2,3)) = Strct(idx).map;
else
    Map = inf * ones(Strct(idx).size);
    Map(Strct(idx).boundingbox(1,1):Strct(idx).boundingbox(2,1),...
        Strct(idx).boundingbox(1,2):Strct(idx).boundingbox(2,2),...
        Strct(idx).boundingbox(1,3):Strct(idx).boundingbox(2,3)) = Strct(idx).map;
end

function [GDP,tgt_pnt] = VTC3D_GDPathSearch3D(U,str_pt,end_pnt,IGVref,offset,endpointThrshld)

Unewmax = max(U(~isnan(U) & ~isinf(U))) + eps;
U(abs(U) == inf | isnan(U)) = Unewmax;
pU = permute(U,[3 1 2]);

[Ux,Uy,Uz] = gradient3Dcpp_mex(U); % EXTERNAL MEX FILE!
pUx = permute(Ux,[3 1 2]);
pUy = permute(Uy,[3 1 2]);
pUz = permute(Uz,[3 1 2]);

xG = 1:size(U,1); yG = 1:size(U,2); zG = 1:size(U,3);

GAMMA_ORIG = 0.1 * mean([max(abs(Ux(:))),max(abs(Uy(:))),max(abs(Uz(:)))]);

adjacentThrshld = 1/3;
maxiter = 10 * round( (3*( sqrt(sum(size(U).^2))/2) ) /adjacentThrshld );

% Because of the Matlab XY flip
str_pt = str_pt(:,[2,1,3]);
end_pnt = end_pnt([2,1,3]);
offset = offset([2,1,3]);

tgt_pnt = true(size(str_pt,1),1);

GDP = struct([]);

for el = 1 : size(str_pt,1)

    gdp = str_pt(el,:);
    
    % Remember: do not start the while cycle if the str_pnt(any) == end_pnt !!!
    if sum( (gdp(end,:) - end_pnt).^2, 2 ) >= endpointThrshld
        StopSearch = false;
    else
        StopSearch = true;
    end
    
    itr = 0;
    tempIGV = 0;
    GAMMA = GAMMA_ORIG;
    
    while ~StopSearch
        
        itr = itr + 1;
                          
        nxtstp = gdp(end,:) - GAMMA.*[ linterp_matlab( yG , xG , zG , pUy , gdp(end,2) , gdp(end,1) , gdp(end,3) ) , ...
                                       linterp_matlab( yG , xG , zG , pUx , gdp(end,2) , gdp(end,1) , gdp(end,3) ) , ...
                                       linterp_matlab( yG , xG , zG , pUz , gdp(end,2) , gdp(end,1) , gdp(end,3) )];

        if sum( (gdp(end,:) - nxtstp).^2, 2 ) < adjacentThrshld
            
            gdp = cat(1,gdp,nxtstp);
            
            boost = (sum((gdp(end,:) - gdp(end-1,:)).^2))^(1/2);
            GAMMA = GAMMA*(1+(1/boost));
            
            tempIGV = tempIGV + abs( linterp_matlab( yG , xG , zG, pU , gdp(  end,2) , gdp(  end,1) , gdp(  end,3) ) - ...
                                     linterp_matlab( yG , xG , zG, pU , gdp(end-1,2) , gdp(end-1,1) , gdp(end-1,3) ) );                
            
            if (sum( (gdp(end,:) - end_pnt).^2 , 2 ))^(1/2) < endpointThrshld
                StopSearch = true;
            end
            
        else
            GAMMA = GAMMA/2;
            if GAMMA <= 0
                GAMMA = GAMMA_ORIG;
            end
        end
        
        % Stop over not convergence
        if itr > maxiter || tempIGV > IGVref
            StopSearch = true;
            tgt_pnt(el) = false;
        end
        
    end
    
    GDP(el).sequence = [gdp(:,2) + offset(1) , gdp(:,1) + offset(2), gdp(:,3) + offset(3)];

end

tgt_pnt = logical(prod(tgt_pnt));

function frame = VTC3D_getImageFrame3D(Size)
frame = false(Size);

frame(1,:,:) = true;
frame(end,:,:) = true;
frame(:,1,:) = true;
frame(:,end,:) = true;
frame(:,:,1) = true;
frame(:,:,end) = true;

function [EIG,SSAFW] = VTC3D_get3DImageVoxelTensorLE(Input,varargin)%
% Use this function to compute with the ImageVoxelTensor.
% It is based on a combined Scale-Space approach and EigenDecomposition.
% The Scale-Space approach is similar to that used for Vesselness response.
% The EigenDecomposition makes use of divide-and-conquer techniques to
% estimate the EigenValues and EigenVectors of the Hessian Matrix relative
% to each Voxel of the Image I.
% Averaging and Interpolation of Tensors over different scales follows
% Log-Euclidean Metrics for mathematical computations. 
% Please see: 'Log-Euclidean Metrics for Fast and Simple Calculus on
% Diffusion Tensors' - Vincent Arsigny et Al. 2006
%
% Inputs:
%         Input:   3D Image of size = [h,w,d] (real-valued double voxels)
%                                           OR
%                  Structure containing the Scale-Space Analysis FrameWork.
%        [mask]:   logical mask of size = [h,w,d] with the ROI.
% [noise_sigma]:   amount of additional noise to the image (scalar value)  
%
% Outputs: 
%           EIG:   Structure containing the EigenValues (El1 El2 El3) maps
%                  each of size [h,w,d], and the EigenVectors (Ev1 Ev2 Ev3)
%                  maps each of size [h,w,d,3]. Using these values it is
%                  possible to determine the voxel Tensor.
%         SSAFW:   Structure containing the Scale-Space Analysis FrameWork.
%
% Notice that for a large image I, the algorithm may take up to several
% minutes to compute the whole Tensor.

% Get the Input and identify it.
if isstruct(Input)
    SSAFW = Input;
    I_temp = SSAFW.ImgOrig;
    retrieveSSAFWflag = false;
else
    assert(length(size(Input)) == 3,'Input MUST be a 3D Image Volume OR the Structure containing the Scale-Space Analysis FrameWork!');
    I = Input;
    I_temp = I;
    retrieveSSAFWflag = true;
end

% Get the Mask (if any)
if nargin > 1
    Msk = logical(varargin{1});
    if ~isequal(size(Msk),size(I_temp))
        disp('The Input Mask is NOT Valid! Default Mask will be applied: Whole Image;');
        Msk = true(size(I_temp));
    end
else
    Msk = true(size(I_temp));
end
clear I_temp;

% Apply additional noise (if any)
if length(varargin) > 1 && retrieveSSAFWflag
    
    noise_sigma = varargin{2};
    if ~isscalar(noise_sigma)
        disp('Noise Level must be Scalar! Default value is applied: noise_sigma = 10e-4;');
        noise_sigma = 10e-4;
    end
    % Adding Gaussian Noise (useful for synthetic images or phantoms)
    I = I + noise_sigma*rand(size(I));
    % Normalization of I within [0,1] range
    I = (I-min(I(:)))./(max(I(:))-min(I(:)));
    
elseif length(varargin) > 1 && ~retrieveSSAFWflag
    disp('Noise Level will be Ignored! A Scale-Space Analysis FrameWork has alreaby been provided!');
end

if retrieveSSAFWflag
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if prod(round(size(I)*0.1) >= 3)
        ScaleRange = [0.1:0.1:1.0];
    else
        ScaleRange = linspace(max(3./size(I)),1.0,10);
    end
    VoxelSize = [1 1 1];% This must agree with the Image infos!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % The following function searches first an existing backup file
    % containing the Scale-Space Analysis FrameWork as temporary file in
    % the current directory. If it does not exist, then the function
    % creates it, and stores it as backup file named 'SSAFW_temp.mat'.
    % At the end of the whole processing, the backup file is supposed to be
    % deleted. In this way, the workspace variable 'SSAFW' can be cleared
    % several times during the processing (i.e. the Anisotropic Fast
    % Marching), and loaded again, when necessary, from the backup file.
    SSAFW = VTC3D_retrieveSSAFW3D(I,ScaleRange,VoxelSize);
end

EIG = VTC3D_computeSS3DEigenDecomposition(SSAFW,Msk);

function SSAFW = VTC3D_retrieveSSAFW3D(I,ScaleRange,VoxelSize)

if exist('SSAFW3D_temp.mat','file') == 2
    load SSAFW3D_temp.mat;
    disp('Scale-Space Analysis FrameWork: Loaded');
else
    disp('Scale-Space Analysis FrameWork: Creating...');
    SSAFW = VTC3D_get3DScaleSpaceAnalysisFrameWork(I,ScaleRange,VoxelSize);
end

function SSAFW = VTC3D_get3DScaleSpaceAnalysisFrameWork(I,ScaleRange,VoxelSize)

assert(length(size(I)) == 3 , 'The Input I must be a 3D Scalar Image Volume!');
assert(isvector(ScaleRange) , 'The Input ScaleRange MUST be a 1xn Vector ]0,1]!');
assert(length(VoxelSize) == 3 , 'The Input VoxelSize MUST be a 1x3 Vector!');

SSAFW = struct('ImgOrig',I,'ScaleRange',ScaleRange,'VoxelSize',VoxelSize,...
               'GridOrig',[],'GridDwnSmpl',[],'ImgDwnSmpl',[]);

[GridOrig.X,GridOrig.Y,GridOrig.Z] = ndgrid(1:size(I,1),1:size(I,2),1:size(I,3));

ReductionFactor = 1./ScaleRange;

% Relative VoxelSize
VoxelSize = VoxelSize./(min(VoxelSize));

% Full Width at Half Maximum for Gaussian Distribution of sigma = 1
FWHM = 2*sqrt(2*log(2));

for s = 1 : length(ScaleRange)
    Sigmas = diag( ( diag( (ReductionFactor(s) .* VoxelSize ) ./ FWHM).^2 - diag(VoxelSize ./ FWHM).^2 ).^(1/2) );
    
    Xv = linspace(1,size(I,1),round(size(I,1)*ScaleRange(s)));
    Yv = linspace(1,size(I,2),round(size(I,2)*ScaleRange(s)));
    Zv = linspace(1,size(I,3),round(size(I,3)*ScaleRange(s)));
    [GridDwnSmpl(s).X,GridDwnSmpl(s).Y,GridDwnSmpl(s).Z] = ndgrid(Xv,Yv,Zv);
    clear Xv Yv Zv;
    
    % Gaussian Blur to Avoid Aliasing
    IBlur = VTC3D_gaussBlur3D(I,Sigmas); % Non-constant speed across scales, but accurate.
    clear Sigmas;
    
    % DownSampling the Image by Interpolation
    ImgDwnSmpl(s).Img = interpn( GridOrig.X,GridOrig.Y,GridOrig.Z, IBlur, GridDwnSmpl(s).X,GridDwnSmpl(s).Y,GridDwnSmpl(s).Z , 'linear');
    SSAFW.ImgGaussBlurDwnUpSmpl(s).Img = interpn( GridDwnSmpl(s).X,GridDwnSmpl(s).Y,GridDwnSmpl(s).Z, ImgDwnSmpl(s).Img, GridOrig.X,GridOrig.Y,GridOrig.Z, 'linear');
    
    clear IBlur;
end

SSAFW.GridOrig = GridOrig;
SSAFW.GridDwnSmpl = GridDwnSmpl;
SSAFW.ImgDwnSmpl = ImgDwnSmpl;

save('SSAFW3D_temp.mat','SSAFW');

function VTC3D_deleteSSAFW3D()
if exist('SSAFW3D_temp.mat','file') == 2
    delete('SSAFW3D_temp.mat');
end

function EIG = VTC3D_computeSS3DEigenDecomposition(SSAFW,Msk)
% Initializing Mid-Term Outputs:
T11LE = zeros(size(SSAFW.ImgOrig));
T12LE = zeros(size(SSAFW.ImgOrig));
T13LE = zeros(size(SSAFW.ImgOrig));
T22LE = zeros(size(SSAFW.ImgOrig));
T23LE = zeros(size(SSAFW.ImgOrig));
T33LE = zeros(size(SSAFW.ImgOrig));

Accumulator = zeros(size(SSAFW.ImgOrig));

AnisotropyExpFactor = 1.0; % Enhancement Factor for Anisotropy. 
    
SE = logical(cat(3,[0 0 0;0 1 0;0 0 0],[0 1 0;1 1 1;0 1 0],[0 0 0;0 1 0;0 0 0]));

for s = 1 : length(SSAFW.ScaleRange)
    
    % Dilate 3 times the mask with a 3D cross SE
    MskDwnSmpl = imdilate(...
                  imdilate(...
                   imdilate(...
                            interpn(SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, ...
                                    Msk,...
                                    SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                                    'nearest'), ...
                            SE),...
                           SE),...
                          SE);
    
    % EigenDecomposition
    [El1DwnSmpl,El2DwnSmpl,El3DwnSmpl,Ev1DwnSmpl,Ev2DwnSmpl,Ev3DwnSmpl,~] = VTC3D_get3DImageEIGs(SSAFW.ImgDwnSmpl(s).Img , MskDwnSmpl , true ); % These MUST BE Sorted By Magnitude!
    
    if sum(sum(sum( abs( dot(Ev1DwnSmpl,Ev2DwnSmpl,4) ) > 1e-3))) || sum(sum(sum( abs( dot(Ev1DwnSmpl,Ev3DwnSmpl,4) ) > 1e-3))) || sum(sum(sum( abs( dot(Ev2DwnSmpl,Ev3DwnSmpl,4) ) > 1e-3)))
        disp('Warning: The Eigenvector Bases Computed are not Orthogonal! -- Please Correct this Issue!');
    end
    
    % Defining the Anisotropy Metric related to the Hessian
    % EigenDecomposition (based on Frangi 1998) - Blobness Measures of the
    % 2D/3D Ellipsoid, i.e. the BlobMeasure and the BlobCrossRef and the BlobReference in 3D.
    
    % Blobness Measure (BlobMeasure):
    % As defined, it indicates how squeezes the tensor is in 3D along the vessel direction.
    % The BlobMeasure value is oriented along the direction of the vessel.
    BlobMeasureDwnSmpl = ( abs(El1DwnSmpl)./sqrt( abs( El2DwnSmpl .* El3DwnSmpl ) ) );
    
    % Blobness Cross-Sectional Reference (BlobCrossRef):
    % As defined, it indicates the largest cross-sectional component of the ellipsoid.
    % The BlobCrossRef value is oriented along the First Othogonal direction of the vessel.
    BlobCrossRefDwnSmpl = ( abs(El2DwnSmpl) ./ abs( El3DwnSmpl ) );
    
    % The last Blobness Measures Component is the BlobReference, which will be
    % equal to 1.
    % N.B.: the anisotropic Ratio must be preserved!
    % The BlobReference is oriented along the Second Orthogonal direction of the vessel. 
    BlobReferenceDwnSmpl = ones(size(BlobMeasureDwnSmpl));
    
    % The Blobness Measures themselves represent the anisotropic tensor of
    % the pixel/voxel.
    % Specifically, the Tensor is defined in 2D as:
    % Tensor = (BlobMeasure * (Evec_alongVessel*Evec_alongVessel')) + (BlobCrossRef * (Evec_orthog1Vessel*Evec_orthog1Vessel')) + (BlobReference * (Evec_orthog2Vessel*Evec_orthog2Vessel'))  
    % Given Evec_alongVessel and Evec_orthogVessel actually orthogonal.
    % N.B.: Ev1DwnSmpl = Evec_alongVessel
    %       Ev2DwnSmpl = Evec_orthog1Vessel
    %       Ev3DwnSmpl = Evec_orthog2Vessel
    
    % In order to upsample the Tensor Vector Field (TF), and to average it
    % across scales, Tensors must be converted into the Log-Euclidean
    % Space. (MATRIX LOG, MATRIX EXP)
    
    % Transforming the Euclidean Tensor Components into Log-Euclidean
    % Space. Output: 6 components of the TensorLE (it is symmetric and 3D!)
    %
    %            | T11LE   T12LE   T13LE |
    % TensorLE = | T12LE   T22LE   T23LE |
    %            | T13LE   T23LE   T33LE |
    

    %%%%%%%%%%%%%%%%%%%%%%%% LOG EUCLIDEAN SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %[DwnSmplT11LE,DwnSmplT12LE,DwnSmplT13LE,DwnSmplT22LE,DwnSmplT23LE,DwnSmplT33LE] = convert3DTensorTo6LogEuclidComponents(El1DwnSmpl,El2DwnSmpl,El3DwnSmpl,Ev1DwnSmpl,Ev2DwnSmpl,Ev3DwnSmpl);
    [DwnSmplT11LE,DwnSmplT12LE,DwnSmplT13LE,DwnSmplT22LE,DwnSmplT23LE,DwnSmplT33LE] = VTC3D_convert3DTensorTo6LogEuclidComponents(BlobMeasureDwnSmpl,BlobCrossRefDwnSmpl,BlobReferenceDwnSmpl,Ev1DwnSmpl,Ev2DwnSmpl,Ev3DwnSmpl);
    
    clear El1DwnSmpl El2DwnSmpl El3DwnSmpl Ev1DwnSmpl Ev2DwnSmpl Ev3DwnSmpl;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Up-Sampling the 6 Components %
    
    UpSmplT11LE = interpn(SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                          DwnSmplT11LE,...
                          SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, 'linear', NaN);
                      
    UpSmplT12LE = interpn(SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                          DwnSmplT12LE,...
                          SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, 'linear', NaN);
    
    UpSmplT13LE = interpn(SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                          DwnSmplT13LE,...
                          SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, 'linear', NaN);
                      
    UpSmplT22LE = interpn(SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                          DwnSmplT22LE,...
                          SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, 'linear', NaN);
                      
    UpSmplT23LE = interpn(SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                          DwnSmplT23LE,...
                          SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, 'linear', NaN);
                      
    UpSmplT33LE = interpn(SSAFW.GridDwnSmpl(s).X,SSAFW.GridDwnSmpl(s).Y,SSAFW.GridDwnSmpl(s).Z,...
                          DwnSmplT33LE,...
                          SSAFW.GridOrig.X,SSAFW.GridOrig.Y,SSAFW.GridOrig.Z, 'linear', NaN);
                      
    clear DwnSmplT*;
    
    % Other possible Operations on the 6Components Should be DONE HERE
    % while they are still in the Log-Euclidean Space!
    % ...
    % ... Other Operations Here [Optional] (e.g. Tensor Regularization)
    
    
    
    % Checking:
    % 1) the interpolation only generated positive definite (PD)
    % Tensors (the eigenvalues must be all positive!)
    % 2) the interpolated bases are orthogonal
    [El1_chk,El2_chk,El3_chk,Ev1_chk,Ev2_chk,Ev3_chk] = VTC3D_convert6LogEuclidComponentsTo3DTensor(UpSmplT11LE,UpSmplT12LE,UpSmplT13LE,UpSmplT22LE,UpSmplT23LE,UpSmplT33LE,Msk);
    
    if (sum(sum(sum(El1_chk < 0))) > 0) || (sum(sum(sum(El2_chk < 0))) > 0) || (sum(sum(sum(El3_chk < 0))) > 0)
        disp('Eigenvalues are NEGATIVE after the linear Upsampling in LogEuclidean space!');
    end
    
    if sum(sum(sum( abs( dot(Ev1_chk,Ev2_chk,4) ) > 1e-3))) || sum(sum(sum( abs( dot(Ev1_chk,Ev3_chk,4) ) > 1e-3))) || sum(sum(sum( abs( dot(Ev2_chk,Ev3_chk,4) ) > 1e-3)))
        disp('Warning: The Up-Sampled Eigenvector Bases are not mutually Orthogonal!');
    end
    
    clear El1_chk El2_chk El3_chk Ev1_chk Ev2_chk Ev3_chk;
    
    % Accumulating the 6 Components over Scales
    [UpSmplT11LE,CHK] = VTC3D_switchNaN2zero(UpSmplT11LE);
    [UpSmplT12LE, ~ ] = VTC3D_switchNaN2zero(UpSmplT12LE);
    [UpSmplT13LE, ~ ] = VTC3D_switchNaN2zero(UpSmplT13LE);
    [UpSmplT22LE, ~ ] = VTC3D_switchNaN2zero(UpSmplT22LE);
    [UpSmplT23LE, ~ ] = VTC3D_switchNaN2zero(UpSmplT23LE);
    [UpSmplT33LE, ~ ] = VTC3D_switchNaN2zero(UpSmplT33LE);
    
    %%% Weighted Version: Weight = UpSmpl(DwnSmpl(GaussBlur(Img)))  -- Consistent with the Tensor Estimation %%%
    Accumulator = Accumulator + SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK;
    T11LE = T11LE + (UpSmplT11LE .* SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK );
    T12LE = T12LE + (UpSmplT12LE .* SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK );
    T13LE = T13LE + (UpSmplT13LE .* SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK );
    T22LE = T22LE + (UpSmplT22LE .* SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK );
    T23LE = T23LE + (UpSmplT23LE .* SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK );
    T33LE = T33LE + (UpSmplT33LE .* SSAFW.ImgGaussBlurDwnUpSmpl(s).Img .* CHK );
    
    clear UpSmplT*;
    
end

Accumulator(Accumulator == 0) = 1;

% Estimating the Final EigenValues and EigenVectors from multiple Scales
% Averaging process: (Log-Mean or Log-WeightedMean over the scales)
avgT11LE = VTC3D_switchNaN2zero(T11LE./Accumulator);
avgT12LE = VTC3D_switchNaN2zero(T12LE./Accumulator);
avgT13LE = VTC3D_switchNaN2zero(T13LE./Accumulator);
avgT22LE = VTC3D_switchNaN2zero(T22LE./Accumulator);
avgT23LE = VTC3D_switchNaN2zero(T23LE./Accumulator);
avgT33LE = VTC3D_switchNaN2zero(T33LE./Accumulator);

% Converting back to Euclidean Space: Retrieve only the 6 Components of the
% 3D Tensor useful for the ellipsoid representation
[El1,El2,El3,Ev1,Ev2,Ev3] = VTC3D_convert6LogEuclidComponentsTo3DTensor(avgT11LE,avgT12LE,avgT13LE,avgT22LE,avgT23LE,avgT33LE,Msk);

% FINAL check and possible Warning
if (sum(sum(sum(El1 < 0))) > 0) || (sum(sum(sum(El2 < 0))) > 0) || (sum(sum(sum(El3 < 0))) > 0)
    disp('Some FINAL Tensors are not SPD!');
end

if sum(sum(sum( abs( dot(Ev1,Ev2,4) ) > 1e-3))) || sum(sum(sum( abs( dot(Ev1,Ev3,4) ) > 1e-3))) || sum(sum(sum( abs( dot(Ev2,Ev3,4) ) > 1e-3)))
    disp('Warning: Some of the FINAL Eigenvector Bases are not mutually Orthogonal!');
end

% Correction for possible inconsistent values -> Isotropy
El1(isnan(El1)|isinf(El1)) = 1;
El2(isnan(El2)|isinf(El2)) = 1;
El3(isnan(El3)|isinf(El3)) = 1;

Ev1v1 = Ev1(:,:,:,1);
Ev1v2 = Ev1(:,:,:,2);
Ev1v3 = Ev1(:,:,:,3);
Ev1v1(isnan(Ev1v1)|isinf(Ev1v1)) = 1;
Ev1v2(isnan(Ev1v2)|isinf(Ev1v2)) = 0;
Ev1v3(isnan(Ev1v3)|isinf(Ev1v3)) = 0;
Ev1 = cat(4,Ev1v1,Ev1v2,Ev1v3);

Ev2v1 = Ev2(:,:,:,1);
Ev2v2 = Ev2(:,:,:,2);
Ev2v3 = Ev2(:,:,:,3);
Ev2v1(isnan(Ev2v1)|isinf(Ev2v1)) = 0;
Ev2v2(isnan(Ev2v2)|isinf(Ev2v2)) = 1;
Ev2v3(isnan(Ev2v3)|isinf(Ev2v3)) = 0;
Ev2 = cat(4,Ev2v1,Ev2v2,Ev2v3);

Ev3v1 = Ev3(:,:,:,1);
Ev3v2 = Ev3(:,:,:,2);
Ev3v3 = Ev3(:,:,:,3);
Ev3v1(isnan(Ev3v1)|isinf(Ev3v1)) = 0;
Ev3v2(isnan(Ev3v2)|isinf(Ev3v2)) = 0;
Ev3v3(isnan(Ev3v3)|isinf(Ev3v3)) = 1;
Ev3 = cat(4,Ev3v1,Ev3v2,Ev3v3);

% Resulting Tensors Assignment
EIG.El1 = El1.^AnisotropyExpFactor;
EIG.El2 = El2.^AnisotropyExpFactor;
EIG.El3 = El3.^AnisotropyExpFactor;
EIG.Ev1 = Ev1;
EIG.Ev2 = Ev2;
EIG.Ev3 = Ev3;

clear El* Ev*;

function IBlur = VTC3D_gaussBlur3D(I,sigma)%
% sigma can be scalar, or a vector containing the diagonal of the gaussian
% covariance matrix:
%
%               s1 0 0
% sigma = diag( 0 s2 0 ) = [s1 s2 s3]
%               0 0 s3

if isscalar(sigma)
    sigma = sigma * [1 1 1];
end

if (isvector(sigma) && sigma(1)~= 0 && sigma(2)~= 0 && sigma(3)~= 0)
    
    KernelSize = (round(3*sigma(:)) *2) +1;

    Ipad = padarray(I,[KernelSize(1),KernelSize(2),KernelSize(3)],'replicate');
    Ipadfreq = fftn(Ipad,size(Ipad));
    
    [XGk,YGk,ZGk] =   ndgrid(-floor(KernelSize(1)/2):1:floor(KernelSize(1)/2),...
                             -floor(KernelSize(2)/2):1:floor(KernelSize(2)/2),...
                             -floor(KernelSize(3)/2):1:floor(KernelSize(3)/2));
                         
    gain = 1/(sqrt((2*pi)^length(size(I)) * (prod(sigma.^2))));
    exponent = -( ((XGk).^2)/(2*sigma(1)^2) + ((YGk).^2)/(2*sigma(2)^2) + ((ZGk).^2)/(2*sigma(3)^2) );
    Gkernel = gain.*exp(exponent);
    
    Gkernel = Gkernel./sum(Gkernel(:));
    
    GKfreq = fftn(Gkernel,size(Ipad));
    
    Ifiltpad = real( ifftn( Ipadfreq .* GKfreq ) );

    SpatialDelay = (KernelSize - 1) / 2;
    
    IBlur = Ifiltpad( KernelSize(1) + SpatialDelay(1) + 1 : end - KernelSize(1) + SpatialDelay(1) , ...
                      KernelSize(2) + SpatialDelay(2) + 1 : end - KernelSize(2) + SpatialDelay(2) , ...
                      KernelSize(3) + SpatialDelay(3) + 1 : end - KernelSize(3) + SpatialDelay(3)  );
    
    
else
    IBlur = I; % No blurring applied if (one of) the sigma is equal to 0!
end

function [El1,El2,El3,Ev1,Ev2,Ev3,CHK] = VTC3D_get3DImageEIGs(I,Msk,SortByMagnitude)
% Use this function to compute with feasible complexity the
% EigenDecomposition of a 3D Volumetric Image I.
% Inputs:
%             I:   3D Image of size = [h,w,d] (real-valued double voxels)
%           Msk:   logical mask of size = [h,w,d] with the ROI.
%
% Outputs:
% El1, El2, El3:   Eigenvalue Maps of size [h,w,d] each, sorted s.t. abs(El1) <= abs(El2) <= abs(El3) 
% Ev1, Ev2, Ev3:   Eigenvector Maps of size [h,w,d,3] each, relative to the sorted Eigenvalue Maps. 
%
% Notice that for a large image I, the algorithm may take up to a couple of
% minutes to compute the whole EIG Maps.

if ~isstruct(I)
    Size = size(I);
    % Computing Gradient and Hessian
    [Ix,Iy,Iz] = gradient(I);
    [Ixx,Ixy,Ixz] = gradient(Ix);
    [ ~ ,Iyy,Iyz] = gradient(Iy);
    [ ~ , ~ ,Izz] = gradient(Iz);
    
    % Struct
    IMG.Ixx = Ixx; IMG.Ixy = Ixy; IMG.Ixz = Ixz;
    IMG.Iyy = Iyy; IMG.Iyz = Iyz; IMG.Izz = Izz;
else
    Size = size(I.Ixx);
    IMG = I;
end

idx = find(Msk);

% Output Initialization
El1 = NaN(Size);
Ev1u = NaN(Size);
Ev1v = NaN(Size);
Ev1w = NaN(Size);

El2 = NaN(Size);
Ev2u = NaN(Size);
Ev2v = NaN(Size);
Ev2w = NaN(Size);

El3 = NaN(Size);
Ev3u = NaN(Size);
Ev3v = NaN(Size);
Ev3w = NaN(Size);

CHK = NaN(Size);

[VEval1NS, VEval2NS, VEval3NS, VEvec1NS, VEvec2NS, VEvec3NS, idxReal] = VTC3D_Compute3DEIGs(IMG,idx,SortByMagnitude);

% -- EIG1 --
El1(idx) = VEval1NS;
Ev1u(idx) = VEvec1NS(:,1);
Ev1v(idx) = VEvec1NS(:,2);
Ev1w(idx) = VEvec1NS(:,3);
Ev1 = cat(4,Ev1u,Ev1v,Ev1w);

% -- EIG2 --
El2(idx) = VEval2NS;
Ev2u(idx) = VEvec2NS(:,1);
Ev2v(idx) = VEvec2NS(:,2);
Ev2w(idx) = VEvec2NS(:,3);
Ev2 = cat(4,Ev2u,Ev2v,Ev2w);

% -- EIG3 --
El3(idx) = VEval3NS;
Ev3u(idx) = VEvec3NS(:,1);
Ev3v(idx) = VEvec3NS(:,2);
Ev3w(idx) = VEvec3NS(:,3);
Ev3 = cat(4,Ev3u,Ev3v,Ev3w);

% -- EIGs Consistency --
CHK(idx) = idxReal;

function [VEval1NS,VEval2NS,VEval3NS,VEvec1NS,VEvec2NS,VEvec3NS,idxReal] = VTC3D_Compute3DEIGs(VG,idx,SortByMagnitude)

block_len = 24^3; %linearized length of a 24x24x24 cube of voxels
num_block = ceil(size(idx,1)/block_len);

% Initialization of the linearized EigenValues
VEval1NS = NaN(size(idx));
VEval2NS = NaN(size(idx));
VEval3NS = NaN(size(idx));

% Initialization of the linearized EigenVectors
VEvec1NS = NaN(size(idx,1),3);
VEvec2NS = NaN(size(idx,1),3);
VEvec3NS = NaN(size(idx,1),3);

% Initialization of the linearized ConsistencyBoolean
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
   
   [Evl1 , Evl2 , Evl3 , Evc1 , Evc2 , Evc3 , idx_discard] = VTC3D_processBlock3D(VG,idx_block,SortByMagnitude);
   
   VEval1NS(ini:fin,1) = Evl1;
   VEval2NS(ini:fin,1) = Evl2;
   VEval3NS(ini:fin,1) = Evl3;
   
   VEvec1NS(ini:fin,:) = Evc1;
   VEvec2NS(ini:fin,:) = Evc2;
   VEvec3NS(ini:fin,:) = Evc3;
   
   idxRealCheck(ini:fin) = ~idx_discard; % list of indices where Eigenvectors have to be considered ( i.e. PURE REAL, so Complex and NaN values are excluded)
   
   clear ini fin idx_block Evc*;
   
end

VEvec1NS(~idxRealCheck,:) = NaN;
VEvec2NS(~idxRealCheck,:) = NaN;
VEvec3NS(~idxRealCheck,:) = NaN;
idxReal = idxRealCheck;

function [Evl1,Evl2,Evl3,Evc1,Evc2,Evc3,idx_discard] = VTC3D_processBlock3D(VG,idx,SortByMagnitude)

% Computing EigenValues
VEval(:,1) = real(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3));
VEval(:,2) = real(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) - (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*1i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2);
VEval(:,3) = real(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) + (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*1i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2);

% Computing EigenVectors
% Eigenvector V of the first eigenvalue
VEvec1(:,1) = (VG.Ixz(idx).^2.*VG.Iyz(idx) - VG.Ixy(idx).*VG.Izz(idx).*VG.Ixz(idx) + VG.Iyz(idx).^3 - VG.Iyy(idx).*VG.Izz(idx).*VG.Iyz(idx))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - (VG.Iyz(idx).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).^2)./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) + ((VG.Ixy(idx).*VG.Ixz(idx) + VG.Iyy(idx).*VG.Iyz(idx) + VG.Iyz(idx).*VG.Izz(idx)).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx));
VEvec1(:,2) = (VG.Ixz(idx).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).^2)./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - (VG.Ixz(idx).^3 + VG.Ixz(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Izz(idx).*VG.Ixz(idx) - VG.Ixy(idx).*VG.Izz(idx).*VG.Iyz(idx))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - ((VG.Ixx(idx).*VG.Ixz(idx) + VG.Ixy(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Izz(idx)).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx));
VEvec1(:,3) = ones(size(VG.Ixx(idx)));%!
% Regularization of Values
% Evec1(~isfinite(Evec1)) = 0; % all NaNs and Vnf must be equal to 0!
% Normalization of Evec1 to unit vector (divided by its norm) 
VEvec1norm = (sum((abs(VEvec1).^2),2)).^(1/2);
E1 = VEvec1./repmat(VEvec1norm,[1 3]);
E1ImagCheck = imag(E1);
VEvec1N = real(E1);
NaN1Check = isnan(VEvec1N(:,1));

% Eigenvector V of the second eigenvalue
VEvec2(:,1) = (VG.Ixz(idx).^2.*VG.Iyz(idx) - VG.Ixy(idx).*VG.Izz(idx).*VG.Ixz(idx) + VG.Iyz(idx).^3 - VG.Iyy(idx).*VG.Izz(idx).*VG.Iyz(idx))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - (VG.Iyz(idx).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) - (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*1i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2).^2)./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) + ((VG.Ixy(idx).*VG.Ixz(idx) + VG.Iyy(idx).*VG.Iyz(idx) + VG.Iyz(idx).*VG.Izz(idx)).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) - (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*1i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx));
VEvec2(:,2) = (VG.Ixz(idx).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) - (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*1i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2).^2)./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - (VG.Ixz(idx).^3 + VG.Ixz(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Izz(idx).*VG.Ixz(idx) - VG.Ixy(idx).*VG.Izz(idx).*VG.Iyz(idx))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - ((VG.Ixx(idx).*VG.Ixz(idx) + VG.Ixy(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Izz(idx)).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) - (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*1i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx));
VEvec2(:,3) = ones(size(VG.Ixx(idx)));%!
% Regularization of Values
% Evec2(~isfinite(Evec2)) = 0; % all NaNs and Vnf must be equal to 0!
% Normalization of Evec2 to unit vector (divided by its norm) 
VEvec2norm = (sum((abs(VEvec2).^2),2)).^(1/2);
E2 = VEvec2./repmat(VEvec2norm,[1 3]);
E2ImagCheck = imag(E2);
VEvec2N = real(E2);
NaN2Check = isnan(VEvec2N(:,1));

% Eigenvector V of the third eigenvalue -- This can be derived from the other two!
%VEvec3(:,1) = (VG.Ixz(idx).^2.*VG.Iyz(idx) - VG.Ixy(idx).*VG.Izz(idx).*VG.Ixz(idx) + VG.Iyz(idx).^3 - VG.Iyy(idx).*VG.Izz(idx).*VG.Iyz(idx))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - (VG.Iyz(idx).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) + (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2).^2)./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) + ((VG.Ixy(idx).*VG.Ixz(idx) + VG.Iyy(idx).*VG.Iyz(idx) + VG.Iyz(idx).*VG.Izz(idx)).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) + (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx));
%VEvec3(:,2) = (VG.Ixz(idx).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) + (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2).^2)./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - (VG.Ixz(idx).^3 + VG.Ixz(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Izz(idx).*VG.Ixz(idx) - VG.Ixy(idx).*VG.Izz(idx).*VG.Iyz(idx))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx)) - ((VG.Ixx(idx).*VG.Ixz(idx) + VG.Ixy(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Izz(idx)).*(VG.Ixx(idx)./3 + VG.Iyy(idx)./3 + VG.Izz(idx)./3 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./(2.*((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)) + (3.^(1./2).*(((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3)./((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3) - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)).*i)./2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + (((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^3./27 - (VG.Ixz(idx).^2.*VG.Iyy(idx))./2 - (VG.Ixy(idx).^2.*VG.Izz(idx))./2 - (VG.Ixx(idx).*VG.Iyz(idx).^2)./2 + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^2 - ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).^2./9 + VG.Ixy(idx).^2./3 + VG.Ixz(idx).^2./3 + VG.Iyz(idx).^2./3 - (VG.Ixx(idx).*VG.Iyy(idx))./3 - (VG.Ixx(idx).*VG.Izz(idx))./3 - (VG.Iyy(idx).*VG.Izz(idx))./3).^3).^(1./2) + ((VG.Ixx(idx) + VG.Iyy(idx) + VG.Izz(idx)).*(VG.Ixy(idx).^2 + VG.Ixz(idx).^2 + VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Iyy(idx) - VG.Ixx(idx).*VG.Izz(idx) - VG.Iyy(idx).*VG.Izz(idx)))./6 + VG.Ixy(idx).*VG.Ixz(idx).*VG.Iyz(idx) + (VG.Ixx(idx).*VG.Iyy(idx).*VG.Izz(idx))./2).^(1./3)./2))./(VG.Ixy(idx).*VG.Ixz(idx).^2 - VG.Ixy(idx).*VG.Iyz(idx).^2 - VG.Ixx(idx).*VG.Ixz(idx).*VG.Iyz(idx) + VG.Ixz(idx).*VG.Iyy(idx).*VG.Iyz(idx));
%VEvec3(:,3) = ones(size(VG.Ixx(idx)));%!
VEvec3 = cross(VEvec1N,VEvec2N);
% Regularization of Values
% Evec3(~isfinite(Evec3)) = 0; % all NaNs and Vnf must be equal to 0!
% Normalization of Evec3 to unit vector (divided by its norm) 
VEvec3norm = (sum((abs(VEvec3).^2),2)).^(1/2);
E3 = VEvec3./repmat(VEvec3norm,[1 3]);
E3ImagCheck = imag(E3);
VEvec3N = real(E3);
NaN3Check = isnan(VEvec3N(:,1));

% Due to machine error, I cannot fix the threshold to 0, however it would
% be something similar to 10^-6
ImagThreshold = 10^-6;

E1Cx = sum( abs(E1ImagCheck) >= ImagThreshold ,2) > 0;
E2Cx = sum( abs(E2ImagCheck) >= ImagThreshold ,2) > 0;
E3Cx = sum( abs(E3ImagCheck) >= ImagThreshold ,2) > 0;

idx_discard = E1Cx | NaN1Check | E2Cx | NaN2Check | E3Cx | NaN3Check;

if SortByMagnitude

    % Sorting Eigenvalues (by norm)
    [~,IDX] = sort(abs(VEval),2);
    
    L1l1 = IDX(:,1) == 1;
    L1l2 = IDX(:,2) == 1;
    L1l3 = IDX(:,3) == 1;
    
    L2l1 = IDX(:,1) == 2;
    L2l2 = IDX(:,2) == 2;
    L2l3 = IDX(:,3) == 2;
    
    L3l1 = IDX(:,1) == 3;
    L3l2 = IDX(:,2) == 3;
    L3l3 = IDX(:,3) == 3;
    
    
    % Direction Along the Vessel (minimum intensity variation)
    % EigenValues
    Evl1(L1l1,1) = VEval(L1l1,1);
    Evl1(L2l1,1) = VEval(L2l1,2);
    Evl1(L3l1,1) = VEval(L3l1,3);
    % EigenVectors
    Evc1(L1l1,:) = VEvec1N(L1l1,:);
    Evc1(L2l1,:) = VEvec2N(L2l1,:);
    Evc1(L3l1,:) = VEvec3N(L3l1,:);
    
    % First Orthogonal Base
    % EigenValues
    Evl2(L1l2,1) = VEval(L1l2,1);
    Evl2(L2l2,1) = VEval(L2l2,2);
    Evl2(L3l2,1) = VEval(L3l2,3);
    % EigenVectors
    Evc2(L1l2,:) = VEvec1N(L1l2,:);
    Evc2(L2l2,:) = VEvec2N(L2l2,:);
    Evc2(L3l2,:) = VEvec3N(L3l2,:);
    
    % Second Orthogonal Base
    % EigenValues
    Evl3(L1l3,1) = VEval(L1l3,1);
    Evl3(L2l3,1) = VEval(L2l3,2);
    Evl3(L3l3,1) = VEval(L3l3,3);
    %EigenVectors
    Evc3(L1l3,:) = VEvec1N(L1l3,:);
    Evc3(L2l3,:) = VEvec2N(L2l3,:);
    Evc3(L3l3,:) = VEvec3N(L3l3,:);

else
    
    Evl1 = VEval(:,1);
    Evl2 = VEval(:,2);
    Evl3 = VEval(:,3);
    
    Evc1 = VEvec1N;
    Evc2 = VEvec2N;
    Evc3 = VEvec3N;
    
end

function [A,CHK] = VTC3D_switchNaN2zero(A)
CHK = isnan(A);
A(CHK) = 0;
CHK = double(~CHK);

function [DwnSmplT11LE,DwnSmplT12LE,DwnSmplT13LE,DwnSmplT22LE,DwnSmplT23LE,DwnSmplT33LE] = VTC3D_convert3DTensorTo6LogEuclidComponents(El1DwnSmpl,El2DwnSmpl,El3DwnSmpl,Ev1DwnSmpl,Ev2DwnSmpl,Ev3DwnSmpl)
w = [1,sqrt(2),sqrt(2),1,sqrt(2),1];
    
DwnSmplT11LE =  w(1) * (...
                        Ev1DwnSmpl(:,:,:,1) .* log( El1DwnSmpl ) .* conj(Ev1DwnSmpl(:,:,:,1)) + ... e1x*LOG(l1)*conj(e1x) +  
                        Ev2DwnSmpl(:,:,:,1) .* log( El2DwnSmpl ) .* conj(Ev2DwnSmpl(:,:,:,1)) + ... e2x*LOG(l2)*conj(e2x) +
                        Ev3DwnSmpl(:,:,:,1) .* log( El3DwnSmpl ) .* conj(Ev3DwnSmpl(:,:,:,1))   ... e3x*LOG(l3)*conj(e3x)
                        ); % Equivalent to the Tensor Component (1,1)
    
DwnSmplT12LE =  w(2) * (... 
                        Ev1DwnSmpl(:,:,:,2) .* log( El1DwnSmpl ) .* conj(Ev1DwnSmpl(:,:,:,1)) + ... e1y*LOG(l1)*conj(e1x) + 
                        Ev2DwnSmpl(:,:,:,2) .* log( El2DwnSmpl ) .* conj(Ev2DwnSmpl(:,:,:,1)) + ... e2y*LOG(l2)*conj(e2x) +
                        Ev3DwnSmpl(:,:,:,2) .* log( El3DwnSmpl ) .* conj(Ev3DwnSmpl(:,:,:,1))   ... e3y*LOG(l3)*conj(e3x)
                        ); % Equivalent to the Tensor Component (1,2)
    
DwnSmplT13LE =  w(3) * (...
                        Ev1DwnSmpl(:,:,:,3) .* log( El1DwnSmpl ) .* conj(Ev1DwnSmpl(:,:,:,1)) + ... e1z*LOG(l1)*conj(e1x) + 
                        Ev2DwnSmpl(:,:,:,3) .* log( El2DwnSmpl ) .* conj(Ev2DwnSmpl(:,:,:,1)) + ... e2z*LOG(l2)*conj(e2x) + 
                        Ev3DwnSmpl(:,:,:,3) .* log( El3DwnSmpl ) .* conj(Ev3DwnSmpl(:,:,:,1))   ... e3z*LOG(l3)*conj(e3x)
                        ); % Equivalent to the Tensor Component (1,3)
    
DwnSmplT22LE =  w(4) * (...
                        Ev1DwnSmpl(:,:,:,2) .* log( El1DwnSmpl ) .* conj(Ev1DwnSmpl(:,:,:,2)) + ... e1y*LOG(l1)*conj(e1y) + 
                        Ev2DwnSmpl(:,:,:,2) .* log( El2DwnSmpl ) .* conj(Ev2DwnSmpl(:,:,:,2)) + ... e2y*LOG(l2)*conj(e2y) + 
                        Ev3DwnSmpl(:,:,:,2) .* log( El3DwnSmpl ) .* conj(Ev3DwnSmpl(:,:,:,2))   ... e3y*LOG(l3)*conj(e3y)
                        ); % Equivalent to the Tensor Component (2,2)
                   
DwnSmplT23LE =  w(5) * (...
                        Ev1DwnSmpl(:,:,:,3) .* log( El1DwnSmpl ) .* conj(Ev1DwnSmpl(:,:,:,2)) + ... e1z*LOG(l1)*conj(e1y) + 
                        Ev2DwnSmpl(:,:,:,3) .* log( El2DwnSmpl ) .* conj(Ev2DwnSmpl(:,:,:,2)) + ... e2z*LOG(l2)*conj(e2y) + 
                        Ev3DwnSmpl(:,:,:,3) .* log( El3DwnSmpl ) .* conj(Ev3DwnSmpl(:,:,:,2))   ... e3z*LOG(l3)*conj(e3y)
                        ); % Equivalent to the Tensor Component (2,3)
    
DwnSmplT33LE =  w(6) * (...
                        Ev1DwnSmpl(:,:,:,3) .* log( El1DwnSmpl ) .* conj(Ev1DwnSmpl(:,:,:,3)) + ... e1z*LOG(l1)*conj(e1z) + 
                        Ev2DwnSmpl(:,:,:,3) .* log( El2DwnSmpl ) .* conj(Ev2DwnSmpl(:,:,:,3)) + ... e2z*LOG(l2)*conj(e2z) + 
                        Ev3DwnSmpl(:,:,:,3) .* log( El3DwnSmpl ) .* conj(Ev3DwnSmpl(:,:,:,3))   ... e3z*LOG(l3)*conj(e3z)
                        ); % Equivalent to the Tensor Component (3,3)
    
clear El1DwnSmpl El2DwnSmpl El3DwnSmpl Ev1DwnSmpl Ev2DwnSmpl Ev3DwnSmpl;
 
function [El1,El2,El3,Ev1,Ev2,Ev3] = VTC3D_convert6LogEuclidComponentsTo3DTensor(avgT11LE,avgT12LE,avgT13LE,avgT22LE,avgT23LE,avgT33LE,Msk)

w = 1./[1,sqrt(2),sqrt(2),1,sqrt(2),1];

T.Ixx = w(1) * avgT11LE;
T.Ixy = w(2) * avgT12LE;
T.Ixz = w(3) * avgT13LE;
T.Iyy = w(4) * avgT22LE;
T.Iyz = w(5) * avgT23LE;
T.Izz = w(6) * avgT33LE;
%                                                                               ***********
[El1,El2,El3,Ev1,Ev2,Ev3, ~ ] = VTC3D_get3DImageEIGs( T , Msk , false ); % These MUST !!! NOT !!! BE Sorted By Magnitude! (Those are still in the Log-Euclidean Space! and the order is preserved! - monotonic - )
%                                                                               ***********
El1 = exp(El1); El1(~Msk) = NaN;
El2 = exp(El2); El2(~Msk) = NaN;
El3 = exp(El3); El3(~Msk) = NaN;

Ev1_1 = Ev1(:,:,:,1); Ev1_1(~Msk) = NaN;
Ev1_2 = Ev1(:,:,:,2); Ev1_2(~Msk) = NaN;
Ev1_3 = Ev1(:,:,:,3); Ev1_3(~Msk) = NaN;
Ev1 = cat(4,Ev1_1,Ev1_2,Ev1_3);
clear Ev1_*;

Ev2_1 = Ev2(:,:,:,1); Ev2_1(~Msk) = NaN;
Ev2_2 = Ev2(:,:,:,2); Ev2_2(~Msk) = NaN;
Ev2_3 = Ev2(:,:,:,3); Ev2_3(~Msk) = NaN;
Ev2 = cat(4,Ev2_1,Ev2_2,Ev2_3);
clear Ev2_*;

Ev3_1 = Ev3(:,:,:,1); Ev3_1(~Msk) = NaN;
Ev3_2 = Ev3(:,:,:,2); Ev3_2(~Msk) = NaN;
Ev3_3 = Ev3(:,:,:,3); Ev3_3(~Msk) = NaN;
Ev3 = cat(4,Ev3_1,Ev3_2,Ev3_3);
clear Ev3_*;

function TLIC3D = VTC3D_TLICFromTEIG_3D(EIG3D)

% Inverse of Eigenvalues: from flat disks to elongated blobs
iEl1 = 1./EIG3D.El1;
iEl2 = 1./EIG3D.El2;
iEl3 = 1./EIG3D.El3;

iElsScaleFactor = (iEl1.*iEl2.*iEl3).^(1/3);

iEl1 = iEl1./iElsScaleFactor;
iEl2 = iEl2./iElsScaleFactor;
iEl3 = iEl3./iElsScaleFactor;

TLIC3D(:,:,:,1) = (iEl1 .* EIG3D.Ev1(:,:,:,2).^2) + (iEl2 .* EIG3D.Ev2(:,:,:,2).^2) + (iEl3 .* EIG3D.Ev3(:,:,:,2).^2);

TLIC3D(:,:,:,2) = (iEl1 .* EIG3D.Ev1(:,:,:,2) .* EIG3D.Ev1(:,:,:,1)) + (iEl2 .* EIG3D.Ev2(:,:,:,2) .* EIG3D.Ev2(:,:,:,1)) + (iEl3 .* EIG3D.Ev3(:,:,:,2) .* EIG3D.Ev3(:,:,:,1));

TLIC3D(:,:,:,3) = (iEl1 .* EIG3D.Ev1(:,:,:,2) .* EIG3D.Ev1(:,:,:,3)) + (iEl2 .* EIG3D.Ev2(:,:,:,2) .* EIG3D.Ev2(:,:,:,3)) + (iEl3 .* EIG3D.Ev3(:,:,:,2) .* EIG3D.Ev3(:,:,:,3));

TLIC3D(:,:,:,4) = (iEl1 .* EIG3D.Ev1(:,:,:,1).^2) + (iEl2 .* EIG3D.Ev2(:,:,:,1).^2) + (iEl3 .* EIG3D.Ev3(:,:,:,1).^2);

TLIC3D(:,:,:,5) = (iEl1 .* EIG3D.Ev1(:,:,:,1) .* EIG3D.Ev1(:,:,:,3)) + (iEl2 .* EIG3D.Ev2(:,:,:,1) .* EIG3D.Ev2(:,:,:,3)) + (iEl3 .* EIG3D.Ev3(:,:,:,1) .* EIG3D.Ev3(:,:,:,3));

TLIC3D(:,:,:,6) = (iEl1 .* EIG3D.Ev1(:,:,:,3).^2) + (iEl2 .* EIG3D.Ev2(:,:,:,3).^2) + (iEl3 .* EIG3D.Ev3(:,:,:,3).^2);

function OldValueStr = VTC3D_showProgressCM(ProgressMessage,OldValueStr,Value,Finished)

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

function includeVTrailsLibs

addpath(strcat(VTRootDir,'libs/'));
addpath(strcat(VTRootDir,'libs/0_NIfTI_IO/'));
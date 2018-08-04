%%% VTrails Pipeline %%%
% Following here few hints on how to use the VTrails ToolKit.
% See /test_imgs/ for a DEMO with a couple of examples.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) [NECESSARY] INSTALL VTrails Toolkit
% >> INSTALL_VTrails

% NB. If you are experiencing compiling issues or installation problems,
% please configure the matlab coder, using the most recent gcc (LINUX) or
% clang (MAC) compiler.

%% 0a) [On Request] Create your own Dictionary of Filtering SLoGS Kernels
% >> VTrailsSLoGS3D_MAIN;

% NB. this step requires to further specify the desired Dictionary of
% Filtering Kernels to 'VTrailsFilter3D_MAIN'.
% NB. this function is available *on request* (not included in the
% open-source package).
% Tutorial -- https://youtu.be/f3LahhqzHFc

%% 1) Filter the Image with the SLoGS Filterbank
% >> scales_range = [0.1:0.1:1.0];
% >> JobDumpDirPath = VTrailsFilter3D_MAIN( 'path/ImgFileName.nii' , 'path/MaskFileName.nii' , scales_range );

% NB. Define the scale_range accordingly with the vascular information
% content of the image and with the desired level of detail.
% As rule of thumb: avoid big gaps from subsequent scales.
% NB. VTrailsFilter3D_MAIN imports automatically nifti files (.nii,.nii.gz)

%% 1a) Integrate the Multi-Resolution Filter Responses over Scales
% >> min_scale = min(scales_range);
% >> max_scale = max(scales_range);
% >> MaxMAP = VTF3D_integrateMaxMultiScaleResponses('PathFolder', JobDumpDirPath , 'ScalesRange', [min_scale , max_scale] );

% NB. select the range of scale-integration.

%% 1b) Determine the Organised Seeds for the Geodesic Connectivity Paradigm
% >> arbitraryQuantileThreshold = 0.75;
% >> verboseFLAG = true;
% >> MaxMAP.OS = VTF3D_AlignVesselSeeds3D( double(MaxMAP.CVMsmth) , MaxMAP.US , arbitraryQuantileThreshold , verboseFLAG );

% NB. the arbitraryQuantileThreshold is a scalar in [0,1], being 0 = all
% the points, and 1 = only the max value respectively.
% See refs.[1] and [2] for more insights on the arbitraryQuantileThreshold.

%% 1c) Sort the Organised Seeds as initial Connected Components
% >> [MaxMAP.SortedOS , MaxMAP.ConnectedSegmentsOS , MaxMAP.BrachPointsOS] = VTF3D_SortConnectedComponents( MaxMAP.OS , verboseFLAG );

% NB. Sort and organise the Seeds and Connected Components prior to the
% exhaustive Connectivity Paradigm.

%% 2) Determine the Over-Connected Geodesic Vascular Graph
% >> DATA = VTrailsConnectingGeodesicGraph3D_MAIN( MaxMAP , 'PathFolder' , JobDumpDirPath );

% NB. Connect the initial Seeds with the exhaustive Connectivity Paradigm,
% by following an Anisotropic Level-Set over the Riemannian Vesselness
% potential (CVM, TF).

%% 3) Extract the Vascular Tree(s) as Geodesic Minimum Spanning Tree(s) - optional Pruning and Refinement
% >> [GeodesicMSTs,~] = VTrailsRefine3D_ConnectedGraph2VTrailsTree_MAIN( DATA , 'PruningLengthThresholdMM' , 10 );

% NB. Determine the Acyclic topology as the geodesic Minimum Spanning Tree
% (MST) or a forest of geodesic MSTs, accordingly with the connected
% components in the vascular graph.

%% 4) Visualise Results
% >> VTA3D_visualiseRiemannianVesselness( MaxMAP.CVM , MaxMAP.TF , MaxMAP.TFValidityMSK );
% >> VTA3D_visualiseConnectingGeodesicGraph( DATA.exploredminPath , DATA.GGM );
% >> VTA3D_visualiseGeodesicMinimumSpanningTree( GeodesicMSTs , MaxMAP.CVM );

% NB. Display Results

%% REFERENCES
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
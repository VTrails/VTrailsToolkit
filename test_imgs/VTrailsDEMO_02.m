function VTrailsDEMO_02()

addpath('../main/');

imgFileName = './test_img02.nii.gz';

%% 1) Filter the Image with the SLoGS Filterbank
scales_range = 0.3:0.1:1.0;
JobDumpDirPath = VTrailsFilter3D_MAIN( imgFileName , [] , scales_range , 'SIQthr' , 0.75 );
fprintf('---------------------------------------------------------------------------------------\n\n');

%% 1a) Integrate the Multi-Resolution Filter Responses over Scales
min_scale = min(scales_range);
max_scale = max(scales_range);
MaxMAP = VTF3D_IntegrateMaxMultiScaleResponses('PathFolder', JobDumpDirPath , 'ScalesRange', [min_scale , max_scale]);
fprintf('---------------------------------------------------------------------------------\n\n');

%% 1b) Determine the Organised Seeds for the Geodesic Connectivity Paradigm
arbitraryQuantileThreshold = 0.1;
verboseFLAG = true;
MaxMAP.OS = VTF3D_AlignVesselSeeds3D( double(MaxMAP.CVMsmth) , MaxMAP.US , arbitraryQuantileThreshold , verboseFLAG );
fprintf('-----------------------------------------------\n\n');

%% 1c) Sort the Organised Seeds as initial Connected Components
[MaxMAP.SortedOS,...
 MaxMAP.ConnectedSegmentsOS,...
 MaxMAP.BranchPointsOS] = VTF3D_SortConnectedComponents( MaxMAP.OS , verboseFLAG );
fprintf('-----------------------------------------------\n\n');

%% 2) Determine the Over-Connected Geodesic Vascular Graph
DATA = VTrailsConnectingGeodesicGraph3D_MAIN( MaxMAP , 'PathFolder' , JobDumpDirPath , 'SPDPowerFactor' , 2);
fprintf('---------------------------------------------------------\n\n');

%% 3) Extract the Vascular Tree(s) as Geodesic Minimum Spanning Tree(s) - optional Pruning and Refinement
[GeodesicMSTs,~] = VTrailsRefine3D_ConnectedGraph2VTrailsTree_MAIN( DATA , 'PruningLengthThresholdMM' , 3.75 );
fprintf('---------------------------------------------------------\n\n');

%% 4) Visualise Results
VTA3D_visualiseRiemannianVesselness( MaxMAP.CVM , MaxMAP.TF , MaxMAP.TFValidityMSK );
VTA3D_visualiseConnectingGeodesicGraph( DATA.exploredminPath , DATA.GGM );
VTA3D_visualiseGeodesicMinimumSpanningTree( GeodesicMSTs , MaxMAP.CVM );
fprintf('---------------------------------------------------------\n');
disp('            *** VTrails: DEMO_02 - DONE! **** ');
fprintf('---------------------------------------------------------\n\n');
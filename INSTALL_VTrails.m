function INSTALL_VTrails()%
%
% INSTALL VTrails Toolkit in Matlab from the command shell:
% 
% >> INSTALL_VTrails
%
% This function defines first the VTrailsToolkitRoot directory and
% dependencies for the execution of other functions. Then it will compile
% and install all the necessary libraries and mex-files.
%
%                       [Tested on Matlab R2016b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VTrailsToolkitRootDir = pwd;

if exist(VTrailsToolkitRootDir,'dir')
    
    VTrailsToolkitRootDir = strcat(VTrailsToolkitRootDir,'/');
    
    [exportSuccess,ME1] = exportVTrailsToolkitRootDir(VTrailsToolkitRootDir);
    
    [installSuccess,ME2] = installVTrailsComponents(VTrailsToolkitRootDir);
    
    if exportSuccess && installSuccess
        disp('**************************************************************');
        disp('            VTrails Toolkit Succesfully Installed!             ');
        disp('**************************************************************');
    else
        fprintf(2,'\nVTrails Toolkit has NOT been correctly installed.\nPlease Check the Exceptions (ErrorLog.mat) and\nRe-run the installer prior to using the Tooltik.\n');
        exportErrorLog(VTrailsToolkitRootDir,cat(2,ME1,ME2));
    end
    
else
    fprintf(2,'\nVTrails Toolkit has NOT been installed.\nPlease Re-run the installer prior to using the Tooltik.');
end

function [exportSuccess,ME] = exportVTrailsToolkitRootDir(VTrailsToolkitRootDir)

ME = struct([]);

try
    disp('--------------------------------------------------------------');
    disp('Step: 0.0 - Exporting: VTRootDir: ...');
    FID = fopen( strcat(VTrailsToolkitRootDir,'main/VTRootDir.m') , 'w' );
    
    templateTXT = ['function VTrdir = VTRootDir() \n%%',...
        '%% This function returns the VTrails Toolkit Root Directory\n\n',...
        'VTrdir = %s; \n\n'];
    
    VTRootDir = sprintf('%s%s%s',sprintf( '\''' ),VTrailsToolkitRootDir,sprintf( '\''' ));
    
    fprintf(FID,templateTXT,VTRootDir);
    
    fclose(FID);
    
    exportSuccess = true;
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch ME
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    exportSuccess = false;
end

function [installSuccess,ME] = installVTrailsComponents(VTrailsToolkitRootDir)

%% Compiling Filtering Mex-files

installSuccess = true;

ME = struct([]);

% Step: 1.1
try
    disp('Step: 1.1 - Compiling: VTF3D_computeGradients2ndOrd: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/VTF3D_computeGradients2ndOrd.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 1.2
try
    disp('Step: 1.2 - Compiling: VTF3D_configTFSynthesis: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/VTF3D_configTFSynthesis.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 1.3
try
    disp('Step: 1.3 - Compiling: VTF3D_eigenDecompBlock3D: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/VTF3D_eigenDecompBlock3D.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 1.4
try
    disp('Step: 1.4 - Compiling: VTF3D_synthTF: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/VTF3D_synthTF.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 1.5
try
    disp('Step: 1.5 - Installing: VTF3D_*: ...');
    movefile( strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/*mex*') , strcat(VTrailsToolkitRootDir,'libs/') );
    rmdir( strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/codegen/') , 's' );
    delete( strcat(VTrailsToolkitRootDir,'libs/1_ImageFilter_MEXCompiler/gcGuiReport.mat') );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

%% Compiling Connectivity Mex-files:

% Step: 2.1
try
    disp('Step: 2.1 - Compiling: mxAnisoDistanceTransform4Graph: ...');
    mex( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/mxAnisoDistanceTransform4Graph.cpp'), ...
        strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/AnisotropicFastMarching4Graph.cpp'), ...
        strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/UpdateNeighborhood4Graph.cpp'), ...
        strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/Auxillary4Graph.cpp'), ...
        '-outdir', strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/') ) ;
    fprintf('\nDONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    try
        mex( '-DMX_COMPACT_32' ,...
             strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/mxAnisoDistanceTransform4Graph.cpp'), ...
             strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/AnisotropicFastMarching4Graph.cpp'), ...
             strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/UpdateNeighborhood4Graph.cpp'), ...
             strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/Auxillary4Graph.cpp'), ...
             '-outdir', strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/') ) ;
        fprintf('\nDONE!\n']);
        disp('--------------------------------------------------------------');
    catch me
        ME = cat(2,ME,me);
        fprintf('\b\b\b\b');
        fprintf(2,'FAILED!\n');
        disp('--------------------------------------------------------------');
        installSuccess = false;
    end
end

% Step: 2.2
try
    disp('Step: 2.2 - Installing: mxAnisoDistanceTransform4Graph: ...');
    movefile( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/AnisotropicFastMarching4Graph/*mex*') , strcat(VTrailsToolkitRootDir,'libs/') );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 2.3
try
    disp('Step: 2.3 - Compiling: gradient3Dcpp: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Gradient3D/gradient3Dcpp.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 2.4
try
    disp('Step: 2.4 - Installing: gradient3Dcpp: ...');
    movefile( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Gradient3D/*mex*') , strcat(VTrailsToolkitRootDir,'libs/') );
    rmdir( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Gradient3D/codegen/') , 's' );
    delete( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Gradient3D/gcGuiReport.mat') );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 2.5
try
    disp('Step: 2.5 - Extracting: boost_1_49_0.zip: ...');
    system( ['zip -s0 ',...
             strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/boost_1_49_0_split.zip'),...
             ' --out ',...
             strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/boost_1_49_0.zip'),...
             ' -q'] );
    unzip( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/boost_1_49_0.zip') , ...
           strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/') );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 2.6
try
    disp('Step: 2.6 - Extracting: rncarpio-linterp-ca556a0.zip: ...');
    unzip( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/rncarpio-linterp-ca556a0.zip') , ...
           strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/') );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 2.7
try
    disp('Step: 2.7 - Compiling: linterp_matlab: ...');
    mex( strcat('-I',VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/boost_1_49_0/'),...
        strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/rncarpio-linterp-ca556a0/src/linterp_matlab.cpp'),...
        '-outdir', strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/rncarpio-linterp-ca556a0/src/') );
    fprintf('\nDONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 2.8
try
    disp('Step: 2.8 - Installing: linterp_matlab: ...');
    movefile( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/rncarpio-linterp-ca556a0/src/linterp_matlab.mex*') , strcat(VTrailsToolkitRootDir,'libs/') );
    delete( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/boost_1_49_0.zip') );
    rmdir( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/boost_1_49_0/') , 's' );
    rmdir( strcat(VTrailsToolkitRootDir,'libs/2_GeodesicConnect_MEXCompiler/Linterp/rncarpio-linterp-ca556a0/') , 's' );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

%% Compiling Tree-Extraction and Refinement Mex-files:

% Step: 3.1
try
    disp('Step: 3.1 - Compiling: VTR3D_getPairs2BeEvaluated: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/VTR3D_getPairs2BeEvaluated.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 3.2
try
    disp('Step: 3.2 - Compiling: VTR3D_notUnique: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/VTR3D_notUnique.prj') );
    fprintf([repmat('\b',[1,195]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 3.3
try
    disp('Step: 3.3 - Compiling: VTR3D_repeatedValue: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/VTR3D_repeatedValue.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 3.4
try
    disp('Step: 3.4 - Compiling: VTR3D_repeatedValueSymm: ...');
    coder('-build', strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/VTR3D_repeatedValueSymm.prj') );
    fprintf([repmat('\b',[1,45]),' DONE!\n']);
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

% Step: 3.5
try
    disp('Step: 3.5 - Installing: VTR3D*: ...');
    movefile( strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/*mex*') , strcat(VTrailsToolkitRootDir,'libs/') );
    rmdir( strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/codegen/') , 's' );
    delete( strcat(VTrailsToolkitRootDir,'libs/3_RefineTree_MEXCompiler/gcGuiReport.mat') );
    fprintf('\b\b\b\b DONE!\n');
    disp('--------------------------------------------------------------');
catch me
    ME = cat(2,ME,me);
    fprintf('\b\b\b\b');
    fprintf(2,'FAILED!\n');
    disp('--------------------------------------------------------------');
    installSuccess = false;
end

function exportErrorLog(VTrailsToolkitRootDir,MExcept)

save( strcat(VTrailsToolkitRootDir,'VTrails_Install_ErrorLog.mat') , 'MExcept' , '-v7.3' );

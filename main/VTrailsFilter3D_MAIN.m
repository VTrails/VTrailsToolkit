function JobDumpDirPath = VTrailsFilter3D_MAIN(imgFileNAME,mskFileNAME,Scales,varargin)%
%
%VTrails: Synthesis of 3D Connected Vesselness Map & Tensor Field via SLoGS
%
%%%%%%%%%
% Inputs:
% 
% - imgFileNAME: NIFTY filename of the 3D Angiographic Image Volume
%   NB: It is assumed BRIGHT Vessels on DARK Background!
% - mskFileName: NIFTY filename of the 3D bw mask for the Angiographic
%                Image Volume. It can be empty [] if none.
% - Scales: Scales range interval ]0,1]. E.g. [0.1 : 0.1 : 1.0]
%
%                               [OPTIONAL]
%
% - ['AutoDenoise']: Bool Option for Automatic Denoise
%                    (Mean Filter) - Default: False.
% - ['AutoDenoiseMedian']: Bool Option for Automatic Denoise
%                          (Median Filter) - Default: False.
% - ['BoneStrip']: Approximated Maximal Intensity Value for Blood Vessels,
%                   when Bones are present in the image (e.g. CT Angio.)
%                   This allows an initial Bone-masking or Skull-Stripping.
%                   Estimate [max] relative to the orig image scale range.
%                   Default: [];
% - [UpdateDFKsMO]: Boolean flag to Enable/Disable Update of the Dictionary
%                   of Filtering Kernels already aligned to previous Main
%                   Orientations (DFKsMO). Default: False.
% - [SIQthr]: Seeds Intensity Quantile threshold -- Default: 0.9
% - [SIQthrFIX]: Scalar boolean flag to disable progressive increase of
%                SIQthr at higher resolutions. Default: False.
% - [PowerFactor]: Scalar real value to enhance the original image
%                  dynamic contrast following the law: 
%                  Img_enhanced = Img.^PowerFactor;
%                  Default: 1.0
% - ['SkipTFflag']: Bool Option for Skipping the Tensor Field Synthesis -
%                   Default: False.
% - ['Dictionary']: Complete filename of the User-Defined Disctionary of
%                   SLoGS to be used. - 
%                   Default: '<VTrailsROOTDir>/libs/DFKs_MAIN.mat'
% - [FullDump]: Scalar boolean flag to export the complete set of outputs.
%               Default: False.
%
%%%%%%%%%%
% Outputs:
%   For each selected Scale, the method will export a .mat file in a
%   dedicated folder located in the function directory, which contains the
%   'MAP' structure listing the following fields:
%
% - GRID: Structure with Image Grid and Scales' infos for Resampling.
% - IMG: Header of the input nifty Images.
% - CVM: Scalar Connected Vesselness Map obtained with SLoGS
% - BDM: Scalar Vessel Boundaries Map obtained with \deltaSLoGS
% - BGM: Scalar Vessel Background Map obtained with \nuSLoGS
% - TFLE: Vascular Tensors Field synthesized with SLoGS in the
%         Log-Euclidean Domain
% - [TF]: Vascular Tensors Field synthesized with SLoGS in the
%         Euclidean Domain -- optional with 'FullDump'
% - SSsIDX: list of Un-Organised Seeds TO BE Aligned to the vessels with
%           'VTF3D_AlignVesselSeeds3D.m' (after Multi-Scale Integration.) 
% - US: Logical Volume of Un-organised Seeds
%
%%%%%%%%%%%
% Example Call:
%
%    JobDumpDirPath = VTrailsFilter3D_MAIN('test_img.nii',[],[0.1:0.1:1.0]);
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

OPT = VTF3D_configureINPUTS(imgFileNAME,mskFileNAME,Scales,varargin);

[DATA,MainOrients,DFKsMO] = VTF3D_configureDATA(imgFileNAME,mskFileNAME,Scales,OPT);

for s = 1 : length(DATA.Scales)
    
    tic;
    
        [JobDumpDirPath,MainOrients,DFKsMO] = VTF3D_processImageAtEachScale(DATA,MainOrients,DFKsMO,s);

    toc;
    
end
disp('----------------------------------------------------------');

% Procedure to Merge the MainOrients and DFKsMO
VTF3D_proceedWithMergingMODFKsMO(DATA);

disp(strcat('Data Succesfully Exported in: ',JobDumpDirPath));
disp('');
disp('VTrails: Synthesis of 3D Connected Vesselness Map & Tensor Field via SLoGS -- COMPLETE!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OPT = VTF3D_configureINPUTS(imgFileNAME,mskFileNAME,Scales,optInputs)
% Function to Assert Correctness of INPUTS

if exist(imgFileNAME,'file') ~= 2
    error('Parsed input imgFileNAME does NOT exist!');
end

if ~isempty(mskFileNAME) && exist(mskFileNAME,'file')~=2
    error('Parsed input mskFileNAME does NOT exist!');
end

assert(isnumeric(Scales) && min(Scales(:)) > 0 && max(Scales(:)) <= 1 ,'Scales MUST be a Vector (or single scalar value) with range ]0,1].');

% Initializing OPT with Default Values
OPT.AutoDenoise = false;
OPT.AutoDenoiseMedian = false;
OPT.BoneStrip = [];
OPT.UpdateDFKsMO = false;
OPT.SIQthr = 0.9;
OPT.SIQthrFIX = false;
OPT.PowerFactor = 1.0;
OPT.FullDump = false;
OPT.Dictionary = strcat(VTRootDir,'libs/DFKs_MAIN.mat');
OPT.SkipTFflag = false;

if ~isempty(optInputs)
    for j = 1 : 2 :length(optInputs)
        
        switch upper(optInputs{j})
            case 'AUTODENOISE'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given AutoDenoise flag: -- Default: "false" applied!');
                else
                    OPT.AutoDenoise = optInputs{j+1};
                end
            case 'AUTODENOISEMEDIAN'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given AutoDenoiseMedian flag: -- Default: "false" applied!');
                else
                    OPT.AutoDenoise = optInputs{j+1};
                end
            case 'BONESTRIP'
                if isempty(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || sum(isfinite(optInputs{j+1}))~=numel(optInputs{j+1})
                    disp('Invalid|Not Given BoneStrip (Vascular Intensity Values [max]) array: -- Default: [] ("none") applied!');
                else
                    OPT.BoneStrip = optInputs{j+1};
                end
            case 'UPDATEDFKSMO'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given UpdateDFKsMO flag: -- Default: "false" applied!');
                else
                    OPT.UpdateDFKsMO = optInputs{j+1};
                end
            case 'SIQTHRFIX'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given SIQthrFIX flag: -- Default: "false" applied!');
                else
                    OPT.SIQthrFIX = optInputs{j+1};
                end
            case 'SIQTHR'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || sum(isfinite(optInputs{j+1}))~=numel(optInputs{j+1})
                    disp('Invalid|Not Given SIQthr (Saliency Intensity Quantile Threshold) value: -- Default: 0.9 applied!');
                else
                    OPT.SIQthr = optInputs{j+1};
                end
            case 'POWERFACTOR'
                if isempty(optInputs{j+1}) || ~isnumeric(optInputs{j+1}) || ~isscalar(optInputs{j+1}) || sum(isfinite(optInputs{j+1}))~=numel(optInputs{j+1})
                    disp('Invalid|Not Given PowerFactor for input Image value: -- Default: 1.0 applied!');
                else
                    OPT.PowerFactor = optInputs{j+1};
                end
            case 'FULLDUMP'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given FullDump : -- Default: "false" applied!');
                else
                    OPT.FullDump = optInputs{j+1};
                end
            case 'SKIPTFFLAG'
                if isempty(optInputs{j+1}) || ~islogical(optInputs{j+1}) || ~isscalar(optInputs{j+1})
                    disp('Invalid|Not Given SkipTFflag : -- Default: "false" applied!');
                else
                    OPT.SkipTFflag = optInputs{j+1};
                end
            case 'DICTIONARY'
                if isempty(optInputs{j+1})
                    disp('Invalid|Not Given Dictionary Path for VTRAILS: -- Default Dictionary: Enabled ');
                else
                    if exist(optInputs{j+1},'file') ~= 7
                        OPT.Dictionary = optInputs{j+1};
                    else
                        disp('Invalid|Not Given Dictionary Path for VTRAILS: -- Default Dictionary: Enabled ');
                        OPT.Dictionary = strcat(VTRootDir,'libs/DFKs_MAIN.mat');
                    end
                end
            otherwise
                disp('Given Unknown parameter(s) for OPT. DEFAULT value(s) will be set.');
        end
        
    end
end

function [DATA,MainOrients,DFKsMO] = VTF3D_configureDATA(imgFileNAME,mskFileNAME,Scales,OPT)

% Display Title Of SW
disp('***** VTrails: Synthesis of 3D Connected Vesselness Map & Tensor Field via SLoGS *****')
disp(strcat( datestr(datetime),' @ ', computer ));

DATA.job = strrep(num2str(now),'.','-');

[DATA.imgPATHexp,DATA.imgNAMEexp,DATA.imgEXTexp] = fileparts(imgFileNAME);

DATA.Scales = Scales;

% Loading Image(s)
DATA.IMG = load_nii(imgFileNAME);

if ~isempty(mskFileNAME)
    MSK = load_nii(mskFileNAME);
    DATA.ROI = (MSK.img > 0);
    [DATA.mskPATHexp,DATA.mskNAMEexp,DATA.mskEXTexp] = fileparts(mskFileNAME);
else
    DATA.ROI = true(size(DATA.IMG.img));
    DATA.mskPATHexp = '';
    DATA.mskNAMEexp = '';
    DATA.mskEXTexp = '';
end

I = DATA.IMG.img;
I = double(I);

if ~isempty(OPT.BoneStrip)
    % Applying the SkullStripping by considering the intensity range of
    % vessels of the CTAngio. 
    [I,DATA.BoneMSK] = VTF3D_skullStripCTA(I,OPT.BoneStrip);
    DATA.ROI = DATA.ROI & ~DATA.BoneMSK;
else
    DATA.BoneMSK = [];
end

% Retrieving VoxelSize
voxelsize = DATA.IMG.hdr.dime.pixdim(2:4); 
DATA.physicalvoxelsize = voxelsize;% in mm
% Relative Voxelsize
DATA.voxelsize = (voxelsize./(min(voxelsize)))';
% Full Width at Half Maximum for Gaussian Distribution of sigma = 1
DATA.FWHM = 2*sqrt(2*log(2));

% Intensity Range Re-Scale: [0-1]
I = (I-min(I(:)))/(max(I(:))-min(I(:)));

DATA.VoxelFactor = ((size(I)-1)./size(I))';

DATA.reductionfactor = repmat((1+(1-DATA.VoxelFactor)),1,size(Scales,2))./repmat(Scales,size(DATA.VoxelFactor,1),1);

if OPT.AutoDenoise
    % Optional Median Filter
    if OPT.AutoDenoiseMedian
        % Applying Median Filter 3D
        I = VTF3D_ordfilt3D(I,14,'replicate');% MEDIAN FILTER in 3D to remove possible Salt&Pepper
    end
    
    % Automatic estimation of Gaussian Noise STD for AutoDenoise
    [~,noisesigma] = VTF3D_estimate3DImageGaussianNoise(I);

    % Applying Mean Filter 3D (according to the estimated noisesigma or with minimal stop-band)
    if ~isempty(noisesigma)
        I = VTF3D_MeanFilter3D(I,noisesigma./DATA.voxelsize);
    else
        I = VTF3D_MeanFilter3D(I,diag( ( diag( (DATA.reductionfactor(:,end) .* DATA.voxelsize ) ./ DATA.FWHM).^2 - diag(DATA.voxelsize ./ DATA.FWHM).^2 ).^(1/2) ) );
    end
end

DATA.I = double(double(I).^OPT.PowerFactor);

% Theta in deg for Orientation Clustering
DATA.theta_threshold = 15;

% Load The NEW Matching Template Database
% Change here the Dictionary of Filtering Kernels in case of a new set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if logical(exist(OPT.Dictionary,'file'))
    load(OPT.Dictionary); % WARNING: variable 'DFKs' is loaded!   Do NOT OVERWRITE!
else
    error('SLoGS Dictrionary Of Filtering Kernels (DFKs): NOT FOUND!');
end

if logical(exist('DFKsMO_MAIN.mat','file'))
    load DFKsMO_MAIN.mat; % WARNING: variables 'DFKsMO' and 'MainOrients' are loaded!   Do NOT OVERWRITE!
    
    if ~isequal(SETtheta_threshold,DATA.theta_threshold)
        disp('Available Oriented SLoGS DFKs: Different Theta Threshold! -- Creating New Set');
        DFKsMO = struct([]);
        MainOrients = [];
    end
    
else
    disp('Available Oriented SLoGS DFKs: NOT FOUND! -- Creating New Set');
    DFKsMO = struct([]);
    MainOrients = [];
end

DATA.DFKs = DFKs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Z] = ndgrid(1:size(I,1),1:size(I,2),1:size(I,3));
DATA.X = X;
DATA.Y = Y;
DATA.Z = Z;

DATA.OPT = OPT;

% Configuring the DFKsticks:
% It is ASSUMED that the 1st DFK (i.e. DFKs(1)) is a Gaussian STICK.
% It will be oriented along several MAIN ORIENTATIONS in the 3D space, by
% sampling the unit sphere with an icosphere (2-level subdivision). 
% The DFKsticks will be then used to reconstruct an INITIAL Tubularity Map,
% by filtering the Image with the 3D-oriented linear tubular probes (DFKsticks).
IcosphereSubdivisionLvL = 2;

DATA.DFKOrientedSticks = VTF3D_DFKOrientedSticksFromIcosphere(DATA.DFKs(1),IcosphereSubdivisionLvL); % 1st Kernel is a Gaussian STICK!

function [TempDirPath,MainOrients,DFKsMO] = VTF3D_processImageAtEachScale(DATA,MainOrients,DFKsMO,s)

disp('--------------------------------------------------------------');
disp(['Scale: ',num2str(DATA.Scales(s),'%.2f'), ' -> ',...
    'Size: [ ',num2str(round(size(DATA.I,1)*DATA.Scales(s)),'%d'),' , ',num2str(round(size(DATA.I,2)*DATA.Scales(s)),'%d'),' , ',num2str(round(size(DATA.I,3)*DATA.Scales(s)),'%d'),' ]']);
% Get Saliencies (Seeds' Locations and Orientations) from SLoGS_tube
[Idwn,ED,DwnsSeeds,SSs,GRID] = VTF3D_get3DTubeSaliencies(DATA,s);

% Configure the SLoGS DFKs (Steering according to the Orientations)
[MainOrients,DFKsMO,DFKsMOsel,DFKsFlat] = VTF3D_configureSLoGS(DATA,ED,DwnsSeeds,MainOrients,DFKsMO);

% Compute the Scalar and Tensorial Maps from SLoGS Filtering
MAP = VTF3D_getCVMandTFwithSLoGS(DATA,Idwn,DFKsMOsel,DFKsFlat,GRID,s);

% Update the Seeds
MAP = VTF3D_updateSeeds(DATA,SSs,GRID,MAP);

% Exporting the Generated Maps as Temporary Variables for each Scale AND
% Exporting the Temporary set or MainOrients and DFKsMO to be merged as one
% single Dataset once the algorithm finishes for the specific image.
TempDirPath = VTF3D_exportTempResults(DATA,MAP,MainOrients,DFKsMO);

function [Idwn,ED,DwnsSeeds,SSs,GRID] = VTF3D_get3DTubeSaliencies(DATA,s)

% Determining the Gaussian Blur sigma relative to the scale
GRID.sigma = diag( ( diag( (DATA.reductionfactor(:,s) .* DATA.voxelsize ) ./ DATA.FWHM).^2 - diag(DATA.voxelsize ./ DATA.FWHM).^2 ).^(1/2) );

GRID.Xv = linspace(1,size(DATA.I,1),round(size(DATA.I,1)*DATA.Scales(s)));
GRID.Yv = linspace(1,size(DATA.I,2),round(size(DATA.I,2)*DATA.Scales(s)));
GRID.Zv = linspace(1,size(DATA.I,3),round(size(DATA.I,3)*DATA.Scales(s)));
[GRID.XR,GRID.YR,GRID.ZR] = ndgrid(GRID.Xv,GRID.Yv,GRID.Zv);

% Smoothing the Image with gaussian kernel (Anti-Aliasing)
Igb = VTF3D_GaussBlur3D(DATA.I,GRID.sigma); % Non-constant speed across scales, but accurate.
    
% DownSampling the Image as per the Scales' Scheme (Pyramid)
Idwn = interpn(DATA.X,DATA.Y,DATA.Z,Igb,GRID.XR,GRID.YR,GRID.ZR);
ROIdwn = imerode( interpn(DATA.X,DATA.Y,DATA.Z, imerode(DATA.ROI,true([3 3 3])) ,GRID.XR,GRID.YR,GRID.ZR,'nearest') , true([3 3 3]) ) ;

% Configure Splitting Pipeline for LARGE Images
IdwnMaxBloclSize = [256 256 256];
Idwnblk = VTF3D_divideImgInBlks(Idwn,IdwnMaxBloclSize);
ROIdwnblk = VTF3D_divideImgInBlks(ROIdwn,IdwnMaxBloclSize);

if length(Idwnblk) == 1
% Local Tubes Saliency Map
    clear Idwnblk ROIdwnblk;
    [~,ED,DwnsSeeds,SSs,GRID.SIQthr] = VTF3D_tubeSaliencySticks3D(Idwn,DATA.DFKOrientedSticks,ROIdwn,s,DATA.OPT.SIQthr,DATA.OPT.SIQthrFIX,DATA.Scales);
else
    ProgressMessage = 'Splitting Procedure for LARGE Images: Enabled --';
    ED = struct('El1',[],'El2',[],'El3',[],'Ev1',[],'Ev2',[],'Ev3',[]);
    DwnsSeeds = false(size(Idwn));
    SSs = struct('Sx',[],'Sy',[],'Sz',[]);
    SIQthrs = [];
    
    ProgressValue = 0;
    ProgressValue_step = 100.0/length(Idwnblk);
    ProgressOldString = VTF3D_showProgressCM(ProgressMessage,'',ProgressValue,false);
    
    for blk = 1 : length(Idwnblk)
        [~,ED_blk,DwnsSeeds_blk,SSs_blk,SIQthr_blk] = VTF3D_tubeSaliencySticks3D( Idwnblk(blk).img , DATA.DFKOrientedSticks , ROIdwnblk(blk).img , s , DATA.OPT.SIQthr , DATA.OPT.SIQthrFIX , DATA.Scales );
        
        ED.El1 = cat(1,ED.El1,ED_blk.El1);
        ED.El2 = cat(1,ED.El2,ED_blk.El2);
        ED.El3 = cat(1,ED.El3,ED_blk.El3);
        ED.Ev1 = cat(1,ED.Ev1,ED_blk.Ev1);
        ED.Ev2 = cat(1,ED.Ev2,ED_blk.Ev2);
        ED.Ev3 = cat(1,ED.Ev3,ED_blk.Ev3);
        
        DwnsSeeds( Idwnblk(blk).pos(1,1) : Idwnblk(blk).pos(1,2) , ...
                   Idwnblk(blk).pos(2,1) : Idwnblk(blk).pos(2,2) , ...
                   Idwnblk(blk).pos(3,1) : Idwnblk(blk).pos(3,2) ) = DwnsSeeds_blk;
        
        SSs.Sx = cat(1, SSs.Sx , SSs_blk.Sx + (Idwnblk(blk).pos(1,1) - 1) );
        SSs.Sy = cat(1, SSs.Sy , SSs_blk.Sy + (Idwnblk(blk).pos(2,1) - 1) );
        SSs.Sz = cat(1, SSs.Sz , SSs_blk.Sz + (Idwnblk(blk).pos(3,1) - 1) );
        
        SIQthrs = cat(1,SIQthrs,SIQthr_blk);
        
        % Progress Percentage
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTF3D_showProgressCM([],ProgressOldString,ProgressValue, blk == length(Idwnblk) );
    end
    GRID.SIQthr = mean(SIQthrs);
end

function [MainOrients,DFKsMO,DFKsMOsel,DFKsFlat] = VTF3D_configureSLoGS(DATA,ED,DwnsSeeds,MainOrients,DFKsMO)

% Estimating the Orientation from Seeds
[MainOrientsNEW,OrientsFlagVolume] = VTF3D_ApproxVesselMainOrientations(ED.Ev1,ED.Ev2,DwnsSeeds,DATA.theta_threshold);

[MainOrients,OrientsFlagVolume,MOstart] = VTF3D_UpdateMainOrientations(MainOrientsNEW,OrientsFlagVolume,MainOrients,DATA.theta_threshold);

DFKsMO = VTF3D_configDFKsOnMainOrients(DATA.DFKs,MainOrients,DFKsMO,MOstart);

% Selecting only the Useful Part of DFKsMO, i.e. DFKsMOsel
DFKsMOsel = retrieveSelectedDFKsMO(DFKsMO,OrientsFlagVolume,DATA.DFKs);

% Defining FLAT DFKs (Boundaries and Background)
DFKsFlat = VTF3D_configFlatDFKs;

function MAP = VTF3D_getCVMandTFwithSLoGS(DATA,Idwn,DFKsMOsel,DFKsFlat,GRID,s)

MAP.GRID.Scale = single(DATA.Scales(s));
MAP.GRID.Scales = single(DATA.Scales);
MAP.GRID.Ord = s;

MAP.IMG.hdr = DATA.IMG.hdr;

MAP.IMG.ImgPath = DATA.imgPATHexp;
MAP.IMG.ImgName = DATA.imgNAMEexp;
MAP.IMG.ImgExt = DATA.imgEXTexp;

MAP.IMG.MskPath = DATA.mskPATHexp;
MAP.IMG.MskName = DATA.mskNAMEexp;
MAP.IMG.MskExt = DATA.mskEXTexp;
    
[VSs,BDs,BGs,TFs,TFsLE] = VTF3D_computeAllMaps(Idwn,DFKsMOsel,DFKsFlat,DATA.OPT.FullDump,DATA.OPT.SkipTFflag);
    
MAP.CVM = single(VSs);
MAP.BDM = single(BDs);
MAP.BGM = single(BGs);
if DATA.OPT.FullDump
    MAP.TF = VTF3D_makeTFsingle(TFs);
end
MAP.TFLE = VTF3D_makeTFsLEsingle(TFsLE);

MAP.GRID.Xr = single(GRID.Xv);
MAP.GRID.Yr = single(GRID.Yv);
MAP.GRID.Zr = single(GRID.Zv);
MAP.GRID.Xo = single(size(DATA.I,1));
MAP.GRID.Yo = single(size(DATA.I,2));
MAP.GRID.Zo = single(size(DATA.I,3));
MAP.GRID.sigma = single(GRID.sigma);
MAP.GRID.physicalvoxelsize = single(DATA.physicalvoxelsize);
MAP.GRID.SIQthr = single(GRID.SIQthr);
MAP.GRID.SIQthrFIX = DATA.OPT.SIQthrFIX;
        
function MAP = VTF3D_updateSeeds(DATA,SSs,GRID,MAP)

% Organizing the Initial Skeleton Seeds
SSs.Sx(:,2) = round(GRID.Xv(SSs.Sx));
SSs.Sy(:,2) = round(GRID.Yv(SSs.Sy));
SSs.Sz(:,2) = round(GRID.Zv(SSs.Sz));
SSsIDX = sub2ind(size(DATA.I),SSs.Sx(:,2),SSs.Sy(:,2),SSs.Sz(:,2));

tempVolume = false(size(DATA.ROI));
tempVolume(SSsIDX) = true;
tempVolume = tempVolume & DATA.ROI;

[XO,YO,ZO] = ndgrid(1:size(DATA.I,1),1:size(DATA.I,2),1:size(DATA.I,3));
[XR,YR,ZR] = ndgrid(GRID.Xv,GRID.Yv,GRID.Zv);

tempVolume = interpn(XO,YO,ZO, tempVolume, XR,YR,ZR, 'nearest');

MAP.SSsIDX = single(find(tempVolume));

if DATA.OPT.FullDump
    MAP.US = false(size(MAP.CVM));
    MAP.US(MAP.SSsIDX) = true;
    MAP.OS = VTF3D_AlignVesselSeeds3D(double(MAP.CVM),MAP.US,0,false);
end
    
function TempDirName = VTF3D_exportTempResults(DATA,MAP,MainOrients,DFKsMO)
% Exporting Temporary MAP Dump and MainOrients and DFKsMO Dump for each scale

TempDirName = strcat([DATA.imgPATHexp,'/JobDump_',DATA.imgNAMEexp,'_',DATA.job,'/']);
if ~logical(exist(TempDirName,'dir'))
    mkdir(TempDirName);
end

TempMAPFileName = strcat([TempDirName,'TempDataset',num2str(MAP.GRID.Scale,'%.2f'),'.mat']);

%disp(['Exporting Temporary MAPs to: ',TempMAPFileName,' ...']);

save(TempMAPFileName,'MAP','-v7.3');

if DATA.OPT.UpdateDFKsMO
    TempMODFKsMOFileName = strcat([TempDirName,'TempMODFKsMO',num2str(MAP.GRID.Ord,'%03d'),'.mat']);
    SETtheta_threshold = DATA.theta_threshold;
    save(TempMODFKsMOFileName,'MainOrients','DFKsMO','SETtheta_threshold','-v7.3');
end

if ~isempty(DATA.OPT.BoneStrip)
    TempBONFileName = strcat([TempDirName,'BoneStrip.nii.gz']);
    if ~logical(exist(TempBONFileName,'file'))
        IMGbonestrip = DATA.IMG;
        IMGbonestrip.img = DATA.BoneMSK;
        save_nii(IMGbonestrip,TempBONFileName);
    end
end

function VTF3D_proceedWithMergingMODFKsMO(DATA)

if DATA.OPT.UpdateDFKsMO
    
    if logical(exist('../libs/DFKsMO_MAIN.mat','file'))
        load('../libs/DFKsMO_MAIN.mat'); % WARNING: variables 'DFKsMO', 'MainOrients' and 'SETtheta_threshold' are loaded!   Do NOT OVERWRITE!
        
        if ~isequal(SETtheta_threshold,DATA.theta_threshold)
            disp('The Existing Set of MainOrients and DFKsMO has a DIFFERENT theta_threshold! -- Creating and Merging a New Set');
            old_DFKsMO = struct([]);
            old_MainOrients = [];
            SETtheta_threshold = DATA.theta_threshold;
        else
            disp('The Existing Set of MainOrients and DFKsMO has the SAME theta_threshold! -- Merging data to Existing Set');
            old_DFKsMO = DFKsMO; clear DFKsMO;
            old_MainOrients = MainOrients; clear MainOrients;
        end
        
    else
        disp('Existing Set of MainOrients and SLoGS DFKsMO: NOT FOUND! -- Creating and Merging a New Set');
        old_DFKsMO = struct([]);
        old_MainOrients = [];
        SETtheta_threshold = DATA.theta_threshold;
    end
    
    TempDirName = strcat([DATA.imgPATHexp,'/JobDump_',DATA.imgNAMEexp,'_',DATA.job,'/']);
    MODFKsMOlist = dir([TempDirName,'TempMODFKsMO*.mat']);
    
    if ~isempty(MODFKsMOlist)
        disp('Loading and Merging Dataset:...');
        TOT_new_Entries = 0;
        
        for ll = 1 : size(MODFKsMOlist)
            
            load([TempDirName,MODFKsMOlist(ll).name]);
            
            [old_MainOrients,old_DFKsMO,newEntries] = VTF3D_mergeMainOrientsAndDFKsMO(old_MainOrients,MainOrients,old_DFKsMO,DFKsMO,SETtheta_threshold);
            
            TOT_new_Entries = TOT_new_Entries + newEntries;
            
            clear MainOrients DFKsMO;
        end
        
        disp(['Total New Entries: ',num2str(TOT_new_Entries,'%d')]);
        
        MainOrients = old_MainOrients;
        DFKsMO = old_DFKsMO;
        
        clear old_DFKsMO old_DFKsMO;
        
        disp('Exporting Merged Dataset to: "<VTrailsToolkitRoot>/libs/DFKsMO_MAIN.mat"');
        save('../libs/DFKsMO_MAIN.mat','MainOrients','DFKsMO','SETtheta_threshold');
        
        delete([TempDirName,'TempMODFKsMO*.mat']);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DFKOrientedSticks = VTF3D_DFKOrientedSticksFromIcosphere(stickDFK,IcoSubLvl)

PrincDirects = VTF3D_PrincipalDirectionsFromIcosphere(IcoSubLvl,true); % e.g.  n-Levels of Subdivision

DFKOrientedSticks = VTF3D_rotateDFKStickAlongPrincDirects(stickDFK,PrincDirects);

function PrincDirects = VTF3D_PrincipalDirectionsFromIcosphere(subdivisions,reduceFlag)

% Defining the Unit Sphere Directions for the Linear-Tubular Probe
[verts,~] = VTF3D_icosphere(subdivisions);

PrincDirects = struct([]);
newEntry = 0;

if reduceFlag

    excludeList = [];    
    
    for jj = 1 :size(verts,1)
        
        if isempty( VTF3D_notUnique(jj,excludeList) )
            
            firstfound = false;
            
            for kk = jj + 1 : size(verts,1)
                
                if (abs(sum(verts(jj,:).*verts(kk,:))) < 1e-6) && ~firstfound
                    newEntry = newEntry+1;
                    PrincDirects(newEntry).M = [verts(jj,:);verts(kk,:)];
                    firstfound = true;
                    
                elseif isequal(verts(jj,:), - verts(kk,:))
                    excludeList = cat(1,excludeList,kk);
                end
                
            end
            
        end
        
    end
    
else
    for jj = 1 :size(verts,1)
        for kk = 1 : size(verts,1)
            if (abs(sum(verts(jj,:).*verts(kk,:))) < 1e-6)
                newEntry = newEntry+1;
                PrincDirects(newEntry).M = [verts(jj,:);verts(kk,:)];
            end
        end
    end
end

function [vv,ff] = VTF3D_icosphere(varargin)
%ICOSPHERE Generate icosphere.
% Create a unit geodesic sphere created by subdividing a regular
% icosahedron with normalised vertices.
%
%   [V,F] = ICOSPHERE(N) generates to matrices containing vertex and face
%   data so that patch('Faces',F,'Vertices',V) produces a unit icosphere
%   with N subdivisions.
%
%   FV = ICOSPHERE(N) generates an FV structure for using with patch.
%
%   ICOSPHERE(N) and just ICOSPHERE display the icosphere as a patch on the
%   current axes and does not return anything.
%
%   ICOSPHERE uses N = 3.
%
%   ICOSPHERE(AX,...) plots into AX instead of GCA.
%
%   See also SPHERE.
%
%   Based on C# code by Andres Kahler
%   http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% Parse possible axes input
if nargin > 2
    error('Too many input variables, must be 0, 1 or 2.');
end
[cax,args,nargs] = axescheck(varargin{:});

n = 3; % default number of sub-divisions
if nargs > 0, n = args{1}; end % override based on input

% generate regular unit icosahedron (20 faced polyhedron)
[v,f] = VTF3D_icosphere_icosahedron(); % size(v) = [12,3]; size(f) = [20,3];

% recursively subdivide triangle faces
for gen = 1:n
    f_ = zeros(size(f,1)*4,3);
    for i = 1:size(f,1) % for each triangle
        tri = f(i,:);
        % calculate mid points (add new points to v)
        [a,v] = VTF3D_icosphere_getMidPoint(tri(1),tri(2),v);
        [b,v] = VTF3D_icosphere_getMidPoint(tri(2),tri(3),v);
        [c,v] = VTF3D_icosphere_getMidPoint(tri(3),tri(1),v);
        % generate new subdivision triangles
        nfc = [tri(1),a,c;
               tri(2),b,a;
               tri(3),c,b;
                    a,b,c];
        % replace triangle with subdivision
        idx = 4*(i-1)+1:4*i;
        f_(idx,:) = nfc;
    end
    f = f_; % update 
end

% remove duplicate vertices
[v,~,ix] = unique(v,'rows'); clear b % b dummy / compatibility
% reassign faces to trimmed vertex list and remove any duplicate faces
f = unique(ix(f),'rows');

switch(nargout)
    case 0 % no output
        cax = newplot(cax); % draw to given axis (or gca)
        VTF3D_icosphere_showSphere(cax,f,v);
    case 1 % return fv structure for patch
        vv = struct('Vertices',v,'Faces',f,...
                    'VertexNormals',v,'FaceVertexCData',v(:,3));
    case 2 % return vertices and faces
        vv = v; ff = f;
    otherwise
        error('Too many output variables, must be 0, 1 or 2.');
end

function [i,v] = VTF3D_icosphere_getMidPoint(t1,t2,v)
%GETMIDPOINT calculates point between two vertices
%   Calculate new vertex in sub-division and normalise to unit length
%   then find or add it to v and return index
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% get vertice positions
p1 = v(t1,:); p2 = v(t2,:);
% calculate mid point (on unit sphere)
pm = (p1 + p2) ./ 2;
pm = pm./norm(pm);
% add to vertices list, return index
i = size(v,1) + 1;
v = [v;pm];

function [v,f] = VTF3D_icosphere_icosahedron()
%ICOSAHEDRON creates unit regular icosahedron
%   Returns 12 vertex and 20 face values.
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
      1, t, 0; % v2
     -1,-t, 0; % v3
      1,-t, 0; % v4
      0,-1, t; % v5
      0, 1, t; % v6
      0,-1,-t; % v7
      0, 1,-t; % v8
      t, 0,-1; % v9
      t, 0, 1; % v10
     -t, 0,-1; % v11
     -t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));

% create faces
f = [ 1,12, 6; % f1
      1, 6, 2; % f2
      1, 2, 8; % f3
      1, 8,11; % f4
      1,11,12; % f5
      2, 6,10; % f6
      6,12, 5; % f7
     12,11, 3; % f8
     11, 8, 7; % f9
      8, 2, 9; % f10
      4,10, 5; % f11
      4, 5, 3; % f12
      4, 3, 7; % f13
      4, 7, 9; % f14
      4, 9,10; % f15
      5,10, 6; % f16
      3, 5,12; % f17
      7, 3,11; % f18
      9, 7, 8; % f19
     10, 9, 2];% f20

function VTF3D_icosphere_showSphere(cax,f,v)
%SHOWSPHERE displays icosphere given faces and vertices.
%   Displays a patch surface on the axes, cax, and sets the view.
%   
%   Properties:
%       - vertex normals == vertex vectors
%       - no backface lighting
%       - colour data matches z value of vertex
%       - material properties match default SURF
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% set some axes properties if not held
if ~ishold(cax)
    az = -37.5; el = 30;
    view(az,el)
    grid
end
% create patch object on cax
patch('Faces',f,'Vertices',v,...
    'VertexNormals',v,...
    'LineWidth',0.5,'FaceLighting','phong',...
    'BackFaceLighting','unlit',...
    'AmbientStrength',0.3,'DiffuseStrength',0.6,...
    'SpecularExponent',10,'SpecularStrength',0.9,...
    'FaceColor','flat','CData',v(:,3),...
    'Parent',cax,'Tag','Icosphere');

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

function DFKOrientedSticks = VTF3D_rotateDFKStickAlongPrincDirects(stickDFK,PrincDirects)

GRID = stickDFK.GRID;

DFKOrientedSticks = struct([]);

for ors = 1 : size(PrincDirects,2)
    
    %%% Rotating the Kernel Data according to the MainOrients %%%
    Dkrot = VTF3D_rotate3DKernel(stickDFK.Dk,stickDFK.phi,PrincDirects(ors).M); % SLOW but Working 
    
    %%% Cropping out the Padding %%%
    % Rotated Kernel Tubularity Probe 'Dkrot'
    Dkrot = Dkrot(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
    
    %%% Balancing the Kernel Responses %%%
    Dkrot = Dkrot - (sum(Dkrot(:))/numel(Dkrot)); % sum up to 0
    
    %%% Assignment and Storing Data %%%
    DFKOrientedSticks(ors,1).Dk = Dkrot;

    clear Dkrot;
end

function KernelRot = VTF3D_rotate3DKernel(Kernel,MainOrients,RotatedOrients)

% Determining the Rotation Matrix
R = VTF3D_RotMatrix4Vect(MainOrients(1,:),MainOrients(2,:),RotatedOrients(1,:),RotatedOrients(2,:));

% Rotating the Filtering Kernel
KernelRot = VTF3D_Rot3DwithR(Kernel,R,NaN);

function R = VTF3D_RotMatrix4Vect(x0,y0,x1,y1)
% This function estimates the Rotation Matrix to obtain the final alignment
% x0 to x1 and y0 to y1. x0, y0 and x1, y1 form both orthogonal bases
% respectively. (No need to align also the third unit vector as derived by
% the cross product between those parsed).
% All input are expected to be unit vectors and row vectors of size 1x3.
% The resulting Rotation Matrix will be 3x3.
tol = 1e-3;
if sum(abs(x0-x1)>tol) > 0 % if the inputs x0 and x1 are different
    if sum(abs(abs(x0)-abs(x1))>tol) > 0 % if the inputs' absolute values are different (they are not just opposite -- sign)
        
        vX = cross(x1,x0); % Cross product, taking into account the target FIRST, the source AFTER
        sX = sum(vX.^2).^(1/2); % Norm of the resulting vector (sin)
        cX = dot(x0,x1); % (cos) between vectors
        
        VXmat = [   0  , -vX(3)   vX(2) ;...
                  vX(3),    0  , -vX(1) ;...
                 -vX(2),  vX(1),   0   ];
        
        R1 = (eye(3) + VXmat + (VXmat^2)*((1-cX)/sX^2));
        
    else % in this case x0 = -x1
        R1 = eye(3);
        pos1 = find((x0.*x1) < 0);
        for jj = 1 : length(pos1)
            R1 (pos1(jj),pos1(jj)) = -1;
        end
        clear pos1;
    end
    
else
    R1 = eye(3);
end

y0r = y0*R1;

if sum(abs(y0r-y1)>tol) > 0 % if the inputs y0r and y1 are different
    if sum(abs(abs(y0r)-abs(y1))>tol) > 0 % if the inputs' absolute values are different (they are not just opposite -- sign)
        
        vY = cross(y1,y0r); % Cross product, taking into account the target FIRST, the source AFTER
        sY = sum(vY.^2).^(1/2); % Norm of the resulting vector (sin)
        cY = dot(y0r,y1); % (cos) between vectors
        
        VYmat = [   0  , -vY(3)   vY(2) ;...
                  vY(3),    0  , -vY(1) ;...
                 -vY(2),  vY(1),   0   ];
        
        R2 = (eye(3) + VYmat + (VYmat^2)*((1-cY)/sY^2));
        
    else % in this case y0r = -y1
        R2 = eye(3);
        pos2 = find((y0r.*y1) < 0);
        for kk = 1 : length(pos2)
            R2 (pos2(kk),pos2(kk)) = -1;
        end
        clear pos2;
    end
else
    R2 = eye(3);
end

R = R1*R2;

function Img3DRotated = VTF3D_Rot3DwithR(Img3D,R,ExtraValue)%
% Use this function to rotate a 3D image according to the Rotation Matrix R
% R is a 3x3 matrix.

% Translations
dx = 0;
dy = 0;
dz = 0;

% Dimensions
[nY, nX, nZ] = size(Img3D);

% Center Grids
[Y,X,Z] = ndgrid( ((1:nY)-round(nY/2)), ...
                  ((1:nX)-round(nX/2)), ...
                  ((1:nZ)-round(nZ/2)));

% Complete Transformation Matrix
T = eye(4);
T(1:3,1:3) = R;
T(1,4) = dx;
T(2,4) = dy;
T(3,4) = dz;

% Transformed Grids
Xtr = -T(1,4) + T(1,1)*X + T(1,2)*Y + T(1,3)*Z;
Ytr = -T(2,4) + T(2,1)*X + T(2,2)*Y + T(2,3)*Z;
Ztr = -T(3,4) + T(3,1)*X + T(3,2)*Y + T(3,3)*Z;

% Interpolation
Img3DRotated = interpn(Y,X,Z,Img3D,Ytr,Xtr,Ztr,'linear');  % More Accurate but Slower
Img3DRotated(isnan(Img3DRotated)) = ExtraValue;

function IBlur = VTF3D_GaussBlur3D(I,sigma)%
% sigma can be scalar, or a vector containing the diagonal of the gaussian
% covariance matrix:
%
%               s1 0 0
% sigma = diag( 0 s2 0 ) = [s1 s2 s3]
%               0 0 s3

if isscalar(sigma)
    sigma = double(sigma * [1 1 1]);
end

if (isvector(sigma) && sigma(1)~= 0 && sigma(2)~= 0 && sigma(3)~= 0)
    
    KernelSize = (round(3*sigma(:)) *2) +1;
    
    if ~isequal(KernelSize,[1 1 1]')
        [XGk,YGk,ZGk] =   ndgrid(-floor(KernelSize(1)/2):1:floor(KernelSize(1)/2),...
            -floor(KernelSize(2)/2):1:floor(KernelSize(2)/2),...
            -floor(KernelSize(3)/2):1:floor(KernelSize(3)/2));
        
        gain = 1/(sqrt((2*pi)^length(size(I)) * (prod(sigma.^2))));
        exponent = -( ((XGk).^2)/(2*sigma(1)^2) + ((YGk).^2)/(2*sigma(2)^2) + ((ZGk).^2)/(2*sigma(3)^2) );
        Gkernel = gain.*exp(exponent);
        
        Gkernel = Gkernel./sum(Gkernel(:));
    else
        sigma = 3*sigma;
        Gkernel = reshape(kron(kron([0 1 0]+sigma(1)*[1 1 1],[0 1 0]'+sigma(2)*[1 1 1]'),[0 1 0]'+sigma(3)*[1 1 1]'),[3 3 3]);
        Gkernel = Gkernel./sum(Gkernel(:));
    end
    
    IBlur = imfilter(double(I),double(Gkernel),'replicate','same');    
    
else
    IBlur = I; % No blurring applied if (one of) the sigma is equal to 0!
end

function IMean = VTF3D_MeanFilter3D(I,sigmaKernel)%
% sigmaKernel can be scalar, a vector or a matrix containing the filter

if isscalar(sigmaKernel)
    sigmaKernel = 3 * sigmaKernel;
    SKernel = reshape(kron(kron([0 1 0]+sigmaKernel*[1 1 1],[0 1 0]'+sigmaKernel*[1 1 1]'),[0 1 0]'+sigmaKernel*[1 1 1]'),[3 3 3]);
    SKernel = SKernel./sum(SKernel(:));
elseif isvector(sigmaKernel)
    sigmaKernel = 3 * sigmaKernel;
    SKernel = reshape(kron(kron([0 1 0]+sigmaKernel(1)*[1 1 1],[0 1 0]'+sigmaKernel(2)*[1 1 1]'),[0 1 0]'+sigmaKernel(3)*[1 1 1]'),[3 3 3]);
    SKernel = SKernel./sum(SKernel(:));
else
    SKernel = sigmaKernel;
end

if ~isequal(SKernel,zeros(size(SKernel)))
    
    IMean = imfilter(double(I),double(SKernel),'replicate','same');    
    
else
    IMean = I; % No blurring applied if SKernel is equal to 0!
end

function [Itube,ED,DwnsSeeds,SSs,quantile_thr] = VTF3D_tubeSaliencySticks3D(I,DFKOrientedSticks,ROI,scale,SIQthr,SIQthrFIX,Scales)%
% Tubularity Saliency Map:
% It is obtained by convolving the Image I with a set of kernels
% representig sticks (small tubes) oriented in as many different directions
% as the unit shpere is sampled by an icosphere (3-level subdivision).
%
%N.B.: the Saliency Map is not returned as it is needed only as
%      initialization and seeds identification.

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First Step: Computing Sticks Convolution %
[Itube,~] = VTF3D_filtScalarMap(I,DFKOrientedSticks,0,[0,1],[]);
ItubeM = Itube .* ROI;

%%% Second Step: Determining Seeds

% Computing Gradient and Hessian components for the EigenDecomposition
VG = VTF3D_computeGradients2ndOrd_mex(Itube);

%%% SINKS: Divergence of the Tubularity Map's Gradient
DIV = divergence( VG.Ix, VG.Iy , VG.Iz );
SINKS = (DIV < 0);

%%% HIGH-INTENSITY Tubular Structures
[f,bins] = hist(ItubeM(ItubeM > 0),1000);
f_cumulative = cumsum(f)/sum(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Threshold can affect dramatically the output of the matching
% template filtering.
% Please Tune it WISELY!
% Usually a value greater than 0.75 (even 0.95 should be OK!!!)
% Low quantile_thr => longer matching template processing & Noise Enhancement
% High quantile_thr => faster mathing template processing & Less details in the Connected Vesselness Map
if ~ SIQthrFIX
    increment = (1 - SIQthr)/numel(Scales);
    quantile_thr = SIQthr + (increment * (scale - 1)); % set value: ]0,1[
    quantile_thr(quantile_thr > 0.999) = 0.999;
else
    quantile_thr = SIQthr;
end

% NB: in future releases, this parameter might be scale-dependent.

% computing the respective tolerance interval using the quantile function
% for a normal distribution.
sigma_multiplier = sqrt(2)*erfinv(quantile_thr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

THR1 = bins(find(f_cumulative >= quantile_thr,1));

THR2 = mean(ItubeM(ItubeM > 0)) + sigma_multiplier * std(ItubeM(ItubeM > 0)); %!!! theshold set by standard deviation multiplier

THR = max([THR1,THR2]); % it could be a mean instead of max...

HIGHINT = (Itube >= THR);

%%% ATTRACTORS: EigenDecomposition
% The mask (SINKS & HIGHINT) is used to initialize the search (as the final
% intersection is considered). This would speed-up the process.
% Moreover the EigenDecomposition will be returned as Main Tubular
% Orientation.
[attEl1,attEl2,attEl3,~,~,~,attCHK] = VTF3D_get3DImageEIGs( VG , SINKS & HIGHINT , true , false); % sorting flag enabled - No absolute values! 

% Attractors of the Saliency Map
ATTRACT = ( attEl1 < 0 ) & ( attEl2 < 0 ) & ( attEl3 < 0 ) & attCHK;

% Final Intersection of logical Maps
DwnsSeeds = SINKS & ATTRACT & HIGHINT & ROI;

% Alignment and Reduction to Main Vessel Seeds %%% TEST %%%
DwnsSeeds = VTF3D_AlignVesselSeeds3D(Itube,DwnsSeeds,0,0); % EXTERNAL FUNCTION!
DwnsSeeds = DwnsSeeds & ROI;

[El1,El2,El3,Ev1,Ev2,Ev3,~] = VTF3D_get3DImageEIGs( VG , DwnsSeeds , true , false); % sorting flag enabled - No absolute values! 

% Linear Indexing of the Seeds
OrientidxReal = find(DwnsSeeds);

% Exporting the Linearized EigenDecomposition (ED) with form:
% length(OrientidxReal) x 1    EIGENVALUES
% length(OrientidxReal) x 3    EIGENVECTORS
ED.El1 = El1(OrientidxReal); % length(OrientidxReal) x 1
ED.El2 = El2(OrientidxReal); % length(OrientidxReal) x 1
ED.El3 = El3(OrientidxReal); % length(OrientidxReal) x 1
Ev1_1 = Ev1(:,:,:,1); Ev1_2 = Ev1(:,:,:,2); Ev1_3 = Ev1(:,:,:,3);
ED.Ev1 = cat(2,Ev1_1(OrientidxReal),Ev1_2(OrientidxReal),Ev1_3(OrientidxReal)); % length(OrientidxReal) x 3
clear Ev1_*;
Ev2_1 = Ev2(:,:,:,1); Ev2_2 = Ev2(:,:,:,2); Ev2_3 = Ev2(:,:,:,3);
ED.Ev2 = cat(2,Ev2_1(OrientidxReal),Ev2_2(OrientidxReal),Ev2_3(OrientidxReal)); % length(OrientidxReal) x 3
clear Ev2_*;
Ev3_1 = Ev3(:,:,:,1); Ev3_2 = Ev3(:,:,:,2); Ev3_3 = Ev3(:,:,:,3);
ED.Ev3 = cat(2,Ev3_1(OrientidxReal),Ev3_2(OrientidxReal),Ev3_3(OrientidxReal)); % length(OrientidxReal) x 3
clear Ev3_*;

% Exporting the Subscripts locations of the Seeds
[SSs.Sx,SSs.Sy,SSs.Sz] = ind2sub(size(DwnsSeeds),OrientidxReal);

function [Ifilt,IidxDFK] = VTF3D_filtScalarMap(I,DFK,basevalue,rescalerng,ProgressMessage)

maxblksize = [128 128 128]; % this can be edited to adapt computational performance!

Iblk = VTF3D_divideImgInBlks(I,maxblksize);

Ifilt = zeros(size(I));
IidxDFK = zeros(size(I));

if ~isempty(ProgressMessage)
    totitr = length(Iblk) * size(DFK,1);
    itr = 0;
    ProgressValue = 0;
    ProgressValue_step = 100.0/totitr;
    ProgressOldString = VTF3D_showProgressCM(ProgressMessage,'',ProgressValue,false);
    moditr = 50 + round(totitr/100);
end

for blk = 1 : length(Iblk)
    
    Iblkmax = zeros(size(Iblk(blk).img)) - eps;
    IblkmaxDFK = zeros(size(Iblk(blk).img));
    Iblkout = zeros(size(Iblk(blk).img));
    
    for dfk = 1 : size(DFK,1)
        
        
        
        % filter Response
        Iout = imfilter(Iblk(blk).img,DFK(dfk,1).Dk,'replicate','same');
        Iout(Iout < basevalue) = basevalue;
        
        % max filter response SLoGS indexing for tensors
        maxMSK = Iout > Iblkmax;
        Iblkmax(maxMSK) = Iout(maxMSK);
        IblkmaxDFK(maxMSK) = dfk;
        
        % integral scalar filter response
        Iblkout = Iblkout + Iout;
        
        if ~isempty(ProgressMessage)
            itr = itr + 1;
            % Progress Percentage
            ProgressValue = ProgressValue + ProgressValue_step;
            if mod(itr,moditr) == 1 || itr == totitr
                ProgressOldString = VTF3D_showProgressCM([],ProgressOldString,ProgressValue, itr == totitr );
            end
        end
        
    end
    
    Ifilt(Iblk(blk).pos(1,1):Iblk(blk).pos(1,2),Iblk(blk).pos(2,1):Iblk(blk).pos(2,2),Iblk(blk).pos(3,1):Iblk(blk).pos(3,2)) = Iblkout;
    IidxDFK(Iblk(blk).pos(1,1):Iblk(blk).pos(1,2),Iblk(blk).pos(2,1):Iblk(blk).pos(2,2),Iblk(blk).pos(3,1):Iblk(blk).pos(3,2)) = IblkmaxDFK;
    
end

% Re-scaling Range of values between 'rescalerng'.
if basevalue < 0
    assert(rescalerng(1) < 0,'The RescaleRange Values is not Negative!');
    % Normalizing symmetrically and linearly in rescalerng
    Ifilt(Ifilt<0) = rescalerng(1) * ( Ifilt(Ifilt<0) / min(Ifilt(Ifilt<0)) );
    Ifilt(Ifilt>0) = rescalerng(2) * ( Ifilt(Ifilt>0) / max(Ifilt(Ifilt>0)) );
else
    % Normalizing linearly in rescalerng
    Ifilt = rescalerng(1) +  ( (Ifilt - min(Ifilt(:)))/(max(Ifilt(:))-min(Ifilt(:))) ) * (rescalerng(2) - rescalerng(1));
end

function Iblk = VTF3D_divideImgInBlks(I,maxblocksize)

Size = size(I);

numblks = ceil(Size./maxblocksize);

entry = 0;

Iblk = struct([]);

for rr = 1 : numblks(1)
    for cc = 1 : numblks(2)
        for pp = 1 : numblks(3)
            
            entry = entry+1;
            
            rr_ini = (rr-1) * maxblocksize(1) + 1;
            rr_fin = (rr) * maxblocksize(1); rr_fin(rr_fin > Size(1)) = Size(1);
            cc_ini = (cc-1) * maxblocksize(2) + 1;
            cc_fin = (cc) * maxblocksize(2); cc_fin(cc_fin > Size(2)) = Size(2);
            pp_ini = (pp-1) * maxblocksize(3) + 1;
            pp_fin = (pp) * maxblocksize(3); pp_fin(pp_fin > Size(3)) = Size(3);

            Iblk(entry).pos =  [rr_ini,rr_fin ; cc_ini,cc_fin ; pp_ini,pp_fin];
            Iblk(entry).img = I(rr_ini:rr_fin , cc_ini:cc_fin , pp_ini:pp_fin);
            
        end
    end
end

function [El1,El2,El3,Ev1,Ev2,Ev3,CHK] = VTF3D_get3DImageEIGs(I,Msk,SortByMagnitude,AbsoluteFlag)
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
    IMG = VTF3D_computeGradients2ndOrd_mex(I);
else
    Size = size(I.Ixx);
    IMG = I;
end
%disp('Hessian: Done');

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

CHK = false(Size);

[VEval1NS, VEval2NS, VEval3NS, VEvec1NS, VEvec2NS, VEvec3NS, idxReal] = VTF3D_compute3DEIGs(IMG,idx,SortByMagnitude,AbsoluteFlag);

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

function [VEval1NS,VEval2NS,VEval3NS,VEvec1NS,VEvec2NS,VEvec3NS,idxReal] = VTF3D_compute3DEIGs(VG,idx,SortByMagnitude,AbsoluteFlag)

block_len = 16^3; %linearized length of a 16x16x16 cube of voxels
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
   
    %disp(['Sub-Block: ',num2mstr(j),'/',num2str(num_block)]);
    %tic;
    
    if j < num_block
       
       ini = ((j-1)*block_len)+1;
       fin = (j*block_len);
       
       idx_block = idx( ini : fin );
   else
       
       ini = ((j-1)*block_len)+1;
       fin = size(idx,1);
       
       idx_block = idx( ini : fin );
   end
   
   [Evl1 , Evl2 , Evl3 , Evc1 , Evc2 , Evc3 , idx_discard] = VTF3D_eigenDecompBlock3D_mex(VG,idx_block,SortByMagnitude,AbsoluteFlag);
   
   VEval1NS(ini:fin,1) = Evl1;
   VEval2NS(ini:fin,1) = Evl2;
   VEval3NS(ini:fin,1) = Evl3;
   
   VEvec1NS(ini:fin,:) = Evc1;
   VEvec2NS(ini:fin,:) = Evc2;
   VEvec3NS(ini:fin,:) = Evc3;
   
   idxRealCheck(ini:fin) = ~idx_discard; % list of indices where Eigenvectors have to be considered ( i.e. PURE REAL, so Complex and NaN values are excluded)
   
   clear ini fin idx_block Evc*;
   
   %toc;
   
end

VEvec1NS(~idxRealCheck,:) = NaN;
VEvec2NS(~idxRealCheck,:) = NaN;
VEvec3NS(~idxRealCheck,:) = NaN;
idxReal = idxRealCheck;

function [MainOrients,OrientsFlagVolume] = VTF3D_ApproxVesselMainOrientations(Evect1,Evect2,Seeds,theta_thrshld)
% The Main Orientation Matrix will be n x 3 x 2 (3components + 3components: Orthogonal BASIS in the third dimension)
% theta_thrshld is in deg.

MainOrients = [];
OrientsFlagVolume = zeros(size(Seeds));

SeedsIdx = find(Seeds);

% Control
if length(SeedsIdx) ~= size(Evect1,1)
    error('Mismatch in number of Seeds and associated Eigenvectors!');
end

% Define the 'Visited' Variable
Visited = zeros(size(Evect1,1),1);
Flag1 = zeros(size(Visited));

%wc1 = 0;

% Start a while cycle (till all the EigenVectors have been visited)
while sum(Visited) ~= length(Visited)
    
    %wc1 = wc1 + 1
    
    E1t = Evect1(~Visited,:);
    
    pos = find(~Visited);
    
    % First Level Clustering: Get a Centroid from the unvisited Evect1
    [CE1g,IDX1] = VTF3D_GetGroup(E1t,theta_thrshld);
    
    if isempty(IDX1)
        disp('IDX1 Empty!');
        error('There must be some Bug in the Evect1 Data!');
    end
    
    E2t = Evect2(IDX1,:);
    
    % Second Level Clustering: Given the Centroid, find all the other
    % orientations in Evect2
    [CE2g,Flag2] = VTF3D_ProcessEvects2(E2t,CE1g,IDX1,theta_thrshld);
    
    % Assignment of Values 
    Flag1(pos(IDX1)) = max(Flag1) + Flag2;
    
    C_Block = cat(3,repmat(CE1g,size(CE2g,1),1),CE2g);
    
    MainOrients = cat( 1 , MainOrients , C_Block );
    
    % Updating the Visited Evect1
    Visited(pos(IDX1)) = 1; % Updating Visited
    
    clear E1t CE1g IDX1 CE2g Flag2 C_Block

end

OrientsFlagVolume(SeedsIdx) = Flag1;

function PEvs = VTF3D_ProjectEv(Evs,E)%
% Evs is matrix of EigenVectors: N x dim (es: 7 x 3) in 3D 
% E is the EigenVector considered for projection 

PEvs = Evs.*repmat(E,size(Evs,1),1);

function BProject = VTF3D_BoolProjectEv(PEvs,theta_thrshld)%
% PEvs is matrix of Projected Eigenvectors
% theta_thrshld is the threshold angle

BProject = sum(PEvs,2) >= cos(pi/180*theta_thrshld);

function C = VTF3D_ComputeCentroid(Evs)%
% Evs is the matrix of Eigenvectors

C = mean(Evs,1);
% Normalize to unit Vector
C = VTF3D_UnitVector(C);

function En = VTF3D_UnitVector(E)%
% E is and EigenVector or matrix of Eigenvectors

Norm = (sum(E.^2,2)).^(1/2);
En = E./repmat(Norm,1,size(E,2));

function [CE1g,IDX] = VTF3D_GetGroup(E1t,theta_thrshld)%
% E1t is the matrix of unvisited Eigenvectors
% theta_thrshld is the threshold angle for makeing up a group

% Randomly pick one EigenVector from E1t.
Ridx = randi(size(E1t,1));
RE1sel = E1t(Ridx,:);

% Project all the non-visited EigenVectors in E1t to the selected one.
PE1t = VTF3D_ProjectEv(E1t,RE1sel);

% Evaluate the projection in terms of angle enclosed among the selected one
% and the projected ones - make use of cosine(theta_thrshld)
BPE1t = VTF3D_BoolProjectEv(PE1t,theta_thrshld);

%-- inner while cycle for Centroid Convergence
% initialization
converged = 0;

BPE1t_groupOLD = ones(sum(BPE1t),1);

while ~converged
    
    % All the Eigenvectors that satisfy the relationship: theta < theta_thrshld
    % are considered as a group (or improperly, cluster)
    E1t_group = E1t(BPE1t,:);
    
    % Given the group, the centroid must be computed
    CE1t_group = VTF3D_ComputeCentroid(E1t_group);
    
    if sum(isnan(CE1t_group)) ~= 0
        disp('Centroid NaN!')
    end
    
    % Again, as second check, the group of EigenVectors must be projected on
    % the centroid, and those that still satisfy the relationship will be part
    % of the final cluster.
    PE1t_group = VTF3D_ProjectEv(E1t_group,CE1t_group);
    BPE1t_group = VTF3D_BoolProjectEv(PE1t_group,theta_thrshld);
    
    % Those that do no more satisfy the second check are left out of the group,
    % and can be considered again for further clusters. (if (false) case)
    
    % Those who are part of the cluster, go through another clustering approach
    % taking into account the Evect2. (if (true) case)
    
    if isequal(BPE1t_group,BPE1t_groupOLD)
        converged = 1;
        
        CE1g = CE1t_group;
        IDX = find(BPE1t);
    else
        %disp('DropOut!');
        converged = 0;
        
        BPE1t(BPE1t == 1) = BPE1t_group;
        
        BPE1t_groupOLD = nonzeros(BPE1t_group);
        
        clear E1t_group CE1t_group PE1t_group BPE1t_group
    end
    
end

function [CE2g,Flag2] = VTF3D_ProcessEvects2(E2t,CE1g,IDX1,theta_thrshld)%
% Once found the first Cluster, proces in a similar fashion also the Evect2

% Precondition
PE2t = sum(VTF3D_ProjectEv(E2t,CE1g),2);
E2t = VTF3D_UnitVector(E2t - (repmat(PE2t,1,3).*repmat(CE1g,size(E2t,1),1)));

Visited2 = zeros(size(IDX1));
Flag2 = zeros(size(IDX1));

wc2 = 0;

while sum(Visited2) ~= length(Visited2)
    
    wc2 = wc2 + 1;
    
    E2temp = E2t(~Visited2,:);
    
    pos = find(~Visited2);
    
    [CE2g(wc2,:),IDX2] = VTF3D_GetGroup(E2temp,theta_thrshld);
    
    Flag2(pos(IDX2)) = wc2; % Storing the number of the Cluster
    
    Visited2(pos(IDX2)) = 1; % Updating Visited2
end

function [MainOrients,OrientsFlagVolume,MOstart] = VTF3D_UpdateMainOrientations(MainOrientsNEW,OrientsFlagVolume,MainOrients,theta_thrshld)
% Use this function to compare MainOrientsNEW to existing MainOrients.
% This will avoid re-computing existing MainOrients and therefore it will
% speed up the process.
% After merging/fusing the NEW MainOrients with the EXISTING ones, a
% mapping is performed to ensure complete matching in the OrientsFlagVolume. 

if ~isempty(MainOrients)
    %%% Compare MainOrientsNEW with existing MainOrients
    
    % Initial Size of MainOrients
    MOlength = size(MainOrients,1);
    
    MOstart = MOlength + 1;
    
    % Initialize the counter for NEWentries of MainOrients set.
    newEntry = 0;
    
    for mos = 1 : size(MainOrientsNEW,1)
        
        % Project the 1st bases (Ev1) of MainOrientsNEW to MainOrients
        ProjEv1s = VTF3D_ProjectEv(MainOrients(:,:,1),MainOrientsNEW(mos,:,1));
        
        % Check if ProjEv1s are ABOVE theta_threshold
        BProject1 = VTF3D_BoolProjectEv(ProjEv1s,theta_thrshld);
        
        if sum(BProject1) > 0 % If there is at least one close Ev1
            
            % Among those who are below theta_threshold, project again considering
            % the 2nd bases (Ev2)
            SelMainOrients = MainOrients;
            SelMainOrients(~BProject1,:,:) = NaN;
            
            % Project the 1st bases (Ev1) of MainOrientsNEW to MainOrients
            ProjEv2s = VTF3D_ProjectEv(SelMainOrients(:,:,2),MainOrientsNEW(mos,:,2));

            % Check if ProjEv1s are ABOVE theta_threshold
            BProject2 = VTF3D_BoolProjectEv(ProjEv2s,theta_thrshld);
            
            if sum(BProject2) > 0 % If there is at least one close to both Ev1 and Ev2
                
                % All the NEW Main Orientations that are both Ev1 and Ev2 below
                % theta_threshold can be merged and associated to existing MainOrients.
                % Those that do not match constitute NEW entries to the MainOrients
                % set.
                
                MainOrientsLBL = find(BProject2);
                
                if length(MainOrientsLBL)>1 % Multiple Solutions
                    
                    % Consider the closest
                    [~,idx] = max( mean( [sum( ProjEv1s(MainOrientsLBL,:), 2) , sum( ProjEv2s(MainOrientsLBL,:), 2 ) ] , 2 ) );
                    
                    %%% Update accordingly (pay attention to the Labels!) the
                    %%% OrientsFlagVolume without introducing newEntries.
                    OrientsFlagVolume( OrientsFlagVolume == mos ) = MainOrientsLBL(idx);
                    clear MainOrientsLBL idx;
                    
                else % Unique Solution
                    %%% Update accordingly (pay attention to the Labels!) the
                    %%% OrientsFlagVolume without introducing newEntries.
                    
                    OrientsFlagVolume( OrientsFlagVolume == mos ) = MainOrientsLBL;
                    clear MainOrientsLBL;
                end
                
                
            else
                % The considered element MUST be a NEW entry
                %%% Include NEW Main Orientations in the MainOrients set.
                newEntry = newEntry + 1;
                MainOrients( MOlength + newEntry , : , : ) = MainOrientsNEW( mos , : , : );
            
                %%% Update accordingly (pay attention to the Labels!) the
                %%% OrientsFlagVolume
                OrientsFlagVolume( OrientsFlagVolume == mos ) = MOlength + newEntry;
            end
            
        else
            % The considered element MUST be a NEW entry
            %%% Include NEW Main Orientations in the MainOrients set.
            newEntry = newEntry + 1;
            MainOrients( MOlength + newEntry , : , : ) = MainOrientsNEW( mos , : , : );
            
            %%% Update accordingly (pay attention to the Labels!) the
            %%% OrientsFlagVolume
            OrientsFlagVolume( OrientsFlagVolume == mos ) = MOlength + newEntry;
        end
        
        clear ProjEv* BProject* SelMainOrients;
        
    end
    
else
    MOstart = size(MainOrients,1) + 1;
    MainOrients = MainOrientsNEW;
end

function DFKsMO = VTF3D_configDFKsOnMainOrients(DFKs,MainOrients,DFKsMO,MOstart)

if ~isempty(MainOrients)
    %HERE DFKsMO will be a LINEAR struct [1 x (#templates * #MainOrientations)]
    
    ProgressValue = 0;
    ProgressValue_step = 100.0/size(DFKs,2);
    ProgressMessage = 'Configuring DFKs on Main Vessel Directions:';
    ProgressOldString = VTF3D_showProgressCM(ProgressMessage,'',ProgressValue,false);
    
    for ker = 1 : size(DFKs,2) % FOR ALL KERNELS -- ALWAYS
        
        GRID = DFKs(ker).GRID; % For Easier Readability
        
        for ors = MOstart : size(MainOrients,1) % FOR THE STARTING MainOrient (MOstart), 'till the END
            %%% Rotating the Kernel Data according to the MainOrients %%%
            [krot,...
             Dkrot,...
             EIGLErot,...
             EIGrot]   =  VTF3D_rotate3DKernelTensorVectorField(DFKs(ker).k,...             Kernel's Impulse Response
                                                                DFKs(ker).Dk,...            Kernel's Tubularity Probe (SLOGS)
                                                                DFKs(ker).EIGLE,...         Kernel's Tensor Vector Field (LOG-EUC)
                                                                DFKs(ker).phi,...           Original Orients of the Kernel
                                                                [MainOrients(ors,:,1);...   MainOrients 2x3 matrix
                                                                 MainOrients(ors,:,2)]   );
            
            %%% Cropping out the Padding %%%
            % Rotated Kernel Impulse Response 'krot'
            krot = krot(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            % Rotated Kernel Tubularity Probe 'Dkrot'
            Dkrot = Dkrot(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            % Rotated Log-Euclidean Tensor Vector Field 'EIGLErot'
            EIGLErot.T11 = EIGLErot.T11(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGLErot.T12 = EIGLErot.T12(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGLErot.T13 = EIGLErot.T13(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGLErot.T22 = EIGLErot.T22(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGLErot.T23 = EIGLErot.T23(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGLErot.T33 = EIGLErot.T33(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            % Rotated Linear Tensor Vector Field 'EIGrot'
            EIGrot.El1 = EIGrot.El1(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGrot.El2 = EIGrot.El2(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGrot.El3 = EIGrot.El3(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3));
            EIGrot.Ev1 = EIGrot.Ev1(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3),1:3);
            EIGrot.Ev2 = EIGrot.Ev2(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3),1:3);
            EIGrot.Ev3 = EIGrot.Ev3(GRID.pad(1)+1:end-GRID.pad(1),GRID.pad(2)+1:end-GRID.pad(2),GRID.pad(3)+1:end-GRID.pad(3),1:3);
            
            %%% Balancing the Kernel Responses %%%
            krot = krot./sum(krot(:)); % sum up to 1
            Dkrot = Dkrot - (sum(Dkrot(:))/numel(Dkrot)); % sum up to 0
            %%%
            
            % Fractional Anisotropy
            Elsrot = sort(cat(4,EIGrot.El1,EIGrot.El2,EIGrot.El3),4);
        
            % Fractional Anisotropy - Blob Elongation
            FAblobrot = Elsrot(:,:,:,1)./sqrt(Elsrot(:,:,:,2).*Elsrot(:,:,:,3));
        
            % Fractional Anisotropy - Cross-Sectional
            FAcrosrot = Elsrot(:,:,:,2)./Elsrot(:,:,:,3);
            
            %%% Assignment and Storing Data %%%
            DFKsMO(ors,ker).k = krot;
            DFKsMO(ors,ker).Dk = Dkrot;
            DFKsMO(ors,ker).EIGLE = EIGLErot;
            DFKsMO(ors,ker).EIG = EIGrot;
            DFKsMO(ors,ker).FAblob = FAblobrot;
            DFKsMO(ors,ker).FAcros = FAcrosrot;
            
            clear krot Dkrot EIGLErot EIGrot Elsrot FAblobrot FAcrosrot DFKsMOtemp;
        end
        
        % Progress Percentage
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTF3D_showProgressCM([],ProgressOldString,ProgressValue,ker == size(DFKs,2));
        
    end
    
end

function DFKsFlat = VTF3D_configFlatDFKs()

%%% Introduction of Additional FLAT Kernels to Clean and Remove Artifacts
%%% from the Filtered (Reconstructed) Scene %%% 
% NB: These Filters MUST be Used in conjunction with the NEGATIVE IMAGE
% NB: The associated Tensor Vector Field to these Additional Kernels is
% Isotropic. The TVF Averaging Process can be performed directly using the
% weights (no need of a Log-Euclidean Representation as it is identically
% null).

% BounDary Filter
bd_Dk = ones([3 3 3]); bd_Dk(2,2,2) = -26; bd_Dk = bd_Dk./numel(bd_Dk); % Normalized by the number of elements and Sum up to 0. (Pure Derivative)
bd_k = zeros([3,3,3]); bd_k(2,2,2) = 1;% Impulse Response: Dirac Impulse. Sum up to 1. (without Gain)
% BackGround Filter
bg_Dk = ones([3 3 3]); bg_Dk = bg_Dk./sum(bg_Dk(:));% Sum up to 1. (Flat Filter without Gain)
bg_k = ones([3 3 3]); bg_k = bg_k./sum(bg_k(:));% Sum up to 1. (Flat Filter without Gain)

% BounDary Filter Assignment
DFKsFlat(1,1).k = bd_k;
DFKsFlat(1,1).Dk = bd_Dk;
DFKsFlat(1,1).EIGLE.T11 = zeros([3 3 3]);
DFKsFlat(1,1).EIGLE.T12 = zeros([3 3 3]);
DFKsFlat(1,1).EIGLE.T13 = zeros([3 3 3]);
DFKsFlat(1,1).EIGLE.T22 = zeros([3 3 3]);
DFKsFlat(1,1).EIGLE.T23 = zeros([3 3 3]);
DFKsFlat(1,1).EIGLE.T33 = zeros([3 3 3]);
DFKsFlat(1,1).EIG.El1 = ones([3 3 3]);
DFKsFlat(1,1).EIG.El2 = ones([3 3 3]);
DFKsFlat(1,1).EIG.El3 = ones([3 3 3]);
DFKsFlat(1,1).EIG.Ev1 = cat(4,ones([3 3 3]),zeros([3 3 3]),zeros([3 3 3]));
DFKsFlat(1,1).EIG.Ev2 = cat(4,zeros([3 3 3]),ones([3 3 3]),zeros([3 3 3]));
DFKsFlat(1,1).EIG.Ev3 = cat(4,zeros([3 3 3]),zeros([3 3 3]),ones([3 3 3]));
DFKsFlat(1,1).FAblob = zeros([3 3 3]);
DFKsFlat(1,1).FAcros = zeros([3 3 3]);

% BackGround Filter Assignment
DFKsFlat(2,1).k = bg_k;
DFKsFlat(2,1).Dk = bg_Dk;
DFKsFlat(2,1).EIGLE.T11 = zeros([3 3 3]);
DFKsFlat(2,1).EIGLE.T12 = zeros([3 3 3]);
DFKsFlat(2,1).EIGLE.T13 = zeros([3 3 3]);
DFKsFlat(2,1).EIGLE.T22 = zeros([3 3 3]);
DFKsFlat(2,1).EIGLE.T23 = zeros([3 3 3]);
DFKsFlat(2,1).EIGLE.T33 = zeros([3 3 3]);
DFKsFlat(2,1).EIG.El1 = ones([3 3 3]);
DFKsFlat(2,1).EIG.El2 = ones([3 3 3]);
DFKsFlat(2,1).EIG.El3 = ones([3 3 3]);
DFKsFlat(2,1).EIG.Ev1 = cat(4,ones([3 3 3]),zeros([3 3 3]),zeros([3 3 3]));
DFKsFlat(2,1).EIG.Ev2 = cat(4,zeros([3 3 3]),ones([3 3 3]),zeros([3 3 3]));
DFKsFlat(2,1).EIG.Ev3 = cat(4,zeros([3 3 3]),zeros([3 3 3]),ones([3 3 3]));
DFKsFlat(2,1).FAblob = zeros([3 3 3]);
DFKsFlat(2,1).FAcros = zeros([3 3 3]);

function [Krot,DKrot,EIGLErot,EIGrot] = VTF3D_rotate3DKernelTensorVectorField(K,DK,EIGLE,MainOrients,RotatedOrients)

% Determining the Rotation Matrix
R = VTF3D_RotMatrix4Vect(MainOrients(1,:),MainOrients(2,:),RotatedOrients(1,:),RotatedOrients(2,:));

% Rotating the Filtering Kernel
Krot = VTF3D_Rot3DwithR(K,R,NaN);
DKrot = VTF3D_Rot3DwithR(DK,R,NaN);

% Rotate Each Component of the LogEuclidean Tensor
T11temp = VTF3D_Rot3DwithR(EIGLE.T11,R,0);
T12temp = VTF3D_Rot3DwithR(EIGLE.T12,R,0);
T13temp = VTF3D_Rot3DwithR(EIGLE.T13,R,0);
T22temp = VTF3D_Rot3DwithR(EIGLE.T22,R,0);
T23temp = VTF3D_Rot3DwithR(EIGLE.T23,R,0);
T33temp = VTF3D_Rot3DwithR(EIGLE.T33,R,0);

%%% Converting to Linear Domain %%%
[El1rot,El2rot,El3rot,Ev1,Ev2,Ev3] = VTF3D_convert6LogEuclidComponentsTo3DTensor(T11temp,T12temp,T13temp,T22temp,T23temp,T33temp,true(size(T11temp)));

% Counter-Rotation for the EigenVectors
R = R';

Ev1rot = cat(4, R(1,1)*Ev1(:,:,:,1) + R(1,2)*Ev1(:,:,:,2) + R(1,3)*Ev1(:,:,:,3) ,...
                R(2,1)*Ev1(:,:,:,1) + R(2,2)*Ev1(:,:,:,2) + R(2,3)*Ev1(:,:,:,3) ,...
                R(3,1)*Ev1(:,:,:,1) + R(3,2)*Ev1(:,:,:,2) + R(3,3)*Ev1(:,:,:,3) );

Ev2rot = cat(4, R(1,1)*Ev2(:,:,:,1) + R(1,2)*Ev2(:,:,:,2) + R(1,3)*Ev2(:,:,:,3) ,...
                R(2,1)*Ev2(:,:,:,1) + R(2,2)*Ev2(:,:,:,2) + R(2,3)*Ev2(:,:,:,3) ,...
                R(3,1)*Ev2(:,:,:,1) + R(3,2)*Ev2(:,:,:,2) + R(3,3)*Ev2(:,:,:,3) );

Ev3rot = cat(4, R(1,1)*Ev3(:,:,:,1) + R(1,2)*Ev3(:,:,:,2) + R(1,3)*Ev3(:,:,:,3) ,...
                R(2,1)*Ev3(:,:,:,1) + R(2,2)*Ev3(:,:,:,2) + R(2,3)*Ev3(:,:,:,3) ,...
                R(3,1)*Ev3(:,:,:,1) + R(3,2)*Ev3(:,:,:,2) + R(3,3)*Ev3(:,:,:,3) );

% Assignment
EIGrot.El1 = El1rot;
EIGrot.El2 = El2rot;
EIGrot.El3 = El3rot;
EIGrot.Ev1 = Ev1rot;
EIGrot.Ev2 = Ev2rot;
EIGrot.Ev3 = Ev3rot;

% Re-converting into the LE-domain
[T11rot,T12rot,T13rot,T22rot,T23rot,T33rot] = VTF3D_convert3DTensorTo6LogEuclidComponents(El1rot,El2rot,El3rot,Ev1rot,Ev2rot,Ev3rot);

clear El1* El2* El3* Ev1* Ev2* Ev3*

% Assignment
EIGLErot.T11 = T11rot;
EIGLErot.T12 = T12rot;
EIGLErot.T13 = T13rot;
EIGLErot.T22 = T22rot;
EIGLErot.T23 = T23rot;
EIGLErot.T33 = T33rot;

clear T*
    
function [T11LE,T12LE,T13LE,T22LE,T23LE,T33LE] = VTF3D_convert3DTensorTo6LogEuclidComponents(El1,El2,El3,Ev1,Ev2,Ev3)
w = [1,sqrt(2),sqrt(2),1,sqrt(2),1];
    
T11LE =  w(1) * (...
                        Ev1(:,:,:,1) .* log( El1 ) .* conj(Ev1(:,:,:,1)) + ... e1x*LOG(l1)*conj(e1x) +  
                        Ev2(:,:,:,1) .* log( El2 ) .* conj(Ev2(:,:,:,1)) + ... e2x*LOG(l2)*conj(e2x) +
                        Ev3(:,:,:,1) .* log( El3 ) .* conj(Ev3(:,:,:,1))   ... e3x*LOG(l3)*conj(e3x)
                        ); % Equivalent to the Tensor Component (1,1)
    
T12LE =  w(2) * (... 
                        Ev1(:,:,:,2) .* log( El1 ) .* conj(Ev1(:,:,:,1)) + ... e1y*LOG(l1)*conj(e1x) + 
                        Ev2(:,:,:,2) .* log( El2 ) .* conj(Ev2(:,:,:,1)) + ... e2y*LOG(l2)*conj(e2x) +
                        Ev3(:,:,:,2) .* log( El3 ) .* conj(Ev3(:,:,:,1))   ... e3y*LOG(l3)*conj(e3x)
                        ); % Equivalent to the Tensor Component (1,2)
    
T13LE =  w(3) * (...
                        Ev1(:,:,:,3) .* log( El1 ) .* conj(Ev1(:,:,:,1)) + ... e1z*LOG(l1)*conj(e1x) + 
                        Ev2(:,:,:,3) .* log( El2 ) .* conj(Ev2(:,:,:,1)) + ... e2z*LOG(l2)*conj(e2x) + 
                        Ev3(:,:,:,3) .* log( El3 ) .* conj(Ev3(:,:,:,1))   ... e3z*LOG(l3)*conj(e3x)
                        ); % Equivalent to the Tensor Component (1,3)
    
T22LE =  w(4) * (...
                        Ev1(:,:,:,2) .* log( El1 ) .* conj(Ev1(:,:,:,2)) + ... e1y*LOG(l1)*conj(e1y) + 
                        Ev2(:,:,:,2) .* log( El2 ) .* conj(Ev2(:,:,:,2)) + ... e2y*LOG(l2)*conj(e2y) + 
                        Ev3(:,:,:,2) .* log( El3 ) .* conj(Ev3(:,:,:,2))   ... e3y*LOG(l3)*conj(e3y)
                        ); % Equivalent to the Tensor Component (2,2)
                   
T23LE =  w(5) * (...
                        Ev1(:,:,:,3) .* log( El1 ) .* conj(Ev1(:,:,:,2)) + ... e1z*LOG(l1)*conj(e1y) + 
                        Ev2(:,:,:,3) .* log( El2 ) .* conj(Ev2(:,:,:,2)) + ... e2z*LOG(l2)*conj(e2y) + 
                        Ev3(:,:,:,3) .* log( El3 ) .* conj(Ev3(:,:,:,2))   ... e3z*LOG(l3)*conj(e3y)
                        ); % Equivalent to the Tensor Component (2,3)
    
T33LE =  w(6) * (...
                        Ev1(:,:,:,3) .* log( El1 ) .* conj(Ev1(:,:,:,3)) + ... e1z*LOG(l1)*conj(e1z) + 
                        Ev2(:,:,:,3) .* log( El2 ) .* conj(Ev2(:,:,:,3)) + ... e2z*LOG(l2)*conj(e2z) + 
                        Ev3(:,:,:,3) .* log( El3 ) .* conj(Ev3(:,:,:,3))   ... e3z*LOG(l3)*conj(e3z)
                        ); % Equivalent to the Tensor Component (3,3)

function [El1,El2,El3,Ev1,Ev2,Ev3] = VTF3D_convert6LogEuclidComponentsTo3DTensor(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE,Msk)

w = 1./[1,sqrt(2),sqrt(2),1,sqrt(2),1];
%%%
% Dummy T fields:
T.Ix = zeros(size(T11LE));T.Iy = zeros(size(T11LE));T.Iz = zeros(size(T11LE));
% Real T fields:
T.Ixx = w(1) * T11LE;
T.Ixy = w(2) * T12LE;
T.Ixz = w(3) * T13LE;
T.Iyy = w(4) * T22LE;
T.Iyz = w(5) * T23LE;
T.Izz = w(6) * T33LE;
% Dummy T fields:
T.Inorm = zeros(size(T11LE));T.IxN = zeros(size(T11LE));T.IyN = zeros(size(T11LE));T.IzN = zeros(size(T11LE));
%%%

%                                                                               
[El1,El2,El3,Ev1,Ev2,Ev3, ~ ] = VTF3D_get3DImageEIGs( T , Msk , true , true); % Sort ENABLED and Absolute value ENABLED!
%

El1 = exp(El1); El1(~Msk) = NaN;
El2 = exp(El2); El2(~Msk) = NaN;
El3 = exp(El3); El3(~Msk) = NaN;

Ev1_1 = Ev1(:,:,:,1); Ev1_1(~Msk) = NaN;
Ev1_2 = Ev1(:,:,:,2); Ev1_2(~Msk) = NaN;
Ev1_3 = Ev1(:,:,:,3); Ev1_3(~Msk) = NaN;
Ev1 = cat(4,Ev1_1,Ev1_2,Ev1_3);

Ev2_1 = Ev2(:,:,:,1); Ev2_1(~Msk) = NaN;
Ev2_2 = Ev2(:,:,:,2); Ev2_2(~Msk) = NaN;
Ev2_3 = Ev2(:,:,:,3); Ev2_3(~Msk) = NaN;
Ev2 = cat(4,Ev2_1,Ev2_2,Ev2_3);

Ev3_1 = Ev3(:,:,:,1); Ev3_1(~Msk) = NaN;
Ev3_2 = Ev3(:,:,:,2); Ev3_2(~Msk) = NaN;
Ev3_3 = Ev3(:,:,:,3); Ev3_3(~Msk) = NaN;
Ev3 = cat(4,Ev3_1,Ev3_2,Ev3_3);

function [Vesselness,Boundaries,Background,TF,TFLE] = VTF3D_computeAllMaps(Idwn,DFKsMOsel,DFKsFlat,EuclidTFflag,SkipTFflag)

% Determining the Negative Image 'nIdwn' (Necessary for the FLAT DFKs)
nIdwn = -Idwn;
nIdwn = (nIdwn - min(nIdwn(:)))./(max(nIdwn(:)) - min(nIdwn(:))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the Hann Window for the Distance Weights (DW)
% The 0.5 weight is designed on PURPOSE so that there is perfect
% reconstruction by considering adjacent positions in a kernel of size 
% [5 5 5], which is EQUAL to the size of the DFKsMO!!!
TensorKernelSize = size(DFKsMOsel(1,1).k);
TensorImgPadSize = TensorKernelSize - 1;
DW = reshape(kron(kron(0.5*hann(TensorKernelSize(1)),0.5*hann(TensorKernelSize(2))'),0.5*hann(TensorKernelSize(3))),TensorKernelSize); % N.B. The size of H is EQUAL to the Size of any Rotated Kernel!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filtering with SLoGS filterbank -- Scalar Response
ProgressMessage = 'Image Filtering with SLoGS:';
[FR,FRidxDFK] = VTF3D_filtScalarMap(Idwn,DFKsMOsel,0,[0,1],ProgressMessage);
[FRbd,  ~   ] = VTF3D_filtScalarMap(nIdwn,DFKsFlat(1,1), -inf,[-1,1],[]); % Boundaries
[FRbg,  ~   ] = VTF3D_filtScalarMap(nIdwn,DFKsFlat(2,1),  0  , [0,1],[]); % Background

FRbd = (FRbd - min(FRbd(:)))/(max(FRbd(:))-min(FRbd(:)));

% Filtering with SLoGS filterbank -- Synthesis of the Tensor Field
ProgressMessage1 = 'Synthesising Anisotropic Tensor Field:';

if ~SkipTFflag
    
    TFStruct = VTF3D_initializeTF(size(Idwn),TensorImgPadSize);
    TFStruct.DW = DW;
    
    % Configuring Tensor Synthesis
    abs_maxitr = 64;
    maxitr = max(size(Idwn));
    maxitr(maxitr > abs_maxitr) = abs_maxitr;
    maxitr = maxitr^3;
    
    totcycles = ceil(numel(Idwn)/maxitr);
    vxl_ini = 1;
    
    ProgressValue = 0;
    ProgressValue_step = 100.0/totcycles;
    ProgressOldString = VTF3D_showProgressCM(ProgressMessage1,'',ProgressValue,false);
    
    for cyl = 1 : totcycles
        [TFSynthStruct,vxl_fin] = VTF3D_configTFSynthesis_mex(size(Idwn),maxitr,vxl_ini); % MEX
        TFStruct = VTF3D_synthTF_mex(TFStruct,TFSynthStruct,FR,FRidxDFK,DFKsMOsel,DFKsFlat,vxl_ini-1,0); % MEX
        vxl_ini = vxl_fin;
        % Progress Percentage
        ProgressValue = ProgressValue + ProgressValue_step;
        ProgressOldString = VTF3D_showProgressCM([],ProgressOldString,ProgressValue, false );
    end
    
    TFStruct = VTF3D_removePadTFStruct(TFStruct,TensorImgPadSize);
    
    %%% PROCEDURE to RECOVER the Structure Tensors Field (STF) ARTIFACTS after
    %%% MULTIPLE Integration and Averaging (as in CVM3D_synth3DOLABlk_mex).
    %%% STF: filter responses Weighted Integration (WI) with Boundary and
    %%% Background Compensation (BBC) + Anisotropic Ratio Recovery (ARR) +
    %%% Boundary Restriction Correction (BRC) to avoid flaking effects on STF.
    %%%%%%%%%%%%%%%%%%%%%%%%%%% STF_WI_BBC_ARR_BRC %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [TF,TFLE] = VTF3D_recoverTensorField(FR,TFStruct,EuclidTFflag);
    
else
    disp([ProgressMessage1,' - SKIP Tensor Field Synthesis -']);
    TF = [];
    TFLE = [];
end

Vesselness = FR;
Boundaries = FRbd;
Background = FRbg;

disp(' -- DONE');

function TFStruct = VTF3D_initializeTF(ImgSize,TensorImgPadSize)

Size = ImgSize + (2*TensorImgPadSize);

% Initializing the Structure of Tensor Field (STF)
TFStruct.T11 = zeros( Size );
TFStruct.T12 = zeros( Size );
TFStruct.T13 = zeros( Size );
TFStruct.T22 = zeros( Size );
TFStruct.T23 = zeros( Size );
TFStruct.T33 = zeros( Size );

% Initializing the Weights Accumulator
TFStruct.WA = zeros( Size );

TFStruct.ARelong = zeros( Size );
TFStruct.ARcross = zeros( Size );

function TFStruct = VTF3D_removePadTFStruct(TFStruct,TensorImgPadSize)
% Removing TensorImgPadSize...      

TFStruct.T11 = TFStruct.T11(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.T12 = TFStruct.T12(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.T13 = TFStruct.T13(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.T22 = TFStruct.T22(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.T23 = TFStruct.T23(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.T33 = TFStruct.T33(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.WA  =  TFStruct.WA(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                            TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                            TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
                        
TFStruct.ARelong  =  TFStruct.ARelong(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                                      TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                                      TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);
TFStruct.ARcross  =  TFStruct.ARcross(TensorImgPadSize(1)/2 + 1 : end - TensorImgPadSize(1)*3/2,...
                                      TensorImgPadSize(2)/2 + 1 : end - TensorImgPadSize(2)*3/2,...
                                      TensorImgPadSize(3)/2 + 1 : end - TensorImgPadSize(3)*3/2);

function [TF,TFLE] = VTF3D_recoverTensorField(FR,TFStruct,EuclidTFflag)

% Referencing...
T11 = TFStruct.T11;
T12 = TFStruct.T12;
T13 = TFStruct.T13;
T22 = TFStruct.T22;
T23 = TFStruct.T23;
T33 = TFStruct.T33;
WA = TFStruct.WA;

% Removing Padding from the TVF Weighted Integration
T11r = T11 ./ WA;
T12r = T12 ./ WA;
T13r = T13 ./ WA;
T22r = T22 ./ WA;
T23r = T23 ./ WA;
T33r = T33 ./ WA;

% Correct for Degenerate Cases in the Log-Euclidean Domain
% Correcting for Degenerate Cases (which will produce NaN in the Linear Domain of the TVF)
T11r(abs(T11r)<1e-6 | ~isfinite(T11r)) = 1e-6;
T12r(abs(T12r)<1e-6 | ~isfinite(T12r)) = 1e-6; 
T13r(abs(T13r)<1e-6 | ~isfinite(T13r)) = 1e-6; 
T22r(abs(T22r)<1e-6 | ~isfinite(T22r)) = 1e-6;
T23r(abs(T23r)<1e-6 | ~isfinite(T23r)) = 1e-6; 
T33r(abs(T33r)<1e-6 | ~isfinite(T33r)) = 1e-6;

if EuclidTFflag % if the Euclidean Tensor Filed must be retrieved (FullDump Enabled)
    LogicalMap = (FR >= quantile(FR(FR>0),0.9));
    
    PowerFactor = 3;
    
    % Converting the TVF from the Log-Euclidean Domain to the Linear one.
    [El1r,El2r,El3r,Ev1r,Ev2r,Ev3r] = VTF3D_convert6LogEuclidComponentsTo3DTensor(T11r,T12r,T13r,T22r,T23r,T33r, LogicalMap );
    
    %       Vessel Component
    El1r = (1./El1r).^(PowerFactor);% Smallest!
    El2r = (1./El2r).^(PowerFactor);
    El3r = (1./El3r).^(PowerFactor);
    
    % Final Normalization s.t. the Product of the Eigenvalue Maps is equal to 1
    ElrProd = (El1r.*El2r.*El3r).^(1/3);
    
    % Final Assignment s.t. the determinant is equal to 1.
    El1r = El1r ./ ElrProd;
    El2r = El2r ./ ElrProd; % This Should be the SMALLEST
    El3r = El3r ./ ElrProd;
    
    % Correct for Degenerate Cases in the Euclidean Domain
    [El1r,El2r,El3r,Ev1r,Ev2r,Ev3r] = VTF3D_correctDegenerateTensors(El1r,El2r,El3r,Ev1r,Ev2r,Ev3r);
    
    % Assignment
    TF.El1 = El1r;
    TF.El2 = El2r;
    TF.El3 = El3r;
    TF.Ev1 = Ev1r;
    TF.Ev2 = Ev2r;
    TF.Ev3 = Ev3r;
else
    TF = [];
end

TFLE.T11 = T11r;
TFLE.T12 = T12r;
TFLE.T13 = T13r;
TFLE.T22 = T22r;
TFLE.T23 = T23r;
TFLE.T33 = T33r;

function [El1,El2,El3,Ev1,Ev2,Ev3] = VTF3D_correctDegenerateTensors(El1,El2,El3,Ev1,Ev2,Ev3)

%%% Correcting for possible Degenerate Cases %%%
El1(~isfinite(El1) | El1 == 0) = 1;
El2(~isfinite(El2) | El2 == 0) = 1;
El3(~isfinite(El3) | El3 == 0) = 1;

Ev1_1 = Ev1(:,:,:,1); Ev1_1(~isfinite(Ev1_1)) = 1;
Ev1_2 = Ev1(:,:,:,2); Ev1_2(~isfinite(Ev1_2)) = 0;
Ev1_3 = Ev1(:,:,:,3); Ev1_3(~isfinite(Ev1_3)) = 0;
Ev1 = cat(4,Ev1_1,Ev1_2,Ev1_3);

Ev2_1 = Ev2(:,:,:,1); Ev2_1(~isfinite(Ev2_1)) = 0;
Ev2_2 = Ev2(:,:,:,2); Ev2_2(~isfinite(Ev2_2)) = 1;
Ev2_3 = Ev2(:,:,:,3); Ev2_3(~isfinite(Ev2_3)) = 0;
Ev2 = cat(4,Ev2_1,Ev2_2,Ev2_3);

Ev3_1 = Ev3(:,:,:,1); Ev3_1(~isfinite(Ev3_1)) = 0;
Ev3_2 = Ev3(:,:,:,2); Ev3_2(~isfinite(Ev3_2)) = 0;
Ev3_3 = Ev3(:,:,:,3); Ev3_3(~isfinite(Ev3_3)) = 1;
Ev3 = cat(4,Ev3_1,Ev3_2,Ev3_3);

function DFKsMOsel = retrieveSelectedDFKsMO(DFKsMO,OrientsFlagVolume,DFKs)

OrientsSel = unique(OrientsFlagVolume(OrientsFlagVolume > 0));

if ~isempty(OrientsSel)
    DFKsMOsel = DFKsMO(OrientsSel,:);
    DFKsMOsel = reshape(DFKsMOsel,[numel(DFKsMOsel),1]);
else
    disp('No Main Vessel Orientations were found -- Considering Icosphere');
    IcoSubLvl = 1;
    IcoPrincDirects = VTF3D_PrincipalDirectionsFromIcosphere(IcoSubLvl,false); % e.g.  3-Levels of Subdivision
    MainOrients = zeros([length(IcoPrincDirects),3,2]);
    
    for ipd = 1 : length(IcoPrincDirects)
       MainOrients(ipd,:,1) = IcoPrincDirects(ipd).M(1,:);
       MainOrients(ipd,:,2) = IcoPrincDirects(ipd).M(2,:);
    end
    
    DFKsMOsel = VTF3D_configDFKsOnMainOrients(DFKs,MainOrients,[],1);
    DFKsMOsel = reshape(DFKsMOsel,[numel(DFKsMOsel),1]);
end

function TFs_single = VTF3D_makeTFsingle(TFs)

if ~isempty(TFs)
    TFs_single.El1 = single(TFs.El1);
    TFs_single.El2 = single(TFs.El2);
    TFs_single.El3 = single(TFs.El3);
    TFs_single.Ev1 = single(TFs.Ev1);
    TFs_single.Ev2 = single(TFs.Ev2);
    TFs_single.Ev3 = single(TFs.Ev3);
else
    TFs_single = [];
end

function TFsLE_single = VTF3D_makeTFsLEsingle(TFsLE)

if ~isempty(TFsLE)
    TFsLE_single.T11 = single(TFsLE.T11);
    TFsLE_single.T12 = single(TFsLE.T12);
    TFsLE_single.T13 = single(TFsLE.T13);
    TFsLE_single.T22 = single(TFsLE.T22);
    TFsLE_single.T23 = single(TFsLE.T23);
    TFsLE_single.T33 = single(TFsLE.T33);
else
    TFsLE_single = [];
end

function [MainOrients,DFKsMO,newEntry] = VTF3D_mergeMainOrientsAndDFKsMO(MainOrients,MainOrientsNEW,DFKsMO,DFKsMOnew,theta_thrshld)
% Use this function to compare and merge MainOrientsNEW to existing
% MainOrients. Accordingly, DFKsMOnew will be merged to DFKsMO.
% This will avoid re-computing existing MainOrients (and DFKsMO) therefore,
% the more images processed, the less the un
% it will speed up the process.
% After merging/fusing the NEW MainOrients with the EXISTING ones, a
% mapping is performed to ensure complete matching in the OrientsFlagVolume. 

if ~isempty(MainOrients)
    %%% Compare MainOrientsNEW with existing MainOrients
    
    % Initial Size of MainOrients
    MOlength = size(MainOrients,1);
    
    % Initialize the counter for NEWentries of MainOrients set.
    newEntry = 0;
    
    for mos = 1 : size(MainOrientsNEW,1)
        
        % Project the 1st bases (Ev1) of MainOrientsNEW to MainOrients
        ProjEv1s = VTF3D_ProjectEv(MainOrients(:,:,1),MainOrientsNEW(mos,:,1));
        
        % Check if ProjEv1s are ABOVE theta_threshold
        BProject1 = VTF3D_BoolProjectEv(ProjEv1s,theta_thrshld);
        
        if sum(BProject1) > 0 % If there is at least one close Ev1
            
            % Among those who are below theta_threshold, project again considering
            % the 2nd bases (Ev2)
            SelMainOrients = MainOrients;
            SelMainOrients(~BProject1,:,:) = NaN;
            
            % Project the 1st bases (Ev1) of MainOrientsNEW to MainOrients
            ProjEv2s = VTF3D_ProjectEv(SelMainOrients(:,:,2),MainOrientsNEW(mos,:,2));

            % Check if ProjEv1s are ABOVE theta_threshold
            BProject2 = VTF3D_BoolProjectEv(ProjEv2s,theta_thrshld);
            
            if sum(BProject2) == 0 % There is NO matching or close MainOrients (Ev1 and Ev2) item relative to the MainOrientsNEW element
                
                % The considered element in MainOrientsNEW MUST be a NEW entry
                %%% Include NEW Main Orientations in the MainOrients set.
                newEntry = newEntry + 1;
                MainOrients( MOlength + newEntry , : , : ) = MainOrientsNEW( mos , : , : );
            
                %%% Update accordingly DFKsMO with the element in DFKsMOnew
                DFKsMO( MOlength + newEntry , : ) = DFKsMOnew( mos , : );
            end
            
        else% There is NO matching or close MainOrients (Ev1) item relative to the MainOrientsNEW element
            
            % The considered element MUST be a NEW entry
            %%% Include NEW Main Orientations in the MainOrients set.
            newEntry = newEntry + 1;
            MainOrients( MOlength + newEntry , : , : ) = MainOrientsNEW( mos , : , : );
            
            %%% Update accordingly DFKsMO with the element in DFKsMOnew
            DFKsMO( MOlength + newEntry , : ) = DFKsMOnew( mos , : );
        end
        
        clear ProjEv* BProject* SelMainOrients;
        
    end
    
else
    MainOrients = MainOrientsNEW;
    DFKsMO = DFKsMOnew;
    newEntry = size(MainOrients,1);
end

function [mu,sigma] = VTF3D_estimate3DImageGaussianNoise(I)

rndsamples = round(1e-5 * numel(I));

kernelsize = 9;

halfkernelsize = floor(kernelsize/2);

valsrange = [0 + halfkernelsize + 1 , size(I,1) - (halfkernelsize + 1);
             0 + halfkernelsize + 1 , size(I,2) - (halfkernelsize + 1);
             0 + halfkernelsize + 1 , size(I,3) - (halfkernelsize + 1)];

mus = NaN(1,rndsamples);
sigmas = NaN(1,rndsamples);
       

%tic;
for jj = 1 : rndsamples
   
    rndIDX = [randi([valsrange(1,1),valsrange(1,2)]) , randi([valsrange(2,1),valsrange(2,2)]) , randi([valsrange(3,1),valsrange(3,2)]) ];
    
    block = I( rndIDX(1) - halfkernelsize : rndIDX(1) + halfkernelsize , ...
               rndIDX(2) - halfkernelsize : rndIDX(2) + halfkernelsize , ... 
               rndIDX(3) - halfkernelsize : rndIDX(3) + halfkernelsize );
    
    mus(jj) = nanmean(block(:));
    sigmas(jj) = nanstd(block(:));
    
end

mus = mus(~isnan(mus));
sigmas = sigmas(~isnan(sigmas));

if ~isempty(mus)
    mu = mean(mus);
    sigma = mean(sigmas);
else
    mu = [];
    sigma = [];
end

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

function [Ivess,Iskull] = VTF3D_skullStripCTA(I,maxval)

vals = maxval + abs(min(I(:)));
I = I + abs(min(I(:)));

valsN = vals./max(I(:));
Inorm = I./max(I(:));

% Estimate Background
[ff,bins] = hist(Inorm(Inorm > 0),100);
binspos = find(abs(diff(ff/sum(ff)))<1e-3,1);
binspos(binspos<1) = 1;

bkg_mu = 1.05 * bins(binspos); % +5%
Ibkg = (-1./(1+exp(-100*(Inorm-bkg_mu))))+1;

mu = 1.05 * valsN; % +5%
Imult = (-1./(1+exp(-100*(Inorm-mu))))+1;

% Skull
Iskull = 1-Imult;

ball = false([5,5,5]);
ball(3,3,3) = true;
ball = bwdist(ball);
ball = ball <= 1.5;

disk = false([13,13,3]);
disk(7,7,1) = true;
disk = bwdist(disk);
disk = disk <= 7;

% Dilated and Smoothed
Iskull = imdilate(Iskull,ball);
Ibkg = imdilate(Ibkg,disk);

newVal4Bone = mean(Inorm(Iskull>=0.49 & Iskull<=0.51));
%newVal4Bkg = mean(Inorm(Ibkg>=0.49 & Ibkg<=0.51));

Ivess = (1-Ibkg) .* ( Inorm.*(1-Iskull) + ( ( newVal4Bone + 0.001*randn(size(Inorm)) ).*(Iskull) ) ) + Ibkg.*( newVal4Bone + 0.001*randn(size(Inorm)) );
Ivess = (Ivess-min(Ivess(:)))./(max(Ivess(:))-min(Ivess(:)));

Iskull = Iskull>=0.5 | Ibkg >=0.5;

function Vr = VTF3D_ordfilt3D(V0,ord,padoption)
% ordfilt3D:    Perform 3-D order-statistic filtering on 26 neighbors
%
%   [Vr] = ordfilt3D(V0,ord,padoption)
%          use 26 neighbors
%       ord = 14 <=> median filtering
%       ord = 1 <=> min
%       ord = [1 27] <=> [min max]
%       padoption: same as in padarray
%
% Olivier Salvado, Case Western Reserve University, 16Aug04

if ~exist('padoption','var')
    padoption = 'replicate';
end

% special care for uint8
if isa(V0,'uint8')
    V = uint8(padarray(V0,[1 1 1],padoption));
    S = size(V);
    Vn = uint8(zeros(S(1),S(2),S(3),26));  % all the neighbor
else
    V = single(padarray(V0,[1 1 1],padoption));
    S = size(V);
    Vn = single(zeros(S(1),S(2),S(3),26));  % all the neighbor
end

% build the neighboord
Vn(:,:,:,1) = V;
i = 1:S(1); ip1 = [i(2:end) i(end)]; im1 = [i(1) i(1:end-1)];
j = 1:S(2); jp1 = [j(2:end) j(end)]; jm1 = [j(1) j(1:end-1)];
k = 1:S(3); kp1 = [k(2:end) k(end)]; km1 = [k(1) k(1:end-1)];

% left
Vn(:,:,:,2)     = V(im1    ,jm1    ,km1);
Vn(:,:,:,3)     = V(im1    ,j      ,km1);
Vn(:,:,:,4)     = V(im1    ,jp1    ,km1);

Vn(:,:,:,5)     = V(im1    ,jm1    ,k);
Vn(:,:,:,6)     = V(im1    ,j      ,k);
Vn(:,:,:,7)     = V(im1    ,jp1    ,k);

Vn(:,:,:,8)     = V(im1    ,jm1    ,kp1);
Vn(:,:,:,9)     = V(im1    ,j      ,kp1);
Vn(:,:,:,10)    = V(im1    ,jp1    ,kp1);

% right
Vn(:,:,:,11)    = V(ip1    ,jm1    ,km1);
Vn(:,:,:,12)    = V(ip1    ,j      ,km1);
Vn(:,:,:,13)    = V(ip1    ,jp1    ,km1);

Vn(:,:,:,14)    = V(ip1    ,jm1    ,k);
Vn(:,:,:,15)    = V(ip1    ,j      ,k);
Vn(:,:,:,16)    = V(ip1    ,jp1    ,k);

Vn(:,:,:,17)    = V(ip1    ,jm1    ,kp1);
Vn(:,:,:,18)    = V(ip1    ,j      ,kp1);
Vn(:,:,:,19)    = V(ip1    ,jp1    ,kp1);

% top
Vn(:,:,:,20)    = V(i       ,jm1    ,kp1);
Vn(:,:,:,21)    = V(i       ,j      ,kp1);
Vn(:,:,:,22)    = V(i       ,jp1    ,kp1);

% bottom
Vn(:,:,:,23)    = V(i       ,jm1    ,km1);
Vn(:,:,:,24)    = V(i       ,j      ,km1);
Vn(:,:,:,25)    = V(i       ,jp1    ,km1);

% front
Vn(:,:,:,26)    = V(i       ,jp1    ,k);

% back
Vn(:,:,:,27)    = V(i       ,jm1    ,k);

% perform the processing
Vn = sort(Vn,4);
Vr = Vn(:,:,:,ord);

% remove padding on the 3 first dimensions
Vr = double(Vr(2:end-1,2:end-1,2:end-1,:));

function includeVTrailsLibs

addpath(strcat(VTRootDir,'libs/'));
addpath(strcat(VTRootDir,'libs/0_NIfTI_IO/'));
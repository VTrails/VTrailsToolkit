function MaxMAP = VTF3D_IntegrateMaxMultiScaleResponses(varargin)%
% Use this function to Integrate the Multi-Scales Responses after SLoGS
% filtering.
%
% INPUTS: [OPTIONAL]
%
% - ['PathFolder']: Path of the folder containing the *.mat datasets with
%                   the multi-resolution filter responses from the function
%                   'VTrailsFilter3D_MAIN.m'
% - ['ScalesRange']: 2-Elements Array with the minScale and maxScale for
%                    the integration.
% - ['Alpha']: Modulation coefficient for boosting Boundaries and
%              Background components - Default: 1.0
%              [See refs. for more insights].
% - ['PowerFactor']: Power for integrating the Connected Vesselness Map
%                    over scales - Default: 1.0
% - ['Msk']: Validity Mask for Tensor Field integration - Default: []
%            (empty), i.e. all image domain will be considered.
% - ['ProcessTensor']: Boolean for integrating also the Tensor Field -
%                      Default: true
% - ['RegulariseTFflag']: Boolean for the Smooth Regularisation of the
%                         Synthetic Tensor Field. Default: True.
% - ['QuantileTensorTHR']:
% - [Up2Orig]: Boolean flag to resample up to the original size of the
%                image. Default: False.
%
% OUTPUTS:
%
% - MaxMAP: structure containing the following Fields,
%
%  * IMG: Header of the input nifty Images.
%  * CVM: Integral Scalar Connected Vesselness Map obtained with SLoGS
%  * CVMsmth: As above, but slightly smoothed for local maxima.
%  * BDM: Integral Scalar Vessel Boundaries Map obtained with \deltaSLoGS
%  * BGM: Integral Scalar Vessel Background Map obtained with \nuSLoGS
%  * TFLE: Integral Vascular Tensors Field synthesized with SLoGS in
%          the Log-Euclidean Domain 
%  * TF: Integral Vascular Tensors Field synthesized with SLoGS in the
%        Euclidean Domain 
%  * GRID: Structure with Image Grid and Scales' infos for Resampling.
%  * ACC: Internal Variable for Across-Scales Integration Regularisation
%  * US: Scalar Volume of Un-Organised Seeds TO BE Aligned to the
%        vessels with 'VTF3D_AlignVesselSeeds3D.m'
%
% Example:
%
%   MaxMAP = VTF3D_integrateMaxMultiScaleResponses('PathFolder','<Path/Of/VTrailsFilter/JobDump>','ScalesRange',[0.1,1.0]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, Sept-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('***** VTrails: Integrating Multi-Resolution Filter Responses Across Scales *****');

if nargin < 1
    % Default Parameters for Scale-Integration + Interactive PathFolder
    % Selection
    MaxScale = 1.0;
    MinScale = 0.0;
    alpha = 1;
    Up2Orig = false;
    PowerFactor = 1;
    QuantileTensorTHR = 0.9;
    Msk = [];
    PathFolder = uigetdir(pwd,'Select the Folder with the Multi-Scale SLoGS Responses');
    ProcessTensor = questdlg('Integrate also Tensor Field?','Tensor Field Integration Dialog','Yes','No','Yes');
else
    
    [PathFolder,ScalesRange,alpha,Up2Orig,PowerFactor,QuantileTensorTHR,Msk,ProcessTensor] = configureInputs(varargin);
    
     MinScale = ScalesRange(1);
     MaxScale = ScalesRange(2);
    
end

MaxMAP = [];

if PathFolder
    list = dir(strcat(PathFolder,'/TempDataset*.mat'));
    if ~isempty(list)
        % Loading Data
        disp(strcat('Loading Datasets Folder Path: ',PathFolder));
        
        availableScales = [];
        
        for jj = 1 : length(list)
            availableScales = cat(1,availableScales,str2num(list(jj).name(end-7:end-4)));
        end
        
        scales2load = availableScales >= MinScale & availableScales <= MaxScale;
        
        entry = 0;
        
        processTensorFieldCHK = [];
        
        for jj = 1 : length(list)
            if scales2load(jj)
                load(strcat(PathFolder,'/',list(jj).name));
                entry = entry + 1;
                if entry == 1
                    MAPcat = MAP;
                    processTensorFieldCHK = cat(1,processTensorFieldCHK,~isempty(MAP.TFLE));
                else
                    MAPcat = cat(2,MAPcat,MAP);
                    processTensorFieldCHK = cat(1,processTensorFieldCHK,~isempty(MAP.TFLE));
                end
            end
        end
        
        clear MAP;

        availableScales = availableScales(availableScales >= MinScale & availableScales <= MaxScale);
        availableMinScale = min(availableScales);
        availableMaxScale = max(availableScales);
        
        if ~logical(prod(processTensorFieldCHK)) && strcmpi(ProcessTensor,'Yes')
            warning('Inconsistent Tensor Field: Integration will be skipped for TF!');
            ProcessTensor = 'No';
        end
        
        % Integrating
        disp(strcat('Integrating Available Scales Range: [',num2str(availableMinScale,'%.1f'),' -> ',num2str(availableMaxScale,'%.1f'),'] ...'));
        MaxMAP = VTF3D_integrateMaxResponses(MAPcat,MinScale,MaxScale,alpha,Up2Orig,PowerFactor,QuantileTensorTHR,Msk,ProcessTensor);
        
    else
        disp('No *.mat files in the given PathFolder!');
    end
    
end
disp('VTrails: Integrating Multi-Resolution Filter Responses Across Scales -- COMPLETE!');

function [PathFolder,ScalesRange,alpha,Up2Orig,PowerFactor,QuantileTensorTHR,Msk,ProcessTensor] = configureInputs(Inputs)

% Initializing OPT with Default Values
PathFolder = [];
ScalesRange = [0.0,1.0];
Up2Orig = false;
alpha = 1;
PowerFactor = 1;
QuantileTensorTHR = 0.9;
Msk = [];
ProcessTensor = 'Yes';

if ~isempty(Inputs)
    for j = 1 : 2 :length(Inputs)
        
        switch upper(Inputs{j})
            case 'PATHFOLDER'
                if isempty(Inputs{j+1}) || exist(Inputs{j+1},'file') ~= 7
                    error('Invalid|Not Existing PathFoilder!');
                else
                    PathFolder = Inputs{j+1};
                end
            case 'SCALESRANGE'
                if isempty(Inputs{j+1}) || ~isvector(Inputs{j+1}) || ~isnumeric(Inputs{j+1}) || sum(isfinite(Inputs{j+1}))~=numel(Inputs{j+1}) || length(Inputs{j+1}) ~= 2
                    disp('Invalid|Not Given ScalesRange [min,max] array: -- Default: [0.0,1.0] applied!');
                else
                    assert(isnumeric(ScalesRange) && min(ScalesRange(:)) >= 0 && max(ScalesRange(:)) <= 1 ,'ScalesRange MUST be a 2-Elements Vector with range [0.0,1.0].');
                    ScalesRange = Inputs{j+1};
                end
            case 'ALPHA'
                if isempty(Inputs{j+1}) || ~isscalar(Inputs{j+1}) || ~isnumeric(Inputs{j+1}) || sum(isfinite(Inputs{j+1}))~=numel(Inputs{j+1})
                    disp('Invalid|Not Given alpha value: -- Default: 1.0 applied!');
                else
                    alpha = abs(Inputs{j+1});
                end
            case 'POWERFACTOR'
                if isempty(Inputs{j+1}) || ~isscalar(Inputs{j+1}) || ~isnumeric(Inputs{j+1}) || sum(isfinite(Inputs{j+1}))~=numel(Inputs{j+1})
                    disp('Invalid|Not Given PowerFactor: -- Default: 1.0 applied!');
                else
                    PowerFactor = abs(Inputs{j+1});
                end
            case 'QUANTILETENSORTHR'
                if isempty(Inputs{j+1}) || ~isscalar(Inputs{j+1}) || ~isnumeric(Inputs{j+1}) || sum(isfinite(Inputs{j+1}))~=numel(Inputs{j+1})
                    disp('Invalid|Not Given QuantileTensorTHR: -- Default: 0.9 applied!');
                else
                    QuantileTensorTHR = abs(Inputs{j+1});
                    QuantileTensorTHR(QuantileTensorTHR > 1 || QuantileTensorTHR <= 0) = 0.9;
                end
            case 'UP2ORIG'
                if isempty(Inputs{j+1}) || ~islogical(Inputs{j+1}) || ~isscalar(Inputs{j+1})
                    disp('Invalid|Not Given Up2Orig flag: -- Default: "false" applied!');
                else
                    Up2Orig = Inputs{j+1};
                end
            case 'MSK'
                if isempty(Inputs{j+1}) || ~islogical(Inputs{j+1})
                    disp('Invalid|Not Given Msk: -- Default: [] applied!');
                else
                    Msk = Inputs{j+1};
                end
            case 'PROCESSTENSOR'
                if isempty(Inputs{j+1}) || ~islogical(Inputs{j+1})
                    disp('Invalid|Not Given Msk: -- Default: "Yes" applied!');
                else
                    if Inputs{j+1}
                        ProcessTensor = 'Yes';
                    else
                        ProcessTensor = 'No';
                    end
                end
            otherwise
                disp('Given Unknown parameter(s). DEFAULT value(s) will be set.');
        end
        
    end
end

if isempty(PathFolder) || exist(PathFolder,'file') ~= 7
    error('Invalid|Not Existing PathFoilder!');
end

function MaxMAP = VTF3D_integrateMaxResponses(MAP,MinScale,MaxScale,alpha,Up2Orig,PowerFactor,QuantileTensorTHR,Msk,ProcessTensor)

[MaxMAP,NScales] = VTF3D_integrateRiemannianResponses(MAP,MinScale,MaxScale,alpha,PowerFactor,ProcessTensor);

clear MAP;

MaxMAP = VTF3D_convertTensorField2Euclidean(MaxMAP,NScales,Up2Orig,QuantileTensorTHR,Msk,ProcessTensor);


function [MaxMAP,processed] = VTF3D_integrateRiemannianResponses(MAP,MinScale,MaxScale,alpha,PowerFactor,ProcessTensor)

processed = 0;

for s = 1 : size(MAP,2)
    
    if MAP(s).GRID.Scale >= MinScale && MAP(s).GRID.Scale <= MaxScale
        
        processed = processed + 1;
        
        if processed == 1
            
            MaxMAP.IMG = MAP(s).IMG;
            
            Vmax = (MAP(s).CVM).^PowerFactor;
            
            % Apply Intensity Correction
            correction = MAP(s).BDM .* (1-MAP(s).BGM);
            correction(correction<0) = 0;
            
            dbm4gg = (2*MAP(s).BDM)-1;
            dbm4gg(dbm4gg<0) = 0;
            dbm4gg = dbm4gg.^0.5;
            
            Vmax = Vmax + correction;
            
            Vmax = Vmax .* (dbm4gg.*(MAP(s).BGM.^0.5));
            
            Vmax = (Vmax-min(Vmax(:)))/(max(Vmax(:))-min(Vmax(:)));

            % Tensor
            if strcmpi(ProcessTensor,'YES')
                MaxMAP.TFLE.T11 = MAP(s).CVM .* MAP(s).TFLE.T11;
                MaxMAP.TFLE.T12 = MAP(s).CVM .* MAP(s).TFLE.T12;
                MaxMAP.TFLE.T13 = MAP(s).CVM .* MAP(s).TFLE.T13;
                MaxMAP.TFLE.T22 = MAP(s).CVM .* MAP(s).TFLE.T22;
                MaxMAP.TFLE.T23 = MAP(s).CVM .* MAP(s).TFLE.T23;
                MaxMAP.TFLE.T33 = MAP(s).CVM .* MAP(s).TFLE.T33;
            else
                MaxMAP.TFLE = [];
            end
            
            % Assignment
            MaxMAP.CVM = Vmax;
            MaxMAP.BDM = MAP(s).BDM;
            MaxMAP.BGM = MAP(s).BGM;
            MaxMAP.GRID = MAP(s).GRID;
           
            MaxMAP.ACC = MAP(s).CVM;

            % Unorganised Seeds
            MaxMAP.US = false(size(MAP(s).CVM));
            MaxMAP.US(MAP(s).SSsIDX) = true;
            MaxMAP.GRID.alpha = alpha;
            
        else
            
            Vnew = (MAP(s).CVM).^PowerFactor;
            
            % Apply Intensity Correction
            correction = MAP(s).BDM .* (1-MAP(s).BGM);
            correction(correction<0) = 0;
            
            dbm4gg = (2*MAP(s).BDM)-1;
            dbm4gg(dbm4gg<0) = 0;
            dbm4gg = dbm4gg.^0.5;
            
            Vnew = Vnew + correction;

            Vnew = Vnew .* (dbm4gg.*(MAP(s).BGM.^0.5));
            
            Vnew = (Vnew-min(Vnew(:)))/(max(Vnew(:))-min(Vnew(:)));

            Vnew = alpha * Vnew;
            
            targetGRID = MAP(s).GRID;
            
            % Resampling (upsampling)
            [VMaxUpsmpl,BDupsmpl,BGupsmpl,TFLEMaxUpsmpl,ACCupsmpl] = VTF3D_upsampleMaxMAP(MaxMAP,targetGRID);
            
            Vmax_old = VMaxUpsmpl;
            
            % Evaluating Local Maxima
            [Vmax,~] = max(cat(4,VMaxUpsmpl,Vnew),[],4);
            
            integralcorrection = BDupsmpl .* (1-BGupsmpl);
            integralcorrection(integralcorrection<0) = 0;
            
            if MAP(s).GRID.Scale < 1
                gain = 10 * exp( 1/sqrt( ( 2 * pi * (1-MAP(s).GRID.Scale) ) ) );
            else
                gain = 1.5882e+14;
            end
            
            Vmax = Vmax_old + (gain * Vmax) + integralcorrection ;
            
            % Boundary Integration
            BDM = (2*MAP(s).BDM) - 1; % [-1,1]
            BDM = (1-BGupsmpl) .* (BDM + BDupsmpl);
            BDM = (BDM - min(BDM(:)))./(max(BDM(:)) - min(BDM(:)));

            % Background Integration
            BGM = (MAP(s).BGM - min(MAP(s).BGM(:)))./(max(MAP(s).BGM(:)) - min(MAP(s).BGM(:))); %[0,1]
            BGM = BGM + BGupsmpl;
            BGM = (BGM - min(BGM(:)))./(max(BGM(:)) - min(BGM(:)));

            % Tensor Integration
            if strcmpi(ProcessTensor,'YES')
                TFLEmax.T11 = TFLEMaxUpsmpl.T11 + (MAP(s).CVM .* MAP(s).TFLE.T11 );
                TFLEmax.T12 = TFLEMaxUpsmpl.T12 + (MAP(s).CVM .* MAP(s).TFLE.T12 );
                TFLEmax.T13 = TFLEMaxUpsmpl.T13 + (MAP(s).CVM .* MAP(s).TFLE.T13 );
                TFLEmax.T22 = TFLEMaxUpsmpl.T22 + (MAP(s).CVM .* MAP(s).TFLE.T22 );
                TFLEmax.T23 = TFLEMaxUpsmpl.T23 + (MAP(s).CVM .* MAP(s).TFLE.T23 );
                TFLEmax.T33 = TFLEMaxUpsmpl.T33 + (MAP(s).CVM .* MAP(s).TFLE.T33 );
            else
                TFLEmax = [];
            end

            % Accumulator Update
            ACC = ACCupsmpl + MAP(s).CVM;

            % Seeds
            SSsIDX = find(MaxMAP.US);
            [SSsX,SSsY,SSsZ] = ind2sub(size(MaxMAP.US),SSsIDX);
            
            SSsXup = round(SSsX./(MaxMAP.GRID.Scale / targetGRID.Scale)); SSsXup(SSsXup<1) = 1; SSsXup(SSsXup > size(MAP(s).CVM,1)) = size(MAP(s).CVM,1);
            SSsYup = round(SSsY./(MaxMAP.GRID.Scale / targetGRID.Scale)); SSsYup(SSsYup<1) = 1; SSsYup(SSsYup > size(MAP(s).CVM,2)) = size(MAP(s).CVM,2);
            SSsZup = round(SSsZ./(MaxMAP.GRID.Scale / targetGRID.Scale)); SSsZup(SSsZup<1) = 1; SSsZup(SSsZup > size(MAP(s).CVM,3)) = size(MAP(s).CVM,3);
            
            MaxMAP.US = false(size(MAP(s).CVM));
            MaxMAP.US(sub2ind(size(MaxMAP.US),SSsXup,SSsYup,SSsZup)) = true;
            
            % Assignment
            MaxMAP.CVM = Vmax;
            MaxMAP.TFLE = TFLEmax;
            MaxMAP.GRID = MAP(s).GRID;
            
            MaxMAP.GRID.alpha = alpha;
            
            MaxMAP.BDM = BDM;
            MaxMAP.BGM = BGM;

            MaxMAP.ACC = ACC;

            % Unorganised Seeds
            MaxMAP.US(MAP(s).SSsIDX) = true;
            
        end
        
    end
    
end

function MaxMAP = VTF3D_convertTensorField2Euclidean(MaxMAP,NScales,Up2Orig,QuantileTensorTHR,Msk,ProcessTensor)

if Up2Orig
    MaxMAP = resampleUp2OriginalSize(MaxMAP);
end

MaxMAP.CVM = (MaxMAP.CVM-min(MaxMAP.CVM(:)))./(max(MaxMAP.CVM(:)) - min(MaxMAP.CVM(:)));

% Smoothed Version for Seeds Alignment
MaxMAP.CVMsmth = single( VTF3D_GaussBlur3D( double(MaxMAP.CVM) , double(MaxMAP.GRID.sigma) ) );

% Recover the Euclidean Maximal Tensor Field
if ~isempty(Msk)
    [XMsk,YMsk,ZMsk] = ndgrid(1:size(Msk,1),1:size(Msk,2),1:size(Msk,3));
    [XMskR,YMskR,ZMskR] = ndgrid(MaxMAP.GRID.Xr,MaxMAP.GRID.Yr,MaxMAP.GRID.Zr);
    MskR = imerode(interpn(XMsk,YMsk,ZMsk,Msk,XMskR,YMskR,ZMskR,'nearest'),true([3 3 3]));
    TensorValidTHRSHLD = quantile(MaxMAP.CVM( MskR ),QuantileTensorTHR);
    TensorValidMSK = (MaxMAP.CVM > TensorValidTHRSHLD) & (MskR) ;
    
    MaxMAP.Barriers = ~ MskR; % Useful for Further Connectivity
    
else
    TensorValidTHRSHLD = quantile(MaxMAP.CVM(MaxMAP.CVM > 0),QuantileTensorTHR);
    TensorValidMSK = MaxMAP.CVM > TensorValidTHRSHLD;
    
    MaxMAP.Barriers = []; % Useful for Further Connectivity: NO BARRIERS!
    
end

MaxMAP.CVM = single(MaxMAP.CVM);

% NB: all the other tensors associated with a Maximal Vesselness Map lower
% than VmaxThreshold will be ISOTROPIC!
if strcmpi(ProcessTensor,'YES')
    MaxMAP.TF = VTF3D_retrieveFinalTVFfromAvgTVFLE(MaxMAP.TFLE, MaxMAP.ACC, TensorValidMSK , NScales );
    MaxMAP.TFValidityMSK = TensorValidMSK;
else
    MaxMAP.TF = [];
    MaxMAP.TFValidityMSK = [];
end

if Up2Orig
    MaxMAP.GRID.Scale = 1.0;
    MaxMAP.GRID.Ord = [];
    MaxMAP.GRID.Xr = single(1:1:MaxMAP.GRID.Xo);
    MaxMAP.GRID.Yr = single(1:1:MaxMAP.GRID.Yo);
    MaxMAP.GRID.Zr = single(1:1:MaxMAP.GRID.Zo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
function [VSupsmpl,BDupsmpl,BGupsmpl,TFLEupsmpl,ACCupsmpl] = VTF3D_upsampleMaxMAP(MAP,targetGRID)

[GRID.Xini,GRID.Yini,GRID.Zini] = ndgrid(MAP.GRID.Xr,MAP.GRID.Yr,MAP.GRID.Zr);
[GRID.Xfin,GRID.Yfin,GRID.Zfin] = ndgrid(targetGRID.Xr,targetGRID.Yr,targetGRID.Zr);

VSupsmpl = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.CVM ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');% Up-sampling using LINEAR interpolation
VSupsmpl = (VSupsmpl - min(VSupsmpl(:)))/(max(VSupsmpl(:)) - min(VSupsmpl(:)));

BDupsmpl = interpn(GRID.Xini,GRID.Yini,GRID.Zini, (2 * MAP.BDM) - 1 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');% Up-sampling using LINEAR interpolation
%BDupsmpl = (BDupsmpl - min(BDupsmpl(:)))/(max(BDupsmpl(:)) - min(BDupsmpl(:)));

BGupsmpl = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.BGM ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');% Up-sampling using LINEAR interpolation
BGupsmpl = (BGupsmpl - min(BGupsmpl(:)))/(max(BGupsmpl(:)) - min(BGupsmpl(:)));

if ~isempty(MAP.TFLE)
    TFLEupsmpl.T11 = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.TFLE.T11 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');
    TFLEupsmpl.T12 = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.TFLE.T12 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');
    TFLEupsmpl.T13 = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.TFLE.T13 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');
    TFLEupsmpl.T22 = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.TFLE.T22 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');
    TFLEupsmpl.T23 = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.TFLE.T23 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');
    TFLEupsmpl.T33 = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.TFLE.T33 ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');
else
    TFLEupsmpl = [];
end

ACCupsmpl = interpn(GRID.Xini,GRID.Yini,GRID.Zini, MAP.ACC ,GRID.Xfin,GRID.Yfin,GRID.Zfin,'linear');% Up-sampling using LINEAR interpolation

function TF = VTF3D_retrieveFinalTVFfromAvgTVFLE(TFLE , ACC , LogicalMap , enhanceFactor)
% Converting the Averaged Structure Tensors Field Across Scales into
% EUCLIDEAN Domain

PowerFactor = 2^ceil(log2((3 + enhanceFactor))) -1;

% Spatial Smoothing for Consistency
HSmoothLin = [.25/2 .75 .25/2];

T11r = TFLE.T11 ./ (ACC) ;
T12r = TFLE.T12 ./ (ACC) ;
T13r = TFLE.T13 ./ (ACC) ;
T22r = TFLE.T22 ./ (ACC) ;
T23r = TFLE.T23 ./ (ACC) ;
T33r = TFLE.T33 ./ (ACC) ;

% Smoothing is necessary to enforce spatial consistency of the Composed Multi-Scale Tensor Field
T11r(abs(T11r)<1e-6 | ~isfinite(T11r)) = 1e-6; T11r = VTF3D_MeanFilter3D(T11r,HSmoothLin);
T12r(abs(T12r)<1e-6 | ~isfinite(T12r)) = 1e-6; T12r = VTF3D_MeanFilter3D(T12r,HSmoothLin);
T13r(abs(T13r)<1e-6 | ~isfinite(T13r)) = 1e-6; T13r = VTF3D_MeanFilter3D(T13r,HSmoothLin);
T22r(abs(T22r)<1e-6 | ~isfinite(T22r)) = 1e-6; T22r = VTF3D_MeanFilter3D(T22r,HSmoothLin);
T23r(abs(T23r)<1e-6 | ~isfinite(T23r)) = 1e-6; T23r = VTF3D_MeanFilter3D(T23r,HSmoothLin);
T33r(abs(T33r)<1e-6 | ~isfinite(T33r)) = 1e-6; T33r = VTF3D_MeanFilter3D(T33r,HSmoothLin);

[El1,El2,El3,Ev1,Ev2,Ev3] = VTF3D_convert6LogEuclidComponentsTo3DTensor(T11r,T12r,T13r,T22r,T23r,T33r, LogicalMap );

El1 = ( (1./El1) ) .^ PowerFactor;
El2 = ( (1./El2) ) .^ PowerFactor; 
El3 = ( (1./El3) ) .^ PowerFactor; 

% Final Normalization s.t. the Product of the Eigenvalue Maps is equal to 1
ElrProd = (El1.*El2.*El3).^(1/3);

% Final Assignment s.t. the determinant is equal to 1.
El1 = El1 ./ ElrProd;
El2 = El2 ./ ElrProd;
El3 = El3 ./ ElrProd;

%%% Correcting for possible Degenerate Cases %%%
[El1,El2,El3,Ev1,Ev2,Ev3] = VTF3D_correctDegenerateTensors(El1,El2,El3,Ev1,Ev2,Ev3);

% Returning the Final Structure Tensors Field in EUCLID Domain
TF.El1 = El1;
TF.El2 = El2;
TF.El3 = El3;
TF.Ev1 = Ev1;
TF.Ev2 = Ev2;
TF.Ev3 = Ev3;

TF = VTF3D_makeTFsingle(TF);

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

function [El1,El2,El3,Ev1,Ev2,Ev3] = VTF3D_convert6LogEuclidComponentsTo3DTensor(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE,Msk)

w = 1./[1,sqrt(2),sqrt(2),1,sqrt(2),1];
%%%
% Dummy T fields:
T.Ix = zeros(size(T11LE)); T.Iy = zeros(size(T11LE)); T.Iz = zeros(size(T11LE));
% Real T fields:
T.Ixx = double( w(1) * T11LE );
T.Ixy = double( w(2) * T12LE );
T.Ixz = double( w(3) * T13LE );
T.Iyy = double( w(4) * T22LE );
T.Iyz = double( w(5) * T23LE );
T.Izz = double( w(6) * T33LE );
% Dummy T fields:
T.Inorm = zeros(size(T11LE)); T.IxN = zeros(size(T11LE)); T.IyN = zeros(size(T11LE)); T.IzN = zeros(size(T11LE));
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

if ~isstruct(I)
    Size = size(I);
    IMG = VTF3D_computeGradients2ndOrd_mex(I);
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

block_len = 1024;
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
   
   [Evl1 , Evl2 , Evl3 , Evc1 , Evc2 , Evc3 , idx_discard] = VTF3D_eigenDecompBlock3D_mex(VG,idx_block,SortByMagnitude,AbsoluteFlag);
   
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

function TFs_single = VTF3D_makeTFsingle(TFs)

TFs_single.El1 = single(TFs.El1);
TFs_single.El2 = single(TFs.El2);
TFs_single.El3 = single(TFs.El3);
TFs_single.Ev1 = single(TFs.Ev1);
TFs_single.Ev2 = single(TFs.Ev2);
TFs_single.Ev3 = single(TFs.Ev3);

function MaxMAP = resampleUp2OriginalSize(MaxMAP)

tgtGRID.Xr = 1:MaxMAP.GRID.Xo;
tgtGRID.Yr = 1:MaxMAP.GRID.Yo;
tgtGRID.Zr = 1:MaxMAP.GRID.Zo;
[MaxMAP.CVM,MaxMAP.BDM,MaxMAP.BGM,MaxMAP.TFLE,MaxMAP.ACC] = VTF3D_upsampleMaxMAP(MaxMAP,tgtGRID);

MaxMAP.BDM = (MaxMAP.BDM-min(MaxMAP.BDM(:)))/(max(MaxMAP.BDM(:))-min(MaxMAP.BDM(:)));
MaxMAP.BGM = (MaxMAP.BGM-min(MaxMAP.BGM(:)))/(max(MaxMAP.BGM(:))-min(MaxMAP.BGM(:)));

[GG.Xini,GG.Yini,GG.Zini] = ndgrid(MaxMAP.GRID.Xr,MaxMAP.GRID.Yr,MaxMAP.GRID.Zr);
[GG.Xfin,GG.Yfin,GG.Zfin] = ndgrid(tgtGRID.Xr,tgtGRID.Yr,tgtGRID.Zr);

MaxMAP.US = interpn(GG.Xini,GG.Yini,GG.Zini, MaxMAP.US ,GG.Xfin,GG.Yfin,GG.Zfin,'nearest');

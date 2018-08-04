function VTA3D_visualiseRiemannianVesselness(CVM,TF,TFValidityMSK)%
% Use this function to display the Sintesized Connected Vesselness Map and
% the Associated Tensor Field 
%
% Example:
%
%      VTA3D_showImageVoxelTensor(MAP.TF,MAP.CVM,find(US));   
%  OR
%      VTA3D_showImageVoxelTensor(MaxMAP.TF, MaxMAP.CVM, find(OS) );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' - Rendering: CVM ... [This may take a little while]');

Size = size(CVM);
figure;
set(gcf,'Color','w');

subplot(121);
h_cc = contourslice(CVM.^2,1:2:size(CVM,1),1:2:size(CVM,2),1:2:size(CVM,3));
for hh = 1 : length(h_cc)
    set(h_cc,'EdgeAlpha',0.25);
end
set(gca,'YDir','Reverse');
axis equal;
axis tight;
view(3);
title('Connected Vesselness Map (CVM)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(122);
%%% Setting the whole Domain Range to the 3D Plot
plot3(1,1,1);
hold on;
plot3(size(TF.El1,1),size(TF.El1,2),size(TF.El1,3));
set(gca,'YDir','Reverse');
axis equal;
hold on;
%%%
title('Synthesised Tensor Field (TF)');

disp(' - Rendering: TF ... [This may take a little while]');

seedslistIDX = find(TFValidityMSK);
if length(seedslistIDX) > 1000
    disp(' * Too Many Samples: Reducing TFValidityMSK for efficiency *');
    seedslistIDX = seedslistIDX(randperm(length(seedslistIDX),1000));
end

plot_thrshld = 0.01;

for j = 1 : length(seedslistIDX) % For each Voxel in the list (idx)
   
    [row,col,pln] = ind2sub(Size,seedslistIDX(j)); % Get the indexed position
    
    if CVM(row,col,pln) > plot_thrshld
        if ~isnan(TF.El1(row,col,pln)) && ~isinf(TF.El1(row,col,pln))
            El1 = 1./TF.El1(row,col,pln);
            El2 = 1./TF.El2(row,col,pln);
            El3 = 1./TF.El3(row,col,pln);
            Ev1 = squeeze(TF.Ev1(row,col,pln,:));
            Ev2 = squeeze(TF.Ev2(row,col,pln,:));
            Ev3 = squeeze(TF.Ev3(row,col,pln,:));
            
            % Normalize the Eigenvalues s.t. prod([El1,El2,El3]) = 1
            Elprod = prod([El1,El2,El3]).^(1/3);
            
            El1 = El1/Elprod;
            El2 = El2/Elprod;
            El3 = El3/Elprod;
            
            T = CVM(row,col,pln) * ( El1*(Ev1*Ev1') + El2*(Ev2*Ev2') + El3*(Ev3*Ev3') );
            
            % Plot Ellipsoid
            VTA3D_plotGaussianEllipsoid([col,row,pln],T);
            clear El1 El2 El3 Ev1 Ev2 Ev3 Elprod T;
            hold on;
        end
    end
    
end

function h = VTA3D_plotGaussianEllipsoid(m, C, sdwidth, npts, axh)
% PLOT_GAUSSIAN_ELLIPSOIDS plots 2-d and 3-d Gaussian distributions
%
% H = PLOT_GAUSSIAN_ELLIPSOIDS(M, C) plots the distribution specified by 
%  mean M and covariance C. The distribution is plotted as an ellipse (in 
%  2-d) or an ellipsoid (in 3-d).  By default, the distributions are 
%  plotted in the current axes. H is the graphics handle to the plotted 
%  ellipse or ellipsoid.
%
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD) uses SD as the standard deviation 
%  along the major and minor axes (larger SD => larger ellipse). By 
%  default, SD = 1. Note: 
%  * For 2-d distributions, SD=1.0 and SD=2.0 cover ~ 39% and 86% 
%     of the total probability mass, respectively. 
%  * For 3-d distributions, SD=1.0 and SD=2.0 cover ~ 19% and 73%
%     of the total probability mass, respectively.
%  
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD, NPTS) plots the ellipse or 
%  ellipsoid with a resolution of NPTS (ellipsoids are generated 
%  on an NPTS x NPTS mesh; see SPHERE for more details). By
%  default, NPTS = 50 for ellipses, and 20 for ellipsoids.
%
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD, NPTS, AX) adds the plot to the
%  axes specified by the axis handle AX.
%
% Examples: 
% -------------------------------------------
%  % Plot three 2-d Gaussians
%  figure; 
%  h1 = plot_gaussian_ellipsoid([1 1], [1 0.5; 0.5 1]);
%  h2 = plot_gaussian_ellipsoid([2 1.5], [1 -0.7; -0.7 1]);
%  h3 = plot_gaussian_ellipsoid([0 0], [1 0; 0 1]);
%  set(h2,'color','r'); 
%  set(h3,'color','g');
% 
%  % "Contour map" of a 2-d Gaussian
%  figure;
%  for sd = [0.3:0.4:4],
%    h = plot_gaussian_ellipsoid([0 0], [1 0.8; 0.8 1], sd);
%  end
%
%  % Plot three 3-d Gaussians
%  figure;
%  h1 = plot_gaussian_ellipsoid([1 1  0], [1 0.5 0.2; 0.5 1 0.4; 0.2 0.4 1]);
%  h2 = plot_gaussian_ellipsoid([1.5 1 .5], [1 -0.7 0.6; -0.7 1 0; 0.6 0 1]);
%  h3 = plot_gaussian_ellipsoid([1 2 2], [0.5 0 0; 0 0.5 0; 0 0 0.5]);
%  set(h2,'facealpha',0.6);
%  view(129,36); set(gca,'proj','perspective'); grid on; 
%  grid on; axis equal; axis tight;
% -------------------------------------------
% 
%  Gautam Vallabha, Sep-23-2007, Gautam.Vallabha@mathworks.com
%
%  Revision 1.0, Sep-23-2007
%    - File created
%  Revision 1.1, 26-Sep-2007
%    - NARGOUT==0 check added.
%    - Help added on NPTS for ellipsoids

if ~exist('sdwidth', 'var'), sdwidth = 1; end
if ~exist('npts', 'var'), npts = []; end
if ~exist('axh', 'var'), axh = gca; end

if numel(m) ~= length(m), 
    error('M must be a vector'); 
end
if ~( all(numel(m) == size(C)) )
    error('Dimensionality of M and C must match');
end
if ~(isscalar(axh) && ishandle(axh) && strcmp(get(axh,'type'), 'axes'))
    error('Invalid axes handle');
end

set(axh, 'nextplot', 'add');

switch numel(m)
   case 2, h=VTA3D_show2d(m(:),C,sdwidth,npts,axh);
   case 3, h=VTA3D_show3d(m(:),C,sdwidth,npts,axh);
   otherwise
      error('Unsupported dimensionality');
end

if nargout==0,
    clear h;
end

function h = VTA3D_show2d(means, C, sdwidth, npts, axh)
if isempty(npts), npts=50; end
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
h = plot(bp(1,:), bp(2,:), '-b', 'parent', axh);

function h = VTA3D_show3d(means, C, sdwidth, npts, axh)
if isempty(npts), npts= 8 ; end %20
[x,y,z] = sphere(npts);
ap = [x(:) y(:) z(:)]';
[v,d]=eig(C); 
if any(d(:) < 0)
   fprintf('warning: negative eigenvalues\n');
   d = max(d,0);
end
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
xp = reshape(bp(1,:), size(x));
yp = reshape(bp(2,:), size(y));
zp = reshape(bp(3,:), size(z));
h = surf(axh, xp,yp,zp,'FaceColor',[0 0 1],'FaceAlpha',0.05,'EdgeAlpha',0.05); %0.1,0.1

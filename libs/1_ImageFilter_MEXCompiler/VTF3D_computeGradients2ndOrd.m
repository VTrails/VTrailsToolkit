function VG = VTF3D_computeGradients2ndOrd(img3D)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VG = struct( 'Ix',NaN(size(img3D)), 'Iy',NaN(size(img3D)), 'Iz',NaN(size(img3D)),...
            'Ixx',NaN(size(img3D)),'Ixy',NaN(size(img3D)),'Ixz',NaN(size(img3D)),...
            'Iyy',NaN(size(img3D)),'Iyz',NaN(size(img3D)),'Izz',NaN(size(img3D)),...
            'Inorm',NaN(size(img3D)),'IxN',NaN(size(img3D)),'IyN',NaN(size(img3D)),'IzN',NaN(size(img3D)));

noiseaddgain = 1e-6;

% Computing Gradient and Hessian components for the EigenDecomposition
[Ix,Iy,Iz] = gradient( img3D + noiseaddgain*rand(size(img3D)) ); % Including additional Noise!
[Ixx,Ixy,Ixz] = gradient(Ix);
[ ~ ,Iyy,Iyz] = gradient(Iy);
[ ~ , ~ ,Izz] = gradient(Iz);

% Struct
VG.Ix = Ix; VG.Iy = Iy; VG.Iz = Iz;
VG.Ixx = Ixx; VG.Ixy = Ixy; VG.Ixz = Ixz;
VG.Iyy = Iyy; VG.Iyz = Iyz; VG.Izz = Izz;
% Normalisation
VG.Inorm = sqrt((VG.Ix.^2 + VG.Iy.^2 + VG.Iz.^2));
VG.IxN = VG.Ix ./ VG.Inorm;
VG.IyN = VG.Iy ./ VG.Inorm;
VG.IzN = VG.Iz ./ VG.Inorm;
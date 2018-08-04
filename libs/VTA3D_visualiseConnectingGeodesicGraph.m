function hfig = VTA3D_visualiseConnectingGeodesicGraph(exploredminPath,GGM)%
%
% Visualise the whole Over-Connected Geodesic Vascular Graph Inferred with
% the Connectivity Paradigm
%
% INPUTS:
% 
% - exploredminPath: struct as in DATA.exploredminPath
%                    (see 'VTrailsConnectingGeodesicGraph3D_MAIN.m' @OUTPUTS)
% - GGM: Geodesic Graph Matrix as in DATA.GGM
%        (see 'VTrailsConnectingGeodesicGraph3D_MAIN.m' @OUTPUTS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxval = max(nonzeros(GGM));
minval = min(nonzeros(GGM));

maxLineWidth = 0.5;

hfig = figure;
set(hfig,'Color','w');
title('Over-Connected Geodesic Vascular Graph');
hold on;

for aa = 1 : length(exploredminPath)

    GGMval = GGM(exploredminPath(aa).pairs(1),exploredminPath(aa).pairs(2));
    GGMval = (GGMval - minval)/(maxval-minval);
    GGMval(GGMval<0) = eps;
    
    hp = plot3(exploredminPath(aa).sequence(:,1), exploredminPath(aa).sequence(:,2), exploredminPath(aa).sequence(:,3),'LineWidth',maxLineWidth*((1-GGMval).^2)+eps);
    
    hp.Color = [hp.Color,0.5];
    hold on;
    
end

set(gca,'YDir','Reverse');
axis equal;
axis tight;
view(3);

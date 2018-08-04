function hfig = VTA3D_visualiseGeodesicMinimumSpanningTree(GeodesicMSTs,CVM)%
%
% Visualise the Minimum Spanning Tree(s) of an Over-Connected Geodesic
% Vascular Graph
%
% INPUTS:
% 
% - GeodesicMSTs: output struct generated with:
%                 'VTrailsRefine3D_ConnectedGraph2VTrailsTree_MAIN',
%                 aka VTRAILS.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stefano Moriconi, July-2017, vtrailstoolkit@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hfig = figure;
set(hfig,'Color','w');
title('Geodesic Minimum Spanning Tree(s) of the Vascular Graph - VTrails');
hold on;

for aa = 1 : length(GeodesicMSTs)

    plot3(GeodesicMSTs(aa).CGPathContinuous(:,1),...
          GeodesicMSTs(aa).CGPathContinuous(:,2),...
          GeodesicMSTs(aa).CGPathContinuous(:,3),'LineWidth',3);
    
      hold on;
    
end

h_cc = contourslice(CVM.^2,1:2:size(CVM,1),1:2:size(CVM,2),1:2:size(CVM,3));
for hh = 1 : length(h_cc)
    set(h_cc,'EdgeAlpha',0.1);
end

set(gca,'YDir','Reverse');
axis equal;
axis tight;
view(3);
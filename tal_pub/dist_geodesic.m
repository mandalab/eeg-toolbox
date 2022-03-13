 function [d, fastmarch, subsurf] = dist_geodesic(gsurf, ndx_seed, dest_ndx)
    % Computes geodesic distance (distance along the surface) from one surface point to a specific point
    %
    % USAGE: 
    %   d = dist_geodesic(gsurf, ndx_seed, dest_ndx)
    %
    % Inputs:
    %   gsurf - Gifti surface object (with vertices and faces fields)
    %   ndx_seed - Index into surface mesh
    %
    % Optional Input:
    %   dest_ndx - if you want the distance to a specific node, pass that node here (faster)
    %
    % OUTPUT:
    %   d   - Distance between ndx_seed and dest_ndx (in millimeters)
    %   fastmarch - a struct with the outputs from fastmarching (D, S, and Q)
    %   subsurf - a struct with the outputs from subsurface (sv, sf, si, vidxs, fidxs)
    %
    % Example:
    %   d = dist_geodesic(bp.surfaces.pial_lh, mesh_ndx)
    %   d = dist_geodesic(bp.surfaces.pial_lh, mesh_ndx1, mesh_ndx2)
    %
    %
    % REVISION HISTORY
    %   MST 09/28 - Created

    if nargin < 3
        dest_ndx = [];
        fastmarch = [];
        subsurf = [];
    end

    V = double(gsurf.vertices)';
    F = double(gsurf.faces)';

    if isempty(dest_ndx);
        d = surfing_dijkstradist(V, F, ndx_seed);
    else
        n2f = surfing_nodeidxs2faceidxs_custom(gsurf.faces', gsurf.vertices');

        d_euclid = norm(double(gsurf.vertices([ndx_seed; dest_ndx],:)));
        PADCLIP = 1.05;
        
        % The following clips the surface using euclidean distance (with a 5 percent padding) for speed improvement
        % and then computes geodesic distance with fast-marching
        [mesh_ndxs, dists, ~, ~, fastmarch, subsurf] = surfing_circleROI_custom(V, F, ndx_seed, d_euclid*PADCLIP, 'geodesic',n2f);
        d = dists(mesh_ndxs == dest_ndx);
    end

end
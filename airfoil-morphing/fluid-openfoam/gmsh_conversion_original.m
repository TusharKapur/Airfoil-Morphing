function dattostl (mainParam, flapParam,InputDirectory,OutputDirectory,MAC)
% system(['LD_PRELOAD="/lib64/libstdc++.so.6"']);

% system(['rm ',directions.path.simulation,' airfoil.geo airfoil.msh']);
% system(['touch ',directions.path.simulation,' airfoil.geo']);

%a_points = readtable("bladenlf13.dat")
size_airfoil = "_VARsize-airfoil_";    % Size of mesh near main airfoil
size_airfoil_flap = "_VARsize-airfoil-flap_";    % Size of mesh near the flap airfoil
size_boundary = "_VARsize-boundary_";    % Size of mesh near boundary
extrude_thickness = "_VARextrude-thickness_";    % Thickness of extruded mesh
no_layers = "_VARno-layers_";  % Number of layers in extrude-direction

rot_x_main = "_VARrot-x-main_"; % Centre of rotation for main airfoil (X - coordinate)
rot_y_main = "_VARrot-y-main_"; % Centre of rotation for main airfoil (Y - coordinate)
rot_z_main = "_VARrot-z-main_"; % Centre of rotation for main airfoil (Z - coordinate)
rot_ang_main = "_VARrot-ang-main_";   % Angle of rotation for main airfoil (unit: radians)
rot_x_flap = "_VARrot-x-flap_"; % Centre of rotation for flap airfoil (X - coordinate)
rot_y_flap = "_VARrot-y-flap_"; % Centre of rotation for flap airfoil (Y - coordinate)
rot_z_flap = "_VARrot-z-flap_"; % Centre of rotation for flap airfoil (Z - coordinate)
rot_ang_flap = "_VARrot-ang-flap_";   % Angle of rotation for flap airfoil (unit: radians)

tra_x_main = "_VARtra-x-main_"; % Translation for main airfoil (X - coordinate)
tra_y_main = "_VARtra-y-main_"; % Translation for main airfoil (Y - coordinate)
tra_z_main = "_VARtra-z-main_"; % Translation for main airfoil (Z - coordinate)

tra_x_flap = "_VARtra-x-flap_"; % Translation for flap airfoil (X - coordinate)
tra_y_flap = "_VARtra-y-flap_"; % Translation for flap airfoil (Y - coordinate)
tra_z_flap = "_VARtra-z-flap_"; % Translation for flap airfoil (Z - coordinate)

no_nodes_by_edge = "_VARno-nodes-by-edge_";  % Number of nodes in a particular edge

thresh_max_dist = "_VARthresh-max-dist_";  % Max distance of the threshold
thresh_min_dist = "_VARthresh-min-dist_";  % Min distance of the threshold
thresh_sigmoid = "_VARthresh-sigmoid_"; % Interpolation between mesh inside and outside threshold
thresh_ifield = "_VARthresh-ifield_";

thresh_elem_far = "_VARsize-boundary_";    % Size of elements far from the airfoil
thresh_elem_near = "_VARsize-airfoil_";    % Size of elements near the airfoil

bd_layer_quads = "_VARbd-layer-quads_"; % Quadrilateral boundary layer elements
bd_layer_size_far = "_VARsize-boundary_";  % Size of elements far from the airfoil
bd_layer_size_near = "_VARbd-layer-size-near_";   % Thickness of boundary layer at the airfoil
bd_layer_aniso = "_VARbd-layer-aniso_";
bd_layer_intersect = "_VARbd-layer-intersect_";
bd_layer_growth_ratio = "_VARbd-layer-growth-ratio_";   % Growth ratio of boundary layers

mesh_algo = "_VARmesh-algo_";  % Choice of mesh algorithm

a_points = readtable(InputDirectory+"/"+mainParam.airfoilName+".dat");
a_points_flap = readtable(InputDirectory+"/"+flapParam.airfoilName+".dat");

% a_points = a_points(2:end,:);
airfoil_main_size = numel(a_points)/2
airfoil_flap_size = numel(a_points_flap)/2

system(['cd ',OutputDirectory,' && rm airfoil.geo']);
system(['cd ',OutputDirectory,' && touch airfoil.geo']);
airfoil_coord_file_name = append(OutputDirectory,'/','airfoil.geo');
fid = fopen(airfoil_coord_file_name, 'a')

fprintf(fid, "Geometry.OldNewReg=0; \n");   % For sequential ordering
fprintf(fid, "General.ExpertMode=1; \n");   % To disable warnings
fprintf(fid, "size_airfoil = %s; \n", size_airfoil);
fprintf(fid, "size_airfoil_flap = %s; \n", size_airfoil_flap);
fprintf(fid, "size_boundary = %s; \n", size_boundary);

pt_no = 0;
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {-43, 50, 0, size_boundary}; \n", pt_no);
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {-43, -49, 0, size_boundary}; \n", pt_no);
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {56, -49, 0, size_boundary}; \n", pt_no);
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {56, 50, 0, size_boundary}; \n", pt_no);
% pt_no = pt_no + 1;
% fprintf(fid, "Point(%i) = {-0.3, 0.4, 0.1, size_boundary}; \n", pt_no);
% pt_no = pt_no + 1;
% fprintf(fid, "Point(%i) = {-0.3, -0.4, 0.1, size_boundary}; \n", pt_no);
% pt_no = pt_no + 1;
% fprintf(fid, "Point(%i) = {2.0, -0.4, 0.1, size_boundary}; \n", pt_no);
% pt_no = pt_no + 1;
% fprintf(fid, "Point(%i) = {2.0, 0.4, 0.1, size_boundary}; \n", pt_no);

line_no = 0;
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 1, 2);
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 2, 3);
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 3, 4);
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 4, 1);

% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 5, 6);
% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 6, 7);
% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 7, 8);
% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 8, 5);

% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 1, 5);
% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 2, 6);
% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 3, 7);
% line_no = line_no + 1;
% fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 4, 8);

% line_no = 0;
% spline_main = [];
% spline_flap = [];

checkpoint_airfoil_main = pt_no + 1;
for n = 1:size(a_points(:,1))
    pt_no = pt_no + 1;
    % line_no  = line_no + 1;
    % length(string(a_points(line,:)))
    fprintf(fid, "Point(%i) = {%f, %f, 0, size_airfoil}; \n", pt_no, table2array(a_points(n,1)), table2array(a_points(n,2)));

    % if line_no > 3
        % pt_no = line_no - 3;
        % spline_main = [spline_main pt_no];
    % end
    %spline[1:int(line_no/2)-2] = reversed(spline[1:int(line_no/2)-2]);
end

checkpoint_airfoil_flap = pt_no + 1;
for n = 1:size(a_points_flap(:,1))
    pt_no = pt_no + 1;
    % line_no  = line_no + 1;
    % length(string(a_points(line,:)))
    fprintf(fid, "Point(%i) = {%f, %f, 0, size_airfoil_flap}; \n", pt_no, table2array(a_points_flap(n,1)), table2array(a_points_flap(n,2)));

    % if line_no > 3
        % pt_no = line_no - 3;
        % spline_flap = [spline_flap pt_no];
    % end
    %spline[1:int(line_no/2)-2] = reversed(spline[1:int(line_no/2)-2]);
end

% for n = 1:size(a_points(:,1))
%     line_no  = line_no + 1;
%     fprintf(fid, "Point(%i) = {%f, %f, 0, size_airfoil}; \n", n, table2array(a_points(n,1)), table2array(a_points(n,2)));

%     % if line_no > 3
%         % pt_no = line_no - 3;
%         pt_no = line_no;
%         spline = [spline pt_no];
%     % end
%     %spline[1:int(line_no/2)-2] = reversed(spline[1:int(line_no/2)-2]);
% end

% Rotate
% fprintf(fid, "Rotate {{0, 0, 1}, {-0.1, 0, 0}, Pi/8} {\n");
% for n = 1:size(a_points(:,1))
%     fprintf(fid, "Point{%i}; ", n);
% end
% fprintf(fid, "\n} \n");



pt_no_temp = checkpoint_airfoil_main;
line_airfoil_main = line_no + 1;
for n = 2:size(a_points(:,1))
    line_no = line_no + 1;
    fprintf(fid, "Line(%i) ={%i, %i}; \n", line_no, pt_no_temp, pt_no_temp+1);
    pt_no_temp = pt_no_temp + 1;
end
line_no = line_no + 1;
fprintf(fid, "Line(%i) ={%i, %i}; \n", line_no, pt_no_temp, checkpoint_airfoil_main);

pt_no_temp = checkpoint_airfoil_flap;
line_airfoil_flap = line_no + 1;
for n = 2:size(a_points_flap(:,1))
    line_no = line_no + 1;
    fprintf(fid, "Line(%i) ={%i, %i}; \n", line_no, pt_no_temp, pt_no_temp+1);
    pt_no_temp = pt_no_temp + 1;
end
line_no = line_no + 1;
fprintf(fid, "Line(%i) ={%i, %i}; \n", line_no, pt_no_temp, checkpoint_airfoil_flap);

curve_loop_no = 1;
fprintf(fid, "Curve Loop(%i) = {%i, %i, %i, %i}; \n", curve_loop_no, 1, 2, 3, 4);

curve_loop_no = curve_loop_no + 1;
fprintf(fid, "Curve Loop(%i) = {", curve_loop_no);
line_no_temp = line_airfoil_main;
for n = 2:size(a_points(:,1))
    fprintf(fid, "%i, ", line_no_temp);
    line_no_temp = line_no_temp + 1;
end
fprintf(fid, "%i}; \n", line_no_temp);

curve_loop_no = curve_loop_no + 1;
fprintf(fid, "Curve Loop(%i) = {", curve_loop_no);
line_no_temp = line_airfoil_flap;
for n = 2:size(a_points_flap(:,1))
    fprintf(fid, "%i, ", line_no_temp);
    line_no_temp = line_no_temp + 1;
end
fprintf(fid, "%i}; \n", line_no_temp);

surface_no = 1;
fprintf(fid, "Plane Surface(%i) = {", surface_no);
for n = 1:curve_loop_no-1
    fprintf(fid, "%i, ", n);
end
fprintf(fid, "%i}; \n", n+1);

fprintf(fid, "Extrude {0, 0, %s} {\n Surface{1}; Layers {%s}; Recombine; \n } \n", extrude_thickness, no_layers);

fprintf(fid, 'Physical Surface("inlet", 2) = {3}; \n Physical Surface("outlet", 3) = {5}; \n Physical Surface("topAndBottom", 4) = {6, 4}; \n Physical Surface("back", 5) = {1}; \n Physical Surface("front", 6) = {%i}; \n', 7+airfoil_main_size+airfoil_flap_size);

airfoil_main_surface_list = '';
airfoil_main_curve_list = '';
for n = 7:7+airfoil_main_size-2
    airfoil_main_surface_list = [airfoil_main_surface_list num2str(n,'%i')];
    airfoil_main_surface_list = [airfoil_main_surface_list ', '];
    airfoil_main_curve_list = [airfoil_main_curve_list num2str(n-2,'%i')];
    airfoil_main_curve_list = [airfoil_main_curve_list ', '];
end
airfoil_main_surface_list = [airfoil_main_surface_list num2str(n+1,'%i')];
airfoil_main_curve_list = [airfoil_main_curve_list num2str(n-1,'%i')];

airfoil_flap_surface_list = '';
airfoil_flap_curve_list = '';
for n = 7+airfoil_main_size:7+airfoil_main_size+airfoil_flap_size-2
    airfoil_flap_surface_list = [airfoil_flap_surface_list num2str(n,'%i')];
    airfoil_flap_surface_list = [airfoil_flap_surface_list ', '];
    airfoil_flap_curve_list = [airfoil_flap_curve_list num2str(n-2,'%i')];
    airfoil_flap_curve_list = [airfoil_flap_curve_list ', '];
end
airfoil_flap_surface_list = [airfoil_flap_surface_list num2str(n+1,'%i')];
airfoil_flap_curve_list = [airfoil_flap_curve_list num2str(n-1,'%i')];

% airfoil_main_surface_list
% fprintf('%s \n', airfoil_main_surface_list)

fprintf(fid, 'Physical Surface("bladenlf26_bladenlf26", 7) = {%s}; \n', airfoil_main_surface_list);

fprintf(fid, 'Physical Surface("bladenlf26_bladenlf26FLAP", 8) = {%s}; \n', airfoil_flap_surface_list');

% Transformations (rotation, translation, scaling)
% Rotation
fprintf(fid, "Rotate {{0, 0, 1}, {%s, %s, %s}, %s} {\n Surface{%s};} \n", rot_x_main, rot_y_main, rot_z_main, rot_ang_main, airfoil_main_surface_list);
fprintf(fid, "Rotate {{0, 0, 1}, {%s, %s, %s}, %s} {\n Surface{%s};} \n", rot_x_flap, rot_y_flap, rot_z_flap, rot_ang_flap, airfoil_flap_surface_list);

% Translation
fprintf(fid, "Translate {%s, %s, %s} {\n Surface{%s};} \n", tra_x_main, tra_y_main, tra_z_main, airfoil_main_surface_list);
fprintf(fid, "Translate {%s, %s, %s} {\n Surface{%s};} \n", tra_x_flap, tra_y_flap, tra_z_flap, airfoil_flap_surface_list);

% % Scaling - Currently not working
% fprintf(fid, "Dilate {{0, 0, 0}, {1, 1, 0}} {\n Surface{%s};} \n", airfoil_main_surface_list);
% fprintf(fid, "Dilate {{0, 0, 0}, {1, 1, 0}} {\n Surface{%s};} \n", airfoil_flap_surface_list);

fprintf(fid, 'Physical Volume("domain", 1) = {1}; \n');

fprintf(fid, "Field[1] = Distance; \n");
fprintf(fid, "Field[1].EdgesList = {%s, %s}; \n", airfoil_main_curve_list, airfoil_flap_curve_list);
fprintf(fid, "Field[1].NNodesByEdge = %s; \n \n", no_nodes_by_edge);

% fprintf(fid, "Field[2] = Threshold; \n");   % Refinement close to the airfoils
% fprintf(fid, "Field[2].DistMax = 0.05; \n");
% fprintf(fid, "Field[2].DistMin = 0.05; \n");
% fprintf(fid, "Field[2].IField = 1; \n");
% fprintf(fid, "Field[2].LcMax = size_boundary; \n");
% fprintf(fid, "Field[2].LcMin = size_airfoil; \n \n");

fprintf(fid, "Field[3] = Threshold; \n");   % Refinement close to the airfoils
fprintf(fid, "Field[3].DistMax = %s; \n", thresh_max_dist);
fprintf(fid, "Field[3].DistMin = %s; \n", thresh_min_dist);
fprintf(fid, "Field[3].Sigmoid = %s; \n", thresh_sigmoid);
fprintf(fid, "Field[3].IField = %s; \n", thresh_ifield);
fprintf(fid, "Field[3].LcMax = %s; \n", thresh_elem_far);
fprintf(fid, "Field[3].LcMin = %s; \n", thresh_elem_near);

% fprintf(fid, "Field[4] = Box; \n"); % Box refinement - not working
% fprintf(fid, "Field[4].Thickness = (size_airfoil*10 + size_boundary)/2; \n");
% fprintf(fid, "Field[4].VIn = size_airfoil*10; \n");
% fprintf(fid, "Field[4].VOut = size_boundary; \n");
% fprintf(fid, "Field[4].XMax = 2; \n");
% fprintf(fid, "Field[4].XMin = -0.3; \n");
% fprintf(fid, "Field[4].YMax = 0.4; \n");
% fprintf(fid, "Field[4].YMin = -0.4; \n");
% fprintf(fid, "Field[4].YMax = 0.4; \n");

fprintf(fid, "Field[5] = BoundaryLayer; \n");
fprintf(fid, "Field[5].EdgesList = {%s, %s}; \n", airfoil_main_curve_list, airfoil_flap_curve_list);
fprintf(fid, "Field[5].Quads = %s; \n", bd_layer_quads);
fprintf(fid, "Field[5].hfar = %s; \n", bd_layer_size_far);
fprintf(fid, "Field[5].hwall_n = %s; \n", bd_layer_size_near);
fprintf(fid, "Field[5].AnisoMax = %s; \n", bd_layer_aniso);
fprintf(fid, "Field[5].IntersectMetrics = %s; \n", bd_layer_intersect);
% fprintf(fid, "Field[4].thickness = 0.00029; \n");
fprintf(fid, "Field[5].ratio = %s; \n", bd_layer_growth_ratio);
fprintf(fid, "BoundaryLayer Field = 5; \n");

fprintf(fid, "Field[6] = Min; \n");
fprintf(fid, "Field[6].FieldsList = {3, 5}; \n");
fprintf(fid, "Background Field = 6; \n");

fprintf(fid, "Mesh.Algorithm = %s; \n", mesh_algo);

% system(['gmsh airfoil.geo -3 airfoil.msh -format msh2']);

fclose('all');
end

size_airfoil = 1.5E-00;    % Size of mesh near main airfoil
size_boundary = 1.5E-00;    % Size of mesh near boundary
extrude_thickness = 1.0;    % Thickness of extruded mesh
no_layers = 1;

a_points = readtable("airfoil.dat");

airfoil_size = numel(a_points)/2

system(['mv mesh.geo mesh-previous.geo']);
system(['touch mesh.geo']);
airfoil_coord_file_name = append('mesh.geo');
fid = fopen('mesh.geo', 'a')

fprintf(fid, "Geometry.OldNewReg=0; \n");   % For sequential ordering
fprintf(fid, "General.ExpertMode=1; \n");   % To disable warnings
fprintf(fid, "size_airfoil = %s; \n", size_airfoil);
fprintf(fid, "size_boundary = %s; \n", size_boundary);

pt_no = 0;
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {-20, 10, 0, size_airfoil}; \n", pt_no);
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {-20, -10, 0, size_airfoil}; \n", pt_no);
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {30, -10, 0, size_airfoil}; \n", pt_no);
pt_no = pt_no + 1;
fprintf(fid, "Point(%i) = {30, 10, 0, size_airfoil}; \n", pt_no);

line_no = 0;
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 1, 2);
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 2, 3);
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 3, 4);
line_no = line_no + 1;
fprintf(fid, "Line(%i) = {%i, %i}; \n", line_no, 4, 1);

checkpoint_airfoil = pt_no + 1;
for n = 1:size(a_points(:,1))
    pt_no = pt_no + 1;
    fprintf(fid, "Point(%i) = {%f, %f, 0, size_airfoil}; \n", pt_no, table2array(a_points(n,1)), table2array(a_points(n,2)));
end

pt_no_temp = checkpoint_airfoil;
line_airfoil = line_no + 1;
for n = 2:size(a_points(:,1))
    line_no = line_no + 1;
    fprintf(fid, "Line(%i) ={%i, %i}; \n", line_no, pt_no_temp, pt_no_temp+1);
    pt_no_temp = pt_no_temp + 1;
end
line_no = line_no + 1;
fprintf(fid, "Line(%i) ={%i, %i}; \n", line_no, pt_no_temp, checkpoint_airfoil);

curve_loop_no = 1;
fprintf(fid, "Curve Loop(%i) = {%i, %i, %i, %i}; \n", curve_loop_no, 1, 2, 3, 4);

curve_loop_no = curve_loop_no + 1;
fprintf(fid, "Curve Loop(%i) = {", curve_loop_no);
line_no_temp = line_airfoil;
for n = 2:size(a_points(:,1))
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

fprintf(fid, "Extrude {0, 0, 1.0} {\n Surface{1}; Layers {1}; Recombine; \n } \n");

fprintf(fid, 'Physical Surface("inlet", 2) = {3}; \n Physical Surface("outlet", 3) = {5}; \n Physical Surface("upperWall", 4) = {6}; \n Physical Surface("lowerWall", 5) = {4}; \n Physical Surface("frontAndBack", 6) = {1, %i}; \n', 7+airfoil_size);

airfoil_surface_list = '';
airfoil_curve_list = '';
for n = 7:7+airfoil_size-2
    airfoil_surface_list = [airfoil_surface_list num2str(n,'%i')];
    airfoil_surface_list = [airfoil_surface_list ', '];
    airfoil_curve_list = [airfoil_curve_list num2str(n-2,'%i')];
    airfoil_curve_list = [airfoil_curve_list ', '];
end
airfoil_surface_list = [airfoil_surface_list num2str(n+1,'%i')];
airfoil_curve_list = [airfoil_curve_list num2str(n-1,'%i')];

fprintf(fid, 'Physical Surface("flap", 7) = {%s}; \n', airfoil_surface_list);

% Transformations (rotation, translation, scaling)
% % Rotation
% fprintf(fid, "Rotate {{0, 0, 1}, {%s, %s, %s}, %s} {\n Surface{%s};} \n", rot_x_main, rot_y_main, rot_z_main, rot_ang_main, airfoil_surface_list);

% % Translation
% fprintf(fid, "Translate {%s, %s, %s} {\n Surface{%s};} \n", tra_x_main, tra_y_main, tra_z_main, airfoil_surface_list);

% % Scaling
% fprintf(fid, "Dilate {{0, 0, 0}, {1, 1, 0}} {\n Surface{%s};} \n", airfoil_surface_list);

fprintf(fid, 'Physical Volume("domain", 1) = {1}; \n');

% system(['gmsh mesh.geo -3 mesh.msh -format msh2']);

fclose('all');
% end

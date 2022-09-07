function dattostl (InputParam,InputDirectory,OutputDirectory,MAC)
%%  change format of .dat into format of .stl with openSCAD
%   dattostl (InputParam,InputDirectory,OutputDirectory)
%   changes the inputParam.airfoilName from .dat into .stl
%   Also translates and rotates as in the input file specified
%   Also moves file from InputDirectory into Output Directory

a_points = readtable(InputDirectory+"/"+InputParam.airfoilName+".dat");
% a_points.Var1 = 1000 .* a_points.Var1;
% a_points.Var2 = 1000 * a_points.Var2;     % millimeters or meters? OpenFoam is in meters!  

n = numel(a_points)/2;

airfoilXYZX= zeros(1,n);
airfoilXYZX = string (airfoilXYZX);

fileNameNDirectTXT = append(OutputDirectory,'/',InputParam.airfoilName,'.txt');
fileNameNDirectSCAD = append(OutputDirectory,'/',InputParam.airfoilName,'.scad');
fileNameNDirectSTL = append(OutputDirectory,'/',InputParam.airfoilName,'.stl');

for i=1:n;
    
    if i<(n)
        
        airfoilXYZX(1,(i*2)-1) = "["+a_points.Var1(i) + ", ";
        airfoilXYZX(1,i*2) = a_points.Var2(i)+"], ";
       
    else
        
        airfoilXYZX(1,(i*2)-1) = "["+a_points.Var1(i) + ", ";
        airfoilXYZX(1,i*2) = a_points.Var2(i)+"]";
        
    end
end

%% save it in .txt file

fid = fopen(fileNameNDirectTXT,'a');


  if fid > 0
     fprintf(fid,"airfoil = [");
     fprintf(fid,'%s',airfoilXYZX');
     fprintf(fid,"];");
     fclose(fid);
  end
  
%% Adjust Parameters

% adjust y Translation, so that airfoil Leading Edge stays on x-Axis
if InputParam.stayOnxAxis == 1;
    
    InputParam.yTransl = InputParam.xScale * sin(InputParam.zRotate);
    
end
  
%% add 3D parameters to end of airfoil file
fid = fopen(fileNameNDirectTXT,'a');

  if fid > 0
     fprintf(fid,'\n \n');
     fprintf(fid,"rotate([%d,%d,%d]) \n",InputParam.xRotate, InputParam.yRotate, InputParam.zRotate);
     fprintf(fid,"linear_extrude(height=%i) \n",InputParam.ExtrH);
     fprintf(fid,"scale([%i, %i]) \n",InputParam.xScale, InputParam.yScale);
     fprintf(fid,"translate([%d,%d, %d]) \n",InputParam.xTransl, InputParam.yTransl, InputParam.zTransl);
     fprintf(fid,"scale (%d) \n",InputParam.Scale);
     fprintf(fid,"polygon(points=airfoil);");
     fclose(fid);
  end
  
%% change format from .txt to .scad and rename it into initial name

system(['mv ',fileNameNDirectTXT,' ',fileNameNDirectSCAD]);

%% Run airfoil in openscad and export it into .stl -> only 
if MAC == 1
    openscad = '/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD';
    system([openscad,' -o ',fileNameNDirectSTL,' ',fileNameNDirectSCAD]);
else 
    % if you want to open the openSCAD gui
    %system(['openscad ',directions.toolPath,'/',GeoParam.airfoilName,'.scad']);

    % if you just want to rewrite .scad int .stl
    system(['openscad -o ',fileNameNDirectSTL,' ',fileNameNDirectSCAD]);
end

%and then delete the old file
system(['rm -f ',fileNameNDirectSCAD]);

%% Find and replace OpenSCAD_Model with ascii

searchnreplace(InputParam.airfoilName,OutputDirectory,".stl","OpenSCAD_Model",InputParam.airfoilName);

%% Clear variables

clear a_points airfoil fid i n
end
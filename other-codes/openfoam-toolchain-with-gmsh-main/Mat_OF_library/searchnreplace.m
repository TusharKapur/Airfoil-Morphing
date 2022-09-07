function searchnreplace (fileName, fileDirectory, fileType, oldString, newString)

    % searchnreplace (fileName, fileDirectory, fileType, oldString, newString)
    % search and replace an 'oldString' with a 'newString' (only Strings!) from the
    % specified 'fileName' within the 'fileDirectory' with the 'dataType' ending 
    % 
    % e.g. searchnreplace ('airfoil', '/home/Matlab_Toolchain', '.txt',
    % 'OpenSCAD', 'ASCII')
    %
    % file ending is added within function through fileType, Path without
    % '/ in the ending'
    
    %% prepare fileName and read in data
    % add datatype ending to fileName 
    fileName_full = append(fileDirectory,'/',fileName,fileType);
    
    % open for extracting data
    fid = fopen(fileName_full,'r');
    
    if fid == -1
        error ('Cannot open file %s', fileName_full);
    end
    
    fileData = fileread(fileName_full);
    fclose(fid);
    
    %% replace oldString with newString
    % replace oldString with newString
    findrep = strrep(fileData, oldString, newString);
    
    %% open again with writing permission, discard, write and close
    fid = fopen(fileName_full, 'w');
    % check whether file can be opened
    if fid == -1
        error ('Cannot open file for writing: %s', fileName_full);
    end 
 
    fprintf(fid,'%s',findrep);
    fclose(fid);

%% clear function variables
clear fid fileData fileName_full

end
function [vertices, faces, lines] = ReadMeshFile ( fileName)

    [filepath, name, ext] = fileparts(fileName);

    vertices = [];
    
    faces = [];
    
    lines = [];

    
if strcmp (ext, '.obj')
    
    fid = fopen (fileName);
    
    lstr = fgetl (fid);
    
    while ischar (lstr)
        
        if (size (lstr, 2) > 0) % skip empty lines
            
            strArr = split (lstr);
            
            leadingChar = strArr{1};
            
            if strcmp(leadingChar, 'v')
                
                val = str2double(strArr(2:end))';
                
                vertices = [vertices; val];
                
            end
            
            if strcmp(leadingChar, 'f')
                
                val = str2double(strArr(2:end))';
                
                faces = [faces; val];
                
            end
            
            if strcmp(leadingChar, 'l')
                
                val = str2double(strArr(2:end))';
                
                lines = [lines; val];
                
            end
            
        end
        
        lstr = fgetl (fid);
        
    end
    
    fclose (fid);
    
end



if strcmp (ext, '.ply')
    
    fid = fopen (fileName);
    
    lstr = fgetl (fid);
    
    while ischar (lstr)
        
        if (size (lstr, 2) > 0) % skip empty lines
            
            val = str2num(lstr); % row vector
            
            if length(val) == 3 || length(val) == 6
                
                vertices = [ vertices; val(1), val(2), val(3) ];
                
            end
            
            if length(val) == 4 && val(1) == 3
                
                faces = [ faces; val(2)+1, val(3)+1, val(4)+1 ];
                
            end 
             
        end
        
        lstr = fgetl (fid);
        
    end
    
    fclose (fid);
    
end
    
end


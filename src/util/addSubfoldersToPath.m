function setPath(path_dir)

if nargin < 1
  here = mfilename('fullpath');
  [~, path_here] = strtok(fliplr(here), '/');
  path_here = fliplr(path_here);
  path_dir = path_here;
end

dlist = dir(path_dir);

for j = 1:length(dlist)
  
 if dlist(j).isdir && ~strcmp(dlist(j).name, '.') && ~strcmp(dlist(j).name, '..')
   addpath([path_dir '/' dlist(j).name]) 
 end
 
end

end
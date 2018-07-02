function str = assignFields(str, varargin)

for j = 1:2:length(varargin)
  str = setfield(str, varargin{j}, varargin{j+1});
end
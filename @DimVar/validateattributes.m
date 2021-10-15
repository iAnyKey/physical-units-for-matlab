function validateattributes(a,classes,varargin)
% Categorize DimVar as numeric for purposes of validating attributes.

% % import functions in case if repository has been includen in a package.
% % if not - `import .*` does nothing 
% eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

% handle the case of being callsed out of package
sClassName = 'DimVar';
sPkgName = strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.');
if ~isempty(sPkgName)
    sClassName = [sPkgName '.' sClassName];
end

if any(strcmp('numeric',classes)) || any(strcmp(sClassName,classes))
    % DimVar is allowed. Evaluate against value.
    
    classes = strrep(classes,sClassName,'double');
    % Does not need a 'unique' call since validateattributes is okay
    % with doubling up the classes input (e.g. {'double' 'double'}.
    
    validateattributes(a.value,classes,varargin{:})
else
    builtin('validateattributes',a,classes,varargin{:})
    % Use built-in function to throw appropriate error for
    % non-DimVar input.
end
end
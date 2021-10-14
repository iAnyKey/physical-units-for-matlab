function [tf,ME] = iscompatible(varargin)
% Returns true if all inputs are DimVars with the same units (compatible
% executes successfully) and otherwise returns false.
% 
%   See also u, compatible.

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

try 
    compatible(varargin{:});
    tf = true;
catch ME
    % Capture all-double case.
    if all(cellfun('isclass',varargin,'double'))
        tf = true;
        return
    end

    if strcmp(ME.identifier,'DimVar:incompatibleUnits')
        tf = false;
    else
        rethrow(ME)
    end

end
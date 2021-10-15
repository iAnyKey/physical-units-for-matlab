function v = round(v,varargin)

% % import functions in case if repository has been includen in a package.
% % if not - `import .*` does nothing 
% eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

warning('DimVar:round','Using round with DimVars may yield unexpected results.')

dispVal = displayparser(v);
delta = round(dispVal,varargin{:})./dispVal;

v.value = v.value.*delta;

% v.value = round(v.value,varargin{:});
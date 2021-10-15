function vOut = cumtrapz(v1,v2,varargin)

% % import functions in case if repository has been includen in a package.
% % if not - `import .*` does nothing 
% eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

% handle the case of being callsed out of package
sClassName = 'DimVar';
sPkgName = strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.');
if ~isempty(sPkgName)
    sClassName = [sPkgName '.' sClassName];
end

if ~isa(v2,sClassName) % v1 is the DimVar.
    vOut = v1;
    vOut.value = cumtrapz(v1.value,v2,varargin{:});


elseif ~isa(v1,sClassName) % v2 is the DimVar.
    vOut = v2;
    vOut.value = cumtrapz(v1,v2.value,varargin{:});

else % BOTH v1 and v2 are DimVars.
    vOut = v1;
    vOut.value = cumtrapz(v1.value,v2.value,varargin{:});
    vOut.exponents = v1.exponents + v2.exponents;
    
    vOut = clearcanceledunits(vOut);
end
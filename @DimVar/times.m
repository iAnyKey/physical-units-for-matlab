function v1 = times(v1,v2)

% % import functions in case if repository has been includen in a package.
% % if not - `import .*` does nothing 
% eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

% handle the case of being callsed out of package
sClassName = 'DimVar';
sPkgName = strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.');
if ~isempty(sPkgName)
    sClassName = [sPkgName '.' sClassName];
end

if ~isa(v2,sClassName) % v1 is the only DimVar.
    v1.value = v1.value .* v2;
elseif ~isa(v1,sClassName) % v2 is the only DimVar.
    v2.value = v1 .* v2.value;
    v1 = v2;
else % BOTH v1 and v2 are DimVars.
    v1.value = v1.value .* v2.value;
    v1.exponents = v1.exponents + v2.exponents;
    v1 = clearcanceledunits(v1);
end
function v = mpower(v,y)

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

if isa(y,'DimVar')
    error('For Z = X^y, y may not be a DimVar.');
end
v.value = v.value^y;
v.exponents = y*v.exponents;

v = clearcanceledunits(v);
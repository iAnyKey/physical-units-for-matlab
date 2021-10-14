function v = sqrt(v)

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

v.value = sqrt(v.value);
v.exponents = v.exponents/2;
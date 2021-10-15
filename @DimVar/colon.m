function a = colon(a,b,c)

% % import functions in case if repository has been includen in a package.
% % if not - `import .*` does nothing 
% eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

if nargin == 2
    error('DimVar:colon:incrementRequired',...
        'DimVar vector creation with colon is only defined for three inputs.')
end

compatible(a,b,c)

a.value = a.value:b.value:c.value;
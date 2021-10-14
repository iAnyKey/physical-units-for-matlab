function disp(v)

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

[val,appendStr] = displayparser(v); %#ok<ASGLU>

t = evalc('disp(val)');
if isempty(t)
    t = sprintf('\t[]');
end

switch get(0,'FormatSpacing')
    case 'compact'
        fprintf('%s\t%s\n',deblank(t),appendStr);
    case 'loose'
        fprintf('%s\t%s\n\n',deblank(t),appendStr);
end
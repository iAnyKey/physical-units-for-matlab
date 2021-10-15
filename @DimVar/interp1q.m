function yi = interp1q(x,y,xi)
% See also interp1q.

% % import functions in case if repository has been includen in a package.
% % if not - `import .*` does nothing 
% if isempty(which('compatible')) || isempty(which('unitsOf'))
%     eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));
% end

compatible(x,xi);
outUnits = unitsOf(y);
yi = outUnits*interp1q(double(x),y/outUnits,double(xi));
function v = scale(v, sf)
%  SCALE(input, scaleFactor) performs a dynamic scaling operation by the
%  scale factor on the input variable.
% 
%   To learn about dynamic scaling:
%   https://www.google.com/search?q=dynamic+scaling+aircraft
% 
%   See also u.

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

a = 0*v.exponents;
a(1:7) = [
    1   % Length
    3   % Mass
    1/2 % Time
    9/4 % Current (driven by Ohm scaling with length/area, K/K�)
    0   % Temp
    3   % Amount
    7/2 % cd (driven by matching scaling of watts/steridian)
    ];


v.value = v.value .* sf.^sum(v.exponents.*a);


function d = u2duration(v)
% u2duration  Convert DimVar time to duration data type.
%   D = u2duration(V) converts time value contained in DimVar V to a
%   duration.
% 
%   u2duration inspects the name of the time dimension and recognizes most
%   common string expressions of seconds, minutes, hours, or days. If you
%   have redefined the time dimension name to something that could be
%   ambiguous, use caution.
% 
%   Note that units uses the Julian year (which defines the light-year),
%   rather than the Gregorian year that duration uses.
% 
%   See also duration, seconds, u, u2num.

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

% handle the case of being callsed out of package
sClassName = 'DimVar';
sPkgName = strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.');
if ~isempty(sPkgName)
    sClassName = [sPkgName '.' sClassName];
end

if isa(v/u.s,sClassName)
    error('A pure time DimVar (with exponent of one) is required.')
else
    [~,unitStr] = displayparser(v);
    switch unitStr
        case {'s' 'ss' 'sec' 'secs' 'second' 'seconds'}
            d = seconds(v/u.s);
        case {'m' 'mm' 'min' 'mins' 'minute' 'minutes'}
            d = minutes(v/u.min);
        case {'h' 'hh' 'hr' 'hrs' 'hour' 'hours'}
            d = hours(v/u.hr);
        case {'d' 'dd' 'day' 'days'}
            d = days(v/u.day);
        case {'y' 'yr' 'yrs' 'year' 'years'}
            d = years(v/u.year_Gregorian);
        otherwise
            d = duration(0,0,v/u.s);
    end        
end
% A test script for overloaded DimVar subsasgn method.
%   <a href="matlab:runtests('testScript_subsasgn')">run tests</a>

% Set up by clearing class
clear all

%% assign dim to dim
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1d = u.m*(1:4);
v1d([2,6]) = 9*u.m; 
assert(isequal(v1d,u.m*[1     9     3     4     0     9]));

%% new variable
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
clear v
v(3) = u.furlong; 
assert(isequal(v,u.furlong*[0 0 1]))

%% assign normal to dim
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1d = u.m*(1:4); %#ok<NASGU>
shoulderror('DimVar:subsasgn:invalidAssignment','v1d(3) = 8');

%% remove element, using a variable (errors with doubles)
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1d = u.m*(1:4);
a = [];
shoulderror('v1d(3) = a;');

%% remove element, special [] syntax
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1d = u.m*(1:4);
v1d(1) = []; 
assert(isequal(v1d,u.m*[2 3 4]));

%% assign NaN to DimVar allowed case
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1d = u.m*(1:4);
v1d(3) = NaN; 
assert(isequaln(v1d,u.m*[1 2 NaN 4]));

%% assign NaNdim to DimVar
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1d = u.m*(1:4);
v1d(3) = u.m*NaN; 
assert(isequaln(v1d,u.m*[1 2 NaN 4]));

%% assign to empty DimVar
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
ve = []*u.m;
ve(3) = u.m;
assert(isequal(ve,u.m*[0 0 1]))

%% don't allow shrinking array without consistent units
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
ve = (1:4)*u.m;
shoulderror('DimVar:incompatibleUnits','ve(3) = []*u.kg;');

%% don't allow changing units of all-nan array
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v = nan(3,5)*u.m;
shoulderror('DimVar:incompatibleUnits','v(5) = u.kg');

%% assign to empty, wrong unit
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
ve = []*u.m;
shoulderror('DimVar:incompatibleUnits','ve(3) = u.kg;');

%% expected to fail: assign dim to empty norm (doesn't call overloaded method)
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
ve = [];
ve(3) = u.m;
shoulderror("assert(isequal(ve,u.m*[0 0 1]))");

%% assign dim to empty normal, subsasgn call
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
ve = [];
S.type = '()';
S.subs = {3};
ve = subsasgn(ve,S,9*u.m);
assert(isequal(ve,u.m*[0 0 9]))

%% expected to fail: assign dim to normal (doesn't call overloaded method)
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1 = 1:4;
shoulderror("v1([2,6]) = 9*u.lb;")%FIXME: something is wrong over here - no error has been thrown and this statement works well for me

%% assign dim to normal, subsasgn call
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1 = 1:4;
S.type = '()';
S.subs = {[2 6]};
shoulderror("subsasgn(v1,S,9*u.m)");

%% assign dim to empty normal
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v = [];
v = subsasgn(v,substruct('()',{3}),u.kg);
assert(isequal(v,u.kg*[0 0 1]))

%% expected to fail: assign DimVar to NaN (doesn't call overloaded method)
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
v1n = NaN(1,4);
v1n(3) = u.m;

%% -- local copy of the shoulderror --
% unfortunately a functioncall clears imported namespaces, so we need to
% import everything again
function varargout = shoulderror(varargin)
% shoulderror  Throws an error if input executes successfully. Useful for test
% scripts.
% 
%   Inputs and outputs are the same as for feval or, if only one char or string
%   input is provided, evalin('caller',...).
% 
%   shoulderror returns the MException object if the error is thrown, though use
%   caution regarding the input's expectations regarding number of outputs.
% 
%   shoulderror(ID,...), where the first argument is any character vector or
%   string scalar containing a colon character, also verifies that the
%   MException identifier of the error thrown matches ID.
% 
%   Examples:
%     shoulderror(@cos,"pi") % Executes without error, returning MException.
%     shoulderror('MATLAB:UndefinedFunction',@cos,"pi") % Checks error ID.
%     shoulderror('mrdivide',1,3) % Errors.
%     a = shoulderror('1 + 1') % Errors; a is not assigned.
%     
%     Be careful with number of expected outputs. These yield different results:
%               shoulderror('MATLAB:deal:narginNargoutMismatch',@deal,1,2) 
%       [~,~] = shoulderror('MATLAB:deal:narginNargoutMismatch',@deal,1,2)
% 
%   See also shouldalert, shouldwarn, runtests, feval, MException.

% Copyright 2018 Sky Sartorius. All rights reserved.
% Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715 

eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 

ID = varargin{1};
testID = false;
if nargin > 1 ...
        && (ischar(ID) || isStringScalar(ID)) ...
        && ~isempty(regexp(ID,'\w+:\w+','once'))
    varargin = varargin(2:end);
    testID = true;
end
nArgs = numel(varargin);

errored = false;
try 
    if nArgs == 1 && (ischar(varargin{1}) || isstring(varargin{1}))
        [varargout{1:nargout}] = evalin('caller',varargin{1}); %#ok<*NASGU>
    else
        [varargout{1:nargout}] = feval(varargin{:});
    end
catch ME
    % Good, it should have thrown an error.
    varargout = [{ME} cell(1,nargout-1)];
    errored = true;
    if testID && ~strcmp(ME.identifier,ID)
        error('AlertChecking:ShouldError:IncorrectIdentifier',...
            'This should have thrown an error with identifier:\n%s',ID);
    end
end

if ~errored
    error('AlertChecking:ShouldError','This should have thrown an error.')
end

end
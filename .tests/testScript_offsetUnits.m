% A test script for OffsetDimVar.
%   <a href="matlab:runtests('testScript_offsetUnits')">run tests</a>

% Set up by clearing class
clear all

%%   Set a value with multiplication (before or after scalar): 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
assert(20*u.degC==293.15*u.K)
assert(u.degC*40== 313.15*u.K)

%% Set with str2u
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
assert(str2u('20 degC')==293.15*u.K)
assert(str2u('20 °C')==293.15*u.K)

%% Special str2u case
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
assert(isa(str2u('degC'),'OffsetDimVar'))

%%       Convert units with division:
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
assert(abs(200*u.K/u.degC - (-73.15)) < sqrt(eps))
assert(isequal(20*u.degC/u.K,293.15))
assert(abs(u.degC/u.degF - 1.8) < sqrt(eps))

%% Should errors
%     Pretty much all other use cases should be avoided and mostly throw an
%     error, but not always (usually in cases where it will be interpreted at 1
%     degC = 274.15 K, e.g.), so be careful.
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
shoulderror('DimVar:incompatibleUnits','u.kg*u.degC');
shoulderror('DimVar:incompatibleUnits','u.degC*u.kg');
shoulderror('OffsetDimVar:incompatibleUnits','u.degC*u.degF');
shoulderror('OffsetDimVar:undefined','u.degC/u.K');
shoulderror('OffsetDimVar:undefined',"str2u('20 degC/s')");
shoulderror('u.degC + 5*u.K');
shoulderror('DimVar:incompatibleUnits','unitconversionfactor(u.K,u.degC)');

%       anything else with unitconversionfactor (since it won't be a "factor"),
%       e.g. unitconversionfactor('degC','K').

shoulderror('DimVar:incompatibleUnits',"unitconversionfactor('degC','K')");

%% this syntax is bad practice, but it works
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));% import functions in case if repository has been includen in a package. if not - `import .*` does nothing 
20*u.degC + 20*u.degF;

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

classdef DimVar
% See also u, myUnits, units.

% Copyright (c) 2012-2017, Sky Sartorius.
properties (Access = protected)
    names
    exponents
    value = 1
    dispUnits
end

methods
    %% Core functions.
    % Constructor
    function v = DimVar(dimensionNames,dimensionToCreate,dispUnits)
        % See also u.
        v.names = dimensionNames;
        v.exponents = zeros(size(v.names));
        dimensionIndex = strcmp(v.names,dimensionToCreate);
        v.exponents(dimensionIndex) = 1;
        if nargin == 3
            v.dispUnits = dispUnits;
        end
    end
    
    function v = clearcanceledunits(v)
        % If all DimVar unit exponents are zero, return normal (double)
        % variable. Exponent tolerance is to fifth decimal.
        
        if ~any(round(1e5*v.exponents)) % Seems to be faster than round(x,5).
            v = v.value;
        end
    end
    
    function compatible(v,varargin)
        % compatible(v1, v2, ...) throws an error unless all inputs are DimVar
        % with the same units.
        %
        %   If throwing an error is not desired, use iscompatible.
        %
        %   See also u, iscompatible.
        
        if ~isa(v,'DimVar')
            
            ME = MException('DimVar:incompatibleUnits',...
                'Incompatible units. All inputs must be DimVar.');
            throwAsCaller(ME);
            
        end
        
        vExpos = v.exponents;
        
        for i = 1:numel(varargin)
            
            if ~isa(varargin{i},'DimVar') || ~isequal(vExpos,varargin{i}.exponents)
                
                ME = MException('DimVar:incompatibleUnits',...
                    ['Incompatible units. Cannot perform operation on '...
                    'variables with different units.']);
                throwAsCaller(ME);
                
            end
        end
    end
    
    %% Concatenation.
    function v = cat(dim,v,varargin)
        if ~isa(v,'DimVar')
            ME = MException('DimVar:incompatibleUnits',...
                'Incompatible units. All inputs must be DimVar.');
            throwAsCaller(ME);
        end
        
        vExpos = v.exponents;
        
        for i = 1:numel(varargin)
            vi = varargin{i};
            
            if ~isa(vi,'DimVar') || ~isequal(vExpos,vi.exponents)
                ME = MException('DimVar:incompatibleUnits',...
                    ['Incompatible units. Cannot perform operation on '...
                    'variables with different units.']);
                throwAsCaller(ME);
            end
            
            v.value = cat(dim,v.value,vi.value);
        end
        
        % Functionality of compatible method is integrated for sake of speed.
    end
    function vOut = horzcat(varargin)
        vOut = cat(2,varargin{:});
    end
    function vOut = vertcat(varargin)
        vOut = cat(1,varargin{:});
    end
    
    %% Validation functions (mustBe__).
    function mustBeGreaterThan(v1,v2)
        compatible(v1,v2);
        mustBeGreaterThan(v1.value,v2.value);
    end
    function mustBeGreaterThanOrEqual(v1,v2)
        compatible(v1,v2);
        mustBeGreaterThanOrEqual(v1.value,v2.value);
    end
    function mustBeLessThan(v1,v2)
        compatible(v1,v2);
        mustBeLessThan(v1.value,v2.value);
    end
    function mustBeLessThanOrEqual(v1,v2)
        compatible(v1,v2);
        mustBeLessThanOrEqual(v1.value,v2.value);
    end
    function mustBeNegative(v);     mustBeNegative(v.value);    end
    function mustBeNonnegative(v);  mustBeNonnegative(v.value); end
    function mustBeNonpositive(v);  mustBeNonpositive(v.value); end
    function mustBeNonzero(v);      mustBeNonzero(v.value);     end
    function mustBePositive(v);     mustBePositive(v.value);    end
    
    %% is__ functions.
    function result = isempty(v);   result = isempty(v.value);      end
    function result = isfinite(v);  result = isfinite(v.value);     end
    function result = isinf(v);     result = isinf(v.value);        end
    function result = isnan(v);     result = isnan(v.value);        end
    function result = isnumeric(v); result = isnumeric(v.value);    end
    function result = isreal(v);    result = isreal(v.value);       end
    
    %% Logical operators (>, <, ==, ~, etc.).
    function result = eq(v1,v2)
        compatible(v1,v2);
        result = v1.value == v2.value;
    end
    function result = ge(v1,v2)
        compatible(v1,v2);
        result = v1.value >= v2.value;
    end
    function result = gt(v1,v2)
        compatible(v1,v2);
        result = v1.value > v2.value;
    end
    function result = le(v1,v2)
        compatible(v1,v2);
        result = v1.value <= v2.value;
    end
    function result = lt(v1,v2)
        compatible(v1,v2);
        result = v1.value < v2.value;
    end
    function result = ne(v1,v2)
        compatible(v1,v2);
        result = v1.value ~= v2.value;
    end
    function result = not(v)
        result = ~v.value;
    end
    
    %% Class conversions.
    function result = double(v)
        % DimVar.double(V)  Returns value property (not diplayingvalue) of V.
        %
        % See also u2num, displayingvalue, displayparser.
        result = v.value;
    end
    function result = logical(v)
        result = logical(v.value);
    end
    function s = string(v)
        s = string(cellfun(@num2str,num2cell(v),'UniformOutput',false));
    end
    
    %% Simple functions that return DimVar.
    function v = abs(v);            v.value = abs(v.value);                 end
    function v = circshift(v,varargin)
        v.value = circshift(v.value,varargin{:});
    end
    function v = conj(v);           v.value = conj(v.value);                end
    function v = ctranspose(v);     v.value = v.value';                     end
    function v = cumsum(v,varargin);v.value = cumsum(v.value,varargin{:});  end
    function v = diag(v,varargin);  v.value = diag(v.value,varargin{:});    end
    function v = diff(v,varargin);  v.value = diff(v.value,varargin{:});    end
    function v = full(v);           v.value = full(v.value);                end
    function v = imag(v);           v.value = imag(v.value);                end
    function v = mean(v,varargin);  v.value = mean(v.value,varargin{:});    end
    function v = median(v,varargin);v.value = median(v.value,varargin{:});  end
    function v = norm(v,varargin);  v.value = norm(v.value,varargin{:});    end
    function v = permute(v,varargin);v.value = permute(v.value,varargin{:});end
    function v = real(v);           v.value = real(v.value);                end
    function v = reshape(v,varargin);v.value = reshape(v.value,varargin{:});end
    function [v,I] = sort(v,varargin)
        [v.value, I] = sort(v.value,varargin{:});
    end
    function v = std(v,varargin);   v.value = std(v.value,varargin{:});     end
    function v = subsref(v,varargin);v.value = subsref(v.value,varargin{:});end
    function v = sum(v,varargin);   v.value = sum(v.value,varargin{:});     end
    function v = trace(v);          v.value = trace(v.value);               end
    function v = transpose(v);      v.value = v.value.';                    end
    function v = uminus(v);         v.value = -v.value;                     end
    function v = uplus(v);                                                  end
    
    %% Functions that require compatibility check.
    function out = atan2(v1,v2)
        compatible(v1,v2);
        out = atan2(v1.value, v2.value);
    end
    function v1 = hypot(v1,v2)
        compatible(v1,v2);
        v1.value = hypot(v1.value,v2.value);
    end
    function v1 = minus(v1,v2)
        compatible(v1,v2);
        v1.value = v1.value - v2.value;
    end
    function v1 = plus(v1,v2)
        compatible(v1,v2);
        v1.value = v1.value + v2.value;
    end
    
    %% Other simple functions.
    function n = length(v);     n = length(v.value);    end
    function n = ndims(v);      n = ndims(v.value);     end
    function n = nnz(v);        n = nnz(v.value);       end
    function n = numel(v);      n = numel(v.value);     end
    function out = sign(v);     out = sign(v.value);    end
    function varargout = size(x,varargin)
        [varargout{1:nargout}] = size(x.value,varargin{:});
    end

end
methods (Static)
    [listString,list] = unitslist(varargin)
    [cTo,cInverse] = strconv(sFrom,sTo,varargin)

    function displayunitsstruct(u)
        % Display the struct created by calling units.
        %
        %   See also units.
        if nargin < 1
            u = units;
        end
        uNames = fieldnames(u);
        nFields = length(uNames);
        lastUnits = 'wR!2V!17r@OenDg6wH&t';
        for nf = 1:nFields
            name = uNames{nf};
            [~,currentUnits] = unitsOf(u.(name));
            if(~strcmp(currentUnits,lastUnits))
                fprintf('\n');
                lastUnits = currentUnits;
            end

            fprintf('%16s: %-s\n',name,num2str(u.(name)));
        end
    end
end
end
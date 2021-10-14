function endIndex = end(vIn,indexToFindEnd,numberOfIndexes)

% import functions in case if repository has been includen in a package.
% if not - `import .*` does nothing 
eval(sprintf('import %s.*', strjoin(regexp(mfilename('fullpath'), '(?<=+)\w*', 'match'), '.')));

sizeVector = size(vIn.value);

if(numberOfIndexes == 1)
    endIndex = max(sizeVector);
else
    endIndex = sizeVector(indexToFindEnd);
end

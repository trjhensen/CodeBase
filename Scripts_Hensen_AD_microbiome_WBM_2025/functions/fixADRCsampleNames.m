function fixedNames = fixADRCsampleNames(namesToFix)
% Problem statement: In the gut microbiome samples names, trailing 
% zeros have been removed from some of the samples in the gut 
% microbiome data. The sample IDs should have eight digits. However, 
% a sample ID such as 42336000 is pruned to 42336 in the microbiome 
% samples. To repair these samples, we need to find sample IDs that 
% have less than 8 digits and append trailing zeros until eight digits 
% are reached again. 

% input: Cell array of names to fix

namesToFix = cellstr(namesToFix);

% Check if a leading X is present
rmLeadingX = any(~cellfun(@isempty,regexp(namesToFix,'^X','match')));

if rmLeadingX == true
    % Remove leading X if present in the sample names
    namesToFix = regexprep(namesToFix,'^X','');
end

% Test if is indeed a maximum of eight characters
if max(cellfun(@length, namesToFix))>8
    error('The sample names must have a maximum of eight characters')
end

% the following line find for each cell entry the number of missing
% trailing zeros and then adds the exact number of missing trailing zeros
% so that the character length of the sample names is exactly eight.
fixedNames = cellfun(@(x) [x repmat('0', 1, 8 - length(x))], namesToFix, 'UniformOutput', false);

% The samples 5827089 and 5946764, however, need to have 7 characters to map
% onto the patient metadata, so lets remove the last trailing zero if
% present.
fixedNames(matches(fixedNames,{'58270890','59467640'})) = {'5827089','5946764'};

% Add back the leading X if present
if rmLeadingX == true
    % Remove leading X if present in the sample names
    fixedNames = append('X',fixedNames);
end

end
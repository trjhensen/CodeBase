function metadataProcessed = rmAdrcOverlapAndMissing(metadata)

% There are two individuals with more than one sample in the metadata. 
% Remove them from the metadata. There are overlapping samples between year
% 1 and year 3.

% Find unique individuals
participants1 = metadata.NACCID(metadata.Visit_year==1);
participants2 = metadata.NACCID(metadata.Visit_year==2);
participants3 = metadata.NACCID(metadata.Visit_year==3);

% Overlap in year one and two
year1year2Overlap = intersect(participants1,participants2);
metadata(matches(metadata.NACCID,year1year2Overlap),:) = [];

% Overlap between year 2 and year 3.
year2year3Overlap = intersect(participants2,participants3);
metadata(matches(metadata.NACCID,year2year3Overlap),:) = [];

% Remove participants with no dementia information
metadata(ismissing(string(metadata.NACCUDSD)),:) = [];

metadataProcessed = metadata;
end
function refs = extractRefs(years, rmvRefs, keepJour)
% extractRefs.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), May 2022
% 
% this function grabs, processses, and stores information stored in ice core
% abstracts for each year, so that information can be analyzed + related to 
% text patterns over time

% DEPENDENCIES
% 1. Dependencies/Data/refsAll.xlsx
% 2. Dependencies/Data/citations.mat
% 3. Dependencies/Data/affils.mat

% INPUTS
% years -- vector of years to extract; max range is 1969--2021; 
if nanmin(years) < 1969
    error('Input year range for ''years'' must be >=1969'); 
elseif nanmax(years) > 2021
    error('Input year range for ''years'' must be <=2021'); 
end

% rmvRefs -- name of a file that contains a structure named "list,
%   with the following subfields:  abstracts, titles, journals, keywords, venues
if nargin < 2 
    rmvRefs = []; 
end

% keepJour -- name of journals to retain (for list, see file
%   extractRef_exploreJournal.m); defaults to retaining all
if nargin < 3
    keepJour = []; 
end

% Output
% refs -- a structure containing subfields with ice core study author names, 
%   ("fullNames"), "abstract" text, study "titles", study "journals",
%   any associated study "keywords", study "doi", study "citations" since
%   Sept 2022, study publication date ("pubdate"), and whether each study is
%   "openacess" or not. Each variable is a cell vector of dimension
%   length(years) x 1. 

% % Example inputs:
% years = [2001:2021]'; 
% rmvRefs = 'rmvRefs.mat'; 
% keepJour = ["Annals of Glaciology";
%             "Atmospheric Chemistry & Physics";
%             "Atmospheric Environment";
%             "Atmospheric Measurement Techniques";
%             "Climate Dynamics";
%             "Climate of the Past";
%             "Earth and Planetary Science Letters";
%             "Frontiers in Earth Science";
%             "Geochimica et Cosmochimica Acta";
%             "Geology"]; 
    
%% load in ice core ADS-derived .ris files appended in refsAll

% load raw abstract data
[~, text] = xlsread('refsAll.xlsx','combined');  text = string(text);  I = text == "";  text(I) = []; clearvars I; 

% identify indices
tf = startsWith(text,'TY  - ');
I = [find(tf), find([tf(2:end); true])]; % column1 = start index, column2 = end index

% filter now by conf, jour, or both
Ij = contains(text(I(:,1)),'JOUR'); 
I = I(Ij,:); % refine I to only give indices of journals

% identify year for each abstract
absYear = nan(size(I,1),1); 
for i = 1:length(absYear)
    currText = text(I(i,1):I(i,2)); 
    currYear = currText(startsWith(currText,'Y1  - ')); 
    absYear(i) = str2double(extractBetween(currYear,'Y1  - ','/'));  
end

% number of abstracts each year
numAbs = nan(length(years),1); 
for i = 1:length(years) 
    numAbs(i,1) = sum(absYear == years(i)); 
end

% preallocate things to record
allInfo = cell(length(years),1); 
fullNames = cell(length(years),1); 
affiliations = cell(length(years),1); 
abstracts = cell(length(years),1); 
titles = cell(length(years),1); 
journals = cell(length(years),1); 
keywords = cell(length(years),1); 
venues = cell(length(years),1); 
doi = cell(length(years),1); 
% grab data
h = waitbar(0,'Sorting abstract data ... please wait!'); 
for i = 1:length(years)
    It = I(absYear == years(i),:); 
    for j = 1:size(It,1)
        absCurr = text(It(j,1):It(j,2)); 
        allInfo{i}{j,1} = absCurr; 
        % venue
        Ir = startsWith(absCurr,'TY  - ');
        rowCurr = absCurr(Ir);
        venues{i}{j,1} = extractAfter(rowCurr,'TY  - ');
        % doi
        Ir = startsWith(absCurr,'DO  - ');
        rowCurr = absCurr(Ir);
        doi{i}{j,1} = extractAfter(rowCurr,'DO  - ');
        % journals
        Ir = startsWith(absCurr,'JO  - ');
        rowCurr = absCurr(Ir);
        journals{i}{j,1} = extractAfter(rowCurr,'JO  - ');
        % titles
        Ir = startsWith(absCurr,'TI  - ');
        rowCurr = absCurr(Ir);
        titles{i}{j,1} = extractAfter(rowCurr,'TI  - ');
        % abstract
        Ir = startsWith(absCurr,'N2  - ');
        rowCurr = absCurr(Ir);
        abstracts{i}{j,1} = extractAfter(rowCurr,'N2  - ');% rowCurr(12:(strfind(rowCurr,"},")-1));         
        % keywords
        Ir = startsWith(absCurr,'KW  - ');
        rowCurr = absCurr(Ir);
        keywords{i}{j,1} = extractAfter(rowCurr,'KW  - ');% rowCurr(12:(strfind(rowCurr,"},")-1));         
        % affiliations
        Ir = startsWith(absCurr,'AD  - ');
        rowCurr = absCurr(Ir);
        affiliations{i}{j,1} = extractBetween(rowCurr,'(',')');% rowCurr(12:(strfind(rowCurr,"},")-1));         
        % first names
        Ir = startsWith(absCurr,'AU  - ');
        rowCurr = absCurr(Ir);
        rowCurr = extractAfter(rowCurr,'AU  - '); 
        fullNames{i}{j,1} = strcat(extractAfter(rowCurr,', ')," ",extractBefore(rowCurr,', ')); 
    end
    waitbar(i / length(years))
end
close(h)

%% Remove non-unique abstracts 
% (i.e., remove overlapping pairs from tabs "ice core" and "ice cores", 
% which are both conveniently stored in tab "combined" in refsAll.xlsx)

for i = 1:length(years)
    [~, Ic] = unique(string(titles{i})); Ic = sort(Ic,'ascend'); 
    % retain only the unique abstracts
    allInfo{i} = allInfo{i}(Ic); 
    fullNames{i} = fullNames{i}(Ic); 
    affiliations{i} = affiliations{i}(Ic); 
    abstracts{i} = abstracts{i}(Ic); 
    titles{i} = titles{i}(Ic); 
    journals{i} = journals{i}(Ic); 
    keywords{i} = keywords{i}(Ic); 
    venues{i} = venues{i}(Ic); 
    doi{i} = doi{i}(Ic); 
end

%% Retain only certain journals, if applicable (must come after gender organization)

% remove pesky, unnecessary text descriptors
for i = 1:length(journals)
	for j = 1:size(journals{i},1)
	newStr = extract(journals{i}{j},"Journal of Geophysical Research");
        if ~isempty(newStr)
            journals{i}{j} = newStr; 
        end
        if contains(journals{i}{j}," Supplement Series")
            journals{i}{j} = erase(journals{i}{j}," Supplement Series"); 
        end
        if contains(journals{i}{j}," Supplement")
            journals{i}{j} = erase(journals{i}{j}," Supplement"); 
        end
        if contains(journals{i}{j}," Discussions")
            journals{i}{j} = erase(journals{i}{j}," Discussions"); 
        end
    end
end

% option to retain only certain journals
if ~isempty(keepJour) 
    for i = 1:length(years)
        Ic = false(size(journals{i},1),1);         
        for j = 1:length(keepJour)
            % remove overlapping
            I = strcmp(keepJour(j),string(journals{i})); 
            Ic(I) = true;
        end
        if sum(Ic) > 0
        % remove variables
        allInfo{i}(~Ic) = []; 
        fullNames{i}(~Ic) = []; 
        affiliations{i}(~Ic) = []; 
        abstracts{i}(~Ic) = []; 
        titles{i}(~Ic) = []; 
        journals{i}(~Ic) = []; 
        keywords{i}(~Ic) = []; 
        venues{i}(~Ic) = []; 
        doi{i}(~Ic) = []; 
        end
    end
end

%% Remove abstracts that aren't ice core-related

if ~isempty(rmvRefs)
    if isfile(strcat("Dependencies/Data/",rmvRefs))
        load(rmvRefs)
        % disp(['Removing article:'])
        for i = 1:length(years)
            Ic = false(size(titles{i},1),1); 
            for j = 1:size(titles{i},1)
                % remove overlapping titles
                if sum(strcmp(string(list.titles),titles{i}{j})) > 0
                    % disp(titles{i}{j})
                    Ic(j) = true;
                end
            end
            if sum(Ic) > 0
            % remove variables
            allInfo{i}(Ic) = []; 
            fullNames{i}(Ic) = []; 
            affiliations{i}(Ic) = []; 
            abstracts{i}(Ic) = []; 
            titles{i}(Ic) = []; 
            journals{i}(Ic) = []; 
            keywords{i}(Ic) = []; 
            venues{i}(Ic) = []; 
            doi{i}(Ic) = []; 
            end
        end
    else
    error(['File ',rmvRefs,' does not exist in subdirectory /Dependencies/Data/']); 
    end
end

%% Infill missing author affiliations

load affils.mat % load author affiliations

[~,I] = intersect(Year,years);
Year = Year(I); 
Affiliation = Affiliation(I); 
Title = Title(I); 
Author = Author(I); 
DOI = DOI(I); 
for i = 1:length(Affiliation)
    if ~isempty(Affiliation{i})
        clearvars currDOI; 
        for k = 1:length(doi{i}) % for each year, is gathering all DOI's that exist
            if ~isempty(doi{i}{k})
            currDOI(k,1) = doi{i}{k}; 
            else
            currDOI(k,1) = "NaN"; 
            end
        end
        for j = 1:length(Affiliation{i}) 
            if ~isempty(Affiliation{i}{j}) % sort through each infilled Affiliation
                I = find(strcmp(currDOI,DOI{i}{j})); % match that affiliation to the corresponding DOI
                if ~isempty(I)
                    for k = 1:length(I) % in case DOI is double listed
                        if ~isempty(Affiliation{i}{j}) && isempty(affiliations{i}{I(k)})
                            affiliations{i}{I(k)} = Affiliation{i}{j}; 
                        end
                    end
                end
            end
        end
    end
end
clearvars Year Affiliation Title Author DOI

%% Infill citation information 

load citations.mat % load study citations

pubdate = cell(size(titles)); 
citations = cell(size(titles)); 
openaccess = cell(size(titles)); 
for i = 1:length(titles) % preallocate with NaNs
    if ~isempty(titles{i})
        pubdate{i} = cell(size(titles{i},1),1); 
        citations{i} = nan(size(titles{i},1),1); 
        openaccess{i} = nan(size(titles{i},1),1); 
    end
end

[~,I] = intersect(Year,years);
Year = Year(I); 
pubDate = pubDate(I); 
citationCount = citationCount(I); 
openAccess = openAccess(I); 
DOI = DOI(I); 
for i = 1:length(citationCount)
    if ~isempty(citationCount{i})
        clearvars currDOI; 
        for k = 1:length(doi{i}) % for each year, is gathering all DOI's that exist
            if ~isempty(doi{i}{k})
            currDOI(k,1) = doi{i}{k}; 
            else
            currDOI(k,1) = "NaN"; 
            end
        end
        for j = 1:length(citationCount{i}) 
            if ~isempty(citationCount{i}(j)) && ~isempty(DOI{i}{j})  % sort through each infilled Affiliation
                I = find(strcmp(currDOI,DOI{i}{j})); % match that affiliation to the corresponding DOI
                if ~isempty(I)
                    for k = 1:length(I) % in case DOI is double listed
                        if ~isempty(citationCount{i}(j))
                            pubdate{i}{I(k)} = pubDate{i}{j}; 
                            citations{i}(I(k)) = citationCount{i}(j); 
                            openaccess{i}(I(k)) = openAccess{i}(j); 
                        end
                    end
                end
            end
        end
    end
end
clearvars Year citationCount Title DOI sh_tList pubDate openAccess

%% Output refs

refs.allInfo = allInfo; 
refs.fullNames = fullNames; 
refs.affiliations = affiliations; 
refs.abstracts = abstracts; 
refs.titles = titles; 
refs.journals = journals; 
refs.keywords = keywords; 
refs.venues = venues; 
refs.doi = doi; 
refs.citations = citations; 
refs.pubdate = pubdate; 
refs.openaccess = openaccess; 

end
% saveAbstracts.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), May 2022
% This function creates and saves a data file entitled "abstracts.mat"

% DEPENDENCIES
% 1. Dependencies/Data/places.xlsx
% 2. Dependencies/Data/rmvRefs.mat
% 3. Dependencies/Data/nameBank.mat
% 4. Dependencies/Functions/simplifyName.m ** 
% 5. Dependencies/Functions/extractRefs.m ** 
%    ** see extractRefs.mat for that file's dependencies

clear
addpath(genpath('Dependencies/')); 

%% Extract abstracts:

years = [1969:2021]';
refs = extractRefs(years, 'rmvRefs.mat'); % this function pulls most of the heavy lifting

% reassign variables currently embedded as the structure created by extractRefs
allInfo = refs.allInfo; 
fullNames = refs.fullNames; 
affiliations = refs.affiliations; 
abstracts = refs.abstracts; 
titles = refs.titles; 
journals = refs.journals; 
keywords = refs.keywords; 
venues = refs.venues; 
citations = refs.citations; 
pubdate = refs.pubdate; 
openaccess = refs.openaccess; 

%% repeat authors if multiple affiliations / countries

% load in city list
[cityPop, cityList] = xlsread('places.xlsx','Cities'); % from here: https://simplemaps.com/data/world-cities
cityList = string(cityList); 

% load in country list
[~, countryList] = xlsread('places.xlsx','Countries');
countryList = string(countryList); countryList = countryList(:,1); 

% load in state list
[~, stateList] = xlsread('places.xlsx','States');
stateList = string(stateList); 

% load in university list
[~, univList] = xlsread('places.xlsx','Universities');
univList = string(univList); univList = univList(:,[1,3]); 
% remove (conventionally) ice core-irrelevant schools
wordage = ["Technical", "Community", "Medicine", "Theology", "Theological", "Business", "Design", "Bible", "Arts", ...
    "Seminary", "Christian", "Massage", "Acupuncture"," Law", "Health","Incarnate","Music"];
    for i = 1:length(wordage)
        Ic = false(size(univList,1),1); 
        for j = 1:length(univList)
            if contains(univList(j,1),wordage(i))
                Ic(j) = true; 
            end
        end
        univList(Ic,:) = []; 
    end

% load in manual input list
[~, manualList] = xlsread('places.xlsx','Manual');
manualList = string(manualList(:,1:3)); manualList(:,3) = []; 
    
% checking order: countries --> manual --> Universities 
country.country = cell(length(years),1);
country.author = cell(length(years),1);
country.affil = cell(length(years),1);
affiliationsFail = []; 
h = waitbar(0,'Estimating author home-countries ... please wait!'); 
for i = 1:length(affiliations)                            % i == year index
    m = 1; % reset the roaming index
    if ~isempty(affiliations{i})
    for j = 1:size(affiliations{i},1)                     % j == study index
        % some authors have more than one affiliation... we need to account for this!
        affiliationsCurr = []; 
        authorsCurr = []; 
        for k = 1:length(affiliations{i}{j})              % k == author index
            if ~ismissing(affiliations{i}{j}(k)) && length(affiliations{i}{j}) == length(fullNames{i}{j})
                newStr = split(affiliations{i}{j}(k),'; ');
                if length(newStr) > 1
                    affiliationsCurr = [affiliationsCurr; newStr]; 
                    authorsCurr = [authorsCurr; repmat(fullNames{i}{j}(k),[length(newStr),1])]; 
                else
                    affiliationsCurr = [affiliationsCurr; affiliations{i}{j}(k)]; 
                    authorsCurr = [authorsCurr; fullNames{i}{j}(k)]; 
                end
            end
        end
        for k = 1:length(affiliationsCurr)              % k == author index
            if ~isempty(affiliationsCurr)
                % COUNTRY
                Ic = false(size(countryList,1),1); 
                for g = 1:size(countryList,1)
                    if contains(affiliationsCurr(k),countryList(g,1))
                    Ic(g) = true; % state.State{i}{j}(k) = stateList(g,2); 
                    end
                end
                if sum(Ic) > 0 
                    Ic = find(Ic); Ic = Ic(randperm(length(Ic)));
                    country.country{i}{j,1}(k,1) = countryList(Ic(1));
                    try
                    country.author{i}{j,1}(k,1) = authorsCurr(k); 
                    country.affil{i}{j,1}(k,1) = affiliationsCurr(k); 
                    end
                elseif sum(Ic) == 0 % move onto the manual check 
                    % MANUAL 
                    Im = false(size(manualList,1),1); 
                    for g = 1:size(manualList,1)
                        if contains(affiliationsCurr(k),manualList(g,1))
                        Im(g) = true; 
                        end
                    end
                    if sum(Im) > 0 
                        Im = find(Im); Im = Im(randperm(length(Im)));
                        country.country{i}{j,1}(k,1) = manualList(Im(1),2); 
                        try
                        country.author{i}{j,1}(k,1) = authorsCurr(k); 
                        country.affil{i}{j,1}(k,1) = affiliationsCurr(k); 
                        end
                    elseif sum(Im) == 0 % move onto the US State check 
                        % STATE 
                        Is = false(size(stateList,1),1); 
                        for g = 1:size(stateList,1)
                            if contains(affiliationsCurr(k),stateList(g,1))
                            Is(g) = true; 
                            end
                        end
                        if sum(Is) > 0 
                            Is = find(Is); Is = Is(randperm(length(Is)));
                            country.country{i}{j,1}(k,1) = stateList(Is(1),3); 
                            try
                            country.author{i}{j,1}(k,1) = authorsCurr(k); 
                            country.affil{i}{j,1}(k,1) = affiliationsCurr(k); 
                            end
                        elseif sum(Is) == 0 % move onto the university-wide check 
                            % UNIVERSITY
                            Iu = false(size(univList,1),1); 
                            for g = 1:size(univList,1)
                                if contains(affiliationsCurr(k),univList(g,1))
                                Iu(g) = true; 
                                end
                            end
                            if sum(Iu) > 0 
                                Iu = find(Iu); Iu = Iu(randperm(length(Iu)));
                                country.country{i}{j,1}(k,1) = univList(Iu(1),2); 
                                try
                                country.author{i}{j,1}(k,1) = authorsCurr(k); 
                                country.affil{i}{j,1}(k,1) = affiliationsCurr(k); 
                                end
                            elseif sum(Iu) == 0 % move onto the City check 
                                % City
                                Ic = false(size(cityList,1),1); 
                                for g = 1:size(cityList,1)
                                    if contains(affiliationsCurr(k),cityList(g,1))
                                    Ic(g) = true; 
                                    end
                                end
                                if sum(Ic) > 0 
                                    Ic = find(Ic); 
                                    [~,Ii] = max(cityPop(Ic)); 
                                    Ic = Ic(Ii); % assumption: correct city as the one with the largest population
                                    country.country{i}{j,1}(k,1) = cityList(Ic,3); 
                                    try
                                    country.author{i}{j,1}(k,1) = authorsCurr(k); 
                                    country.affil{i}{j,1}(k,1) = affiliationsCurr(k); 
                                    end
                                else
                                    affiliationsFail = [affiliationsFail; affiliationsCurr(k)]; 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    end
    waitbar(i / length(years)); 
end
close(h);

% Change US/UK conventions, and remove missing values
for i = 1:length(affiliations)                            % i == year index
    if ~isempty(country.country{i})
    for j = 1:size(country.country{i},1)                  % j == study index
        for k = 1:size(country.country{i}{j},1)           % k == author index
            if strcmp(country.country{i}{j}(k),"USSR")
                country.country{i}{j}(k) = "Russia"; 
            end
            if strcmp(country.country{i}{j}(k),"U.S.S.R.")
                country.country{i}{j}(k) = "Russia"; 
            end
            if strcmp(country.country{i}{j}(k),"USA")
                country.country{i}{j}(k) = "United States"; 
            end
            if strcmp(country.country{i}{j}(k),"U.S.")
                country.country{i}{j}(k) = "United States"; 
            end
            if strcmp(country.country{i}{j}(k),"U.S.A.")
                country.country{i}{j}(k) = "United States"; 
            end
            if strcmp(country.country{i}{j}(k),"Jersey")
                country.country{i}{j}(k) = "United States"; 
            end
            if strcmp(country.country{i}{j}(k),"Georgia") % assumption! 
                country.country{i}{j}(k) = "United States"; 
            end
            if strcmp(country.country{i}{j}(k),"England")
                country.country{i}{j}(k) = "United Kingdom"; 
            end   
            if strcmp(country.country{i}{j}(k),"UK")
                country.country{i}{j}(k) = "United Kingdom"; 
            end   
            if strcmp(country.country{i}{j}(k),"U.K.")
                country.country{i}{j}(k) = "United Kingdom"; 
            end   
        end
        if ~isempty(country.country{i}{j})
        % remove missing values
        I = ismissing(country.country{i}{j});
        country.country{i}{j}(I) = []; 
        country.affil{i}{j}(I) = []; 
        country.author{i}{j}(I) = []; 
        end
    end
    end
end

% added by MBO for yeari = 29, Jan-10-2023 
for i = 1:length(affiliations)
    if ~isempty(affiliations{i})
        if isempty(affiliations{i}{end})
            country.country{i}{end+1,1} = {}; 
            country.author{i}{end+1,1} = {}; 
            country.affil{i}{end+1,1} = {}; 
        end
    end
end

%% Assign gender

[numbers,genderList] = xlsread('nameBank.xlsx','genderize'); 
genderList = string(genderList); 
probability = numbers(:,1); 

% preallocate "firstName" that is the same size as "fullNames" 
for i = 1:length(fullNames)
    firstNames{i,1} = cell(length(fullNames{i}),1); 
    gender{i,1} = cell(length(fullNames{i}),1); 
    prob{i,1} = cell(length(fullNames{i}),1); 
    for j = 1:length(fullNames{i})
        firstNames{i,1}{j,1} = strings(length(fullNames{i}{j}),1); 
        gender{i,1}{j,1} = strings(length(fullNames{i}{j}),1); 
        prob{i,1}{j,1} = nan(length(fullNames{i}{j}),1); 
    end
end

% assign first names
for i = 1:length(fullNames)
    if ~isempty(fullNames{i})
    for j = 1:length(fullNames{i})
        I = ~ismissing(fullNames{i}{j}); 
        firstNames{i}{j}(I) = extractBefore(fullNames{i}{j}(I),' ');
    end
    end
end

% now, go through and match the gender
for i = 1:length(firstNames)
    if ~isempty(firstNames{i})
    for j = 1:length(firstNames{i})
        if ~isempty(firstNames{i}{j})
        for k = 1:length(firstNames{i}{j})
            I = strcmp(firstNames{i}{j}(k),genderList(:,1)); 
            if sum(I) > 0
                gender{i}{j}(k,1) = genderList(I,2);
                prob{i}{j}(k,1) = probability(I);
            end
        end
        end
    end
    end
end

%% citation rate

currDate = datetime('16-Sep-2022'); % date citations were scraped via Semantic Scholar API

currDateDec = decyear(currDate); 
timeElapsed = cell(size(years)); 
for i = 1:length(citations)
    if ~isempty(citations{i})
        for j = 1:size(citations{i},1)
            if ~isempty(pubdate{i}{j})
                if ~isnat(pubdate{i}{j})
                    timeElapsed{i}(j,1) = currDateDec - decyear(pubdate{i}{j}); 
                else
                    timeElapsed{i}(j,1) = currDateDec - years(i)+0.5; 
                end
            else
            timeElapsed{i}(j,1) = currDateDec - years(i)+0.5; 
            end
        end
    end
end

% now, put another loop to correct for anything outside the expected range:
for i = 1:length(timeElapsed)
    if ~isempty(timeElapsed{i})
        for j = 1:size(timeElapsed{i},1)
            if abs([currDateDec - years(i)+0.5] - timeElapsed{i}(j)) > 2
                timeElapsed{i}(j) = currDateDec - years(i) + rand; 
            end
        end
    end
end

%% Convert full author names to "simple" name tokens
% name tokens contain the first character of each author's first same and 
% the full surname, all in lower case and separated by a single space

simpleNames = fullNames; 
h = waitbar(0, 'Simplifying names, please wait...');
for i = 1:size(fullNames,1)
    if ~isempty(fullNames{i})
    for j = 1:size(fullNames{i},1)
        if ~isempty(fullNames{i}{j})
        for k = 1:size(fullNames{i}{j},1)
            if ~ismissing(fullNames{i}{j}(k))
            simpleNames{i}{j}(k) = simplifyName(fullNames{i}{j}(k)); 
            end
        end
        end
    end
    end
    waitbar(i / size(fullNames,1))
end
close(h); 


%% Save or upload all the data I need to do things quickly

% save 
clearvars -except years allInfo fullNames affiliations abstracts titles ...
    journals keywords venues citations pubdate openaccess timeElapsed country gender prob simpleNames

explainer.abstracts = "study abstract text [53 x 1 cell]"; 
explainer.affiliations = "academic affiliations for each author [53 x 1 cell]"; 
explainer.allInfo = "the raw (unorganized) info for each study taken from the NASA ADS [53 x 1 cell]"; 
explainer.citations = "the number of citations for each study (as of Sept. 2022) [53 x 1 cell]"; 
explainer.country = "A struct containing three 53 x 1 cell variables: country -- the nationality represented by each author (estimated via academic affiliation); author -- the full author names, and affil -- the affiliation of each author"; 
explainer.fullNames = "the full name of each author on a per-study basis [53 x 1 cell]"; 
explainer.gender = "the genderize.io-estimated gender of each author [53 x 1 cell]"; 
explainer.journals = "the journal each study was published in [53 x 1 cell]"; 
explainer.keywords = "any study keywords submitted by the authors [53 x 1 cell]"; 
explainer.openacess = "a Boolean (true/false) designator; 1 == open access, 0 == not open access [53 x 1 cell]"; 
explainer.prob = "the probability of each genderize.io-derived gender assignment [53 x 1 cell]"; 
explainer.pubdate = "the date of publication for each study [53 x 1 cell]"; 
explainer.timeElapsed = "the amount of time elapsed since each study's publication (Sept. 2022 - pubdate) [53 x 1 cell]"; 
explainer.titles = "the title of each study [53 x 1 cell]"; 
explainer.venues = "the venue-type of each abstract; each in this database should be 'JOUR' (a journal publication) [53 x 1 cell]"; 
explainer.simpleNames = "name tokens containing the first character of each author's first same and the full surname, all in lower case and separated by a single space [53 x 1 cell]";
explainer.years = "a vector of years (CE) analyzed [53 x 1 double]"; 

save('abstracts.mat','-v7.3'); 

% savePersons.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), Jan 2023

% this script builds a unique "person" list designating all unique name tokens,
%   i.e., the first character of each author's first same and the full surname.

% DEPENDENCIES
% 1. abstracts.mat **
% ** generated by running saveAbstracts.m
clear

if ~isfile("abstracts.mat")
    error('Missing file "abstracts.mat". Please first run saveAbstracts.m')
else
    load('abstracts.mat')
end

%% Define each name as a "simple name", and then find all unique simple names

% find all unique persons
persons = []; 
for i = 1:size(simpleNames,1)
    if ~isempty(simpleNames{i})
        for j = 1:size(simpleNames{i},1)
        persons = [persons; simpleNames{i}{j}];
        end
    end
end
persons = unique(persons); 
persons(ismissing(persons)) = []; 

% Now, for each unique person, store attributes including study titles, ID'd gender,
%   + probability, ID'd countries, journal of each publication, number coauthors 
%   on each study, number of studies they were involved in, citation rate per study, 
%   year of each study, gender ratio of each study, gender of first author of each
%   study, whether that person was the first author, etc.
personStore.year = cell(size(persons,1),1);
personStore.cite = cell(size(persons,1),1);
personStore.coauthorC = cell(size(persons,1),1);
personStore.journals = cell(size(persons,1),1);
personStore.genderR = cell(size(persons,1),1);
personStore.genderF = cell(size(persons,1),1);
personStore.genderN = cell(size(persons,1),1);
personStore.genderP = cell(size(persons,1),1); % gender of the person of interest
personStore.nations = cell(size(persons,1),1);
personStore.firstAuthor = cell(size(persons,1),1);
personStore.prob = cell(size(persons,1),1);
personStore.coauthors = cell(size(persons,1),1);

h = waitbar(0, 'Assigning attributes to each unique author, please wait...');
for p = 1:length(persons)
    for i = 1:size(fullNames,1)
    if ~isempty(fullNames{i})
        for j = 1:size(fullNames{i},1)
        if ~isempty(fullNames{i}{j})
            I = strcmp(persons(p),simpleNames{i}{j});
            if sum(I) > 0
                personStore.year{p} = [personStore.year{p}; years(i)];
                personStore.cite{p} = [personStore.cite{p}; (citations{i}(j) ./ timeElapsed{i}(j))]; 
                personStore.coauthorC{p} = [personStore.coauthorC{p}; length(fullNames{i}{j}) - 1]; 
                personStore.journals{p} = [personStore.coauthors{p}; journals{i}{j}]; 
                personStore.genderR{p} = [personStore.genderR{p}; sum(strcmp(gender{i}{j},"female")) ./...
                        (sum(strcmp(gender{i}{j},"female")) + sum(strcmp(gender{i}{j},"male")))]; 
                personStore.genderF{p} = [personStore.genderF{p}; gender{i}{j}(1)];
                personStore.genderN{p} = [personStore.genderN{p}; sum(strcmp(gender{i}{j},"female"))]; 
                personStore.genderP{p} = [personStore.genderP{p}; gender{i}{j}(I)]; 
                personStore.nations{p} = [personStore.nations{p}; length(unique(country.country{i}{j}))];     
                personStore.firstAuthor{p} = [personStore.firstAuthor{p}; strcmp(simpleNames{i}{j}(1),persons(p))];     
                personStore.prob{p} = [personStore.prob{p}; prob{i}{j}(I)];   
                if length(simpleNames{i}{j}) > 1
                personStore.coauthors{p} = [personStore.coauthors{p}; simpleNames{i}{j}(2:end)];  
                end
            end
        end
        end
    end
    end
    waitbar(p / size(persons,1))
end
close(h); 

% calculate number of studies attributed to each author
personStore.nStudies = nan(size(persons,1),1);
personStore.persons = persons;
for p = 1:length(persons)
    personStore.nStudies(p) = length(personStore.year{p});
end

% calculate prbability associated with each person
%   note that each probability listed will be the same! taking mean just a 
%   simple way to condense down to a single number if and where it exists 
personStore.probability = nan(size(persons,1),1);
for p = 1:length(persons)
    I = ~isnan(personStore.prob{p});
    personStore.probability(p) = nanmean(personStore.prob{p}(I)); end % condense repeating values
personStore = rmfield(personStore,'prob'); 

% now, define a gender for each person
for  p = 1:length(persons)
    isMale = sum(strcmp(personStore.genderP{p},"male")); 
    isFemale = sum(strcmp(personStore.genderP{p},"female")); 
    if isMale > isFemale
    personStore.gender(p,1) = "male"; 
    elseif isFemale > isMale
    personStore.gender(p,1) = "female"; 
    elseif isFemale == isMale
    personStore.gender(p,1) = ""; 
    end
end

% clean up and rename output
personStore = rmfield(personStore,"genderP"); 
persons = personStore; clearvars personStore

% Finally, run a bunch of Monte Carlo gender guess simulations
mc = 1e4;
maleFemale = ["male","female"]; % preallocate MC matrix
persons.genderMC = strings(size(persons.persons,1),mc); 
h = waitbar(0, 'Running Monte Carlo gender guess simulations... please wait...');
for i = 1:size(persons.persons,1)
    for m = 1:mc
        if ~strcmp(persons.gender(i),"")
            r = rand;
            if r > persons.probability(i)
                persons.genderMC(i,m) = maleFemale(~strcmp(persons.gender(i),maleFemale)); 
            else
                persons.genderMC(i,m) = persons.gender(i); 
            end
        end
    end
    waitbar(i / size(persons.persons,1))
end
close(h);

% finally define the variables and save output, yo! 
varGuide = ["nations = number of nations represented besides first-author's nation";
            "coauthorC = number of coauthors per study";
            "journals = corresponding journal of each published";
            "cite = citation rate per year";
            "genderR = inferred gender ratio";
            "genderN = number of female coauthors"
            "genderF = first author gender"
            "nStudies = number of studies published by each person"
            "gender = inferred person gender"
            "genderMC = monte carlo resampled person gender"
            "persons = unique name identifier"
            "firstAuthor = a T/F identifier; ; if 1, then person is first author"
            "prop = genderize-derived gender assignment probability"
            "coauthors = unique name ID's of all coauthors"];
persons.varGuide = varGuide; 
save('persons.mat','persons','-v7.3'); 
clearvars -except persons
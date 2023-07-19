% EDF2.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), Apr 2023

% this file generates Extended Data Figure 2 of Koffman et al. (in revision at Nature Geoscience)
% showing proportion of studies led by women across leading ice core journals

clear
homeFold = cd; 
cd ../; workFold = cd;
% check on dependency availability
if ~isfile("Dependencies/Data/Fig1_mc.mat")
    error('Missing file "Fig1_mc.mat". Please run Fig1.m first')
end

%% User-defined input parameters:

mc = 1000; % number of monte carlo gender simulations (suggested mc >= 1000)
timeBins = [1969, 2011; 
            2012, 2021]; % input the time bins over which to compare -- make sure the size(timeBins,1) = 2 

%% Load the abstract data

if isfile("abstracts.mat")
    load('abstracts.mat')
else
    error('Missing file "abstracts.mat". Please first run saveAbstracts.m')
end

%% Estimate initialed author genders

if isfile("persons.mat")
    load('persons.mat')
else
    error('Missing file "persons.mat". Please first run savePersons.m')
end
genderP = cell(size(simpleNames)); 
probability = cell(size(simpleNames)); 
h = waitbar(0, 'Matching initialed authors ... please wait...');
for i = 1:size(simpleNames,1)
    if ~isempty(simpleNames{i})
    genderP{i} = cell(size(simpleNames{i})); 
    probability{i} = cell(size(simpleNames{i})); 
    for j = 1:size(simpleNames{i},1)
        if ~isempty(simpleNames{i}{j})
        genderP{i}{j} = strings(size(simpleNames{i}{j})); 
        probability{i}{j} = nan(size(simpleNames{i}{j})); 
        for k = 1:size(simpleNames{i}{j},1)
            if ~strcmp(gender{i}{j}(k),"") % if gender information already exists in 'gender', then use that info!
                genderP{i}{j}(k) = gender{i}{j}(k); 
                probability{i}{j}(k) = prob{i}{j}(k); 
            else % otherwise, try to use prior coauthorship overlap to see if we can infer it
                I = strcmp(simpleNames{i}{j}(k),persons.persons);
                In = 1:size(simpleNames{i}{j},1); In = find(In~=k);
                if ~isempty(In)
                    if sum(I) > 0 % if there's a name match...
                        % find if there exists at least 1 coauthor overlap
                        if ~isempty(intersect(persons.coauthors{I},simpleNames{i}{j}(In)))
                            genderP{i}{j}(k) = persons.gender(I);
                            probability{i}{j}(k) = persons.probability(I);
                        end
                    end
                end
            end
        end
        end
    end
    end
    waitbar(i / size(simpleNames,1))
end
close(h)
gender = genderP; clearvars genderP
prob = probability; clearvars probability

%% Determine metrics for each journal

% place impact factor

keyJour = ["Annals of Glaciology"
        "Atmospheric Chemistry & Physics"
        "Atmospheric Environment"
        "Climate Dynamics"
        "Climate of the Past"
        "Earth and Planetary Science Letters"
        "Environmental Science and Technology"
        "Geochimica et Cosmochimica Acta"
        "Geophysical Research Letters"
        "Global and Planetary Change"
        "Global Biogeochemical Cycles"
        "Journal of Climate"
        "Journal of Geophysical Research"
        "Journal of Glaciology"
        "Nature"
        "Nature Geoscience"
        "Paleoceanography"
        "Proceedings of the National Academy of Science"
        "Quaternary Science Reviews"
        "Science"
        "Tellus Series B Chemical and Physical Meteorology B"];
IF = [2.863
      7.197
      4.798
      4.901
      4.498
      5.785
      11.357
      5.921
      5.58
      5.114
      5.96
      5.148
      4.261
      4.278
      69.5
      21.53
      3.99
      12.78
      4.571
      63.71
      2.279];

refJour = erase(string(keyJour)," "); 
refJour = erase(refJour,"&"); 

for r = 1:length(refJour)
    jour.(refJour(r)).type = [];
    jour.(refJour(r)).year = [];
    jour.(refJour(r)).nations = [];
    jour.(refJour(r)).coauthors = [];
    jour.(refJour(r)).cite = [];
    jour.(refJour(r)).oa = [];
    jour.(refJour(r)).genderR = []; % gender ratio
    jour.(refJour(r)).genderN = []; % number of females 
    jour.(refJour(r)).genderF = []; % gender of first author
    jour.(refJour(r)).prob = []; % probability
end

h = waitbar(0,'Assigning journal metrics ... please wait');
for r = 1:length(refJour)
    jour.(refJour(r)).IF = IF(r);
    for i = 1:length(years)
        for j = 1:length(journals{i})
            if strcmp(journals{i}{j},keyJour(r))
                ca = unique(country.country{i}{j}); 
                if length(ca) == 1 % if the coauthor countries are the same as the first author' country.. 
                    jour.(refJour(r)).type = [jour.(refJour(r)).type; "domestic"]; 
                elseif length(ca) > 1
                    jour.(refJour(r)).type = [jour.(refJour(r)).type; "international"]; 
                else
                    jour.(refJour(r)).type = [jour.(refJour(r)).type; ""]; 
                end
                % year
                jour.(refJour(r)).year = [jour.(refJour(r)).year; years(i)];
                % open access
                jour.(refJour(r)).oa = [jour.(refJour(r)).oa; openaccess{i}(j)];
                % cite rate 
                jour.(refJour(r)).cite = [jour.(refJour(r)).cite; (citations{i}(j) ./ timeElapsed{i}(j))];
                % number of nationalies
                if isempty(ca) % if the coauthor countries are the same as the first author' country.. 
                    jour.(refJour(r)).nations = [jour.(refJour(r)).nations; 0]; 
                else
                    jour.(refJour(r)).nations = [jour.(refJour(r)).nations; length(ca)]; 
                end            
                % number of coauthors
                jour.(refJour(r)).coauthors = [jour.(refJour(r)).coauthors; length(country.country{i}{j})]; 
                % gender
                if ~isempty(gender{i}{j})
                    % gender ratio 
                    jour.(refJour(r)).genderR = [jour.(refJour(r)).genderR; sum(strcmp(gender{i}{j},"female")) ./...
                        (sum(strcmp(gender{i}{j},"female")) + sum(strcmp(gender{i}{j},"male")))]; 
                    % number of females
                    jour.(refJour(r)).genderN = [jour.(refJour(r)).genderN; sum(strcmp(gender{i}{j},"female"))]; 
                    % gender of first author
                    jour.(refJour(r)).genderF = [jour.(refJour(r)).genderF; gender{i}{j}(1)]; 
                    % probability
                    jour.(refJour(r)).prob = [jour.(refJour(r)).prob; prob{i}{j}(1)]; 
                else
                    jour.(refJour(r)).genderR = [jour.(refJour(r)).genderR; NaN]; 
                    jour.(refJour(r)).genderN = [jour.(refJour(r)).genderN; NaN]; 
                    jour.(refJour(r)).genderF = [jour.(refJour(r)).genderF; ""];                 
                    jour.(refJour(r)).prob = [jour.(refJour(r)).prob; NaN];                 
                end
            end
        end
    end
    waitbar(r / length(refJour))
end
close(h);

%% Calculate uncertainty

maleFemale = ["male","female"];
h = waitbar(0,'conducting MC experiments ... please wait');
for g = 1:length(refJour)
    jour.(refJour(g)).genderMC = strings(size(jour.(refJour(g)).genderF,1),mc); % gender first author
    for k = 1:mc
        for i = 1:size(jour.(refJour(g)).genderMC,1)
                % first author
                genderFirst = jour.(refJour(g)).genderF(i);
                if ~strcmp(genderFirst,"")
                    r = rand;
                if r > jour.(refJour(g)).prob(i); genderFirst = maleFemale(~strcmp(genderFirst,maleFemale)); end
            end
            jour.(refJour(g)).genderMC(i,k) = genderFirst; 
        end
    end
    waitbar(g / length(refJour))
end
close(h)

%% Figure -- compare two periods without career stage

cd(workFold)

div = timeBins(1,2); 
barWidth = 0.30; 
K = length(refJour);
[IF,Isort] = sort(IF,'descend');
refJour = refJour(Isort); 

% define figure params
h = figure; hold on;
set(h,'units','centimeters','position',[1.8,0.90,15.5,18.5]);
    ax = gca; ax.Visible = 'off';
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');

% assign colors
cd Figures/cbrewer/
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,5); 
    CT2 = cbrewer('seq','Blues' ,5);
    CT3 = cbrewer('seq','Purples' ,5);    
    CT4 = cbrewer('seq','Greys' ,5);    
    CT5 = cbrewer('qual','Set1' ,length(refJour)); CT5(CT5<0) = 0; CT5(CT5>1) = 1; CT5 = flipud(CT5); 
warning('on','all')
cd ../../

% calculate the proportion of women authors for each journal 

    propFemale = nan(length(refJour),2);
    for r = 1:length(refJour)
        % proportion 
        I = jour.(refJour(r)).year <= div;
        propFemale(r,1) = nansum(strcmp(jour.(refJour(r)).genderF(I),'female')) ./ ...
            (nansum(strcmp(jour.(refJour(r)).genderF(I),'female')) + nansum(strcmp(jour.(refJour(r)).genderF(I),'male')));
        propFemale(r,2) = nansum(strcmp(jour.(refJour(r)).genderF(~I),'female')) ./ ...
            (nansum(strcmp(jour.(refJour(r)).genderF(~I),'female')) + nansum(strcmp(jour.(refJour(r)).genderF(~I),'male')));
        % error range
        currError = nansum(strcmp(jour.(refJour(r)).genderMC(I,:),'female')) ./ ...
            (nansum(strcmp(jour.(refJour(r)).genderMC(I,:),'female'),1) + nansum(strcmp(jour.(refJour(r)).genderMC(I,:),'male'),1));
        jour.(refJour(r)).error(1,:) = [prctile(currError,2.5),prctile(currError,25),prctile(currError,50),prctile(currError,75),prctile(currError,97.5)]; 
        currError = nansum(strcmp(jour.(refJour(r)).genderMC(~I,:),'female')) ./ ...
            (nansum(strcmp(jour.(refJour(r)).genderMC(~I,:),'female'),1) + nansum(strcmp(jour.(refJour(r)).genderMC(~I,:),'male'),1));
        jour.(refJour(r)).error(2,:) = [prctile(currError,2.5),prctile(currError,25),prctile(currError,50),prctile(currError,75),prctile(currError,97.5)]; 
    end

    xLabels = refJour;
    % rename journals manually
    xLabels(3) = "Nature Geoscience"; 
    xLabels(4) = "Proc. Nat. Acad. Sci."; 
    xLabels(5) = "Env. Sci. & Technology"; 
    xLabels(6) = "Atmos. Chem. & Phys."; 
    xLabels(7) = "Glob. & Planetary Change"; 
    xLabels(8) = "Geochem. Cosm. Acta"; 
    xLabels(9) = "Earth & Plan. Sci. Letters"; 
    xLabels(10) = "Geophys. Res. Letters"; 
    xLabels(11) = "J. Climate"; 
    xLabels(12) = "Glob. Biogeochem. Cycles"; 
    xLabels(13) = "Climate Dynamics"; 
    xLabels(14) = "Atmos. Env."; 
    xLabels(15) = "Quat. Sci. Reviews"; 
    xLabels(16) = "Climate of the Past"; 
    xLabels(17) = "J. Glaciology"; 
    xLabels(18) = "J. Geophysical Research"; 
    xLabels(20) = "Annals of Glaciology"; 
    xLabels(21) = "Tellus Series B";     

    ax1 = axes('position',[0.40 0.15 0.20 0.70]); hold on; grid off; 
        set(ax1,'Color','none','linewidth',1.0,'fontsize',10,'Box','off','xaxislocation','bottom','xdir','reverse'); hold on;
        set(ax1,'Xtick',[1:K],'Xticklabel',xLabels,'XTickLabelRotation',0,'Ytick',0:.20:1.00);
        ylabel([{'Proportion studies'},{'woman-led'}]); 
        xlim([0 K+1]); 

    % gender parity period one
    load Fig1_mc.mat
    genAll = nan(size(timeBins,1),1); 
    for i = 1:size(timeBins,1)
        I = years >= min(timeBins(i,:)) & years <= max(timeBins(i,:)) & ~isnan(nanmean(grA,2));
        genAll(i,1) = nansum(numStudies(I) .* nanmean(grA(I,:),2)) ./ sum(numStudies(I));
    end
    p0 = plot([0 K+1],genAll(1).*[1 1],':','linewidth',1.5,'Color',CT1(2,:));
    % gender parity period two
    p0 = plot([0 K+1],genAll(2).*[1 1],':','linewidth',1.5,'Color',CT1(4,:));
    % gender parity
    p0 = plot([0 K+1],[0.5 0.5],':','linewidth',1.5,'Color',CT4(end-2,:));

    for i = 1:length(refJour)
        b1 = bar(ax1,[i-barWidth/2],jour.(refJour(i)).error(1,3),'BarWidth',barWidth,'LineWidth',1); 
            b1(1).FaceColor = CT1(2,:); b1(1).EdgeColor = CT4(end-1,:); 
        p1 = plot(ax1,[i-barWidth/2, i-barWidth/2],[jour.(refJour(i)).error(1,1), jour.(refJour(i)).error(1,end)],'LineWidth',1,'Color','k');  
        b2 = bar(ax1,[i+barWidth/2],jour.(refJour(i)).error(2,3),'stacked','BarWidth',barWidth,'LineWidth',1); 
            b2(1).FaceColor = CT1(4,:); b2(1).EdgeColor = CT4(end-1,:); 
        p2 = plot(ax1,[i+barWidth/2, i+barWidth/2],[jour.(refJour(i)).error(2,1), jour.(refJour(i)).error(2,end)],'LineWidth',1,'Color','k');  

    end

set(gca,'view',[90 -90])

% total number studies per journal

    nStudies = nan(length(refJour),2);
    for r = 1:length(refJour)
        I = jour.(refJour(r)).year <= div;
        nStudies(r,1) = nansum(strcmp(jour.(refJour(r)).genderF(I),'female')) + nansum(strcmp(jour.(refJour(r)).genderF(I),'male')); 
        nStudies(r,2) = nansum(strcmp(jour.(refJour(r)).genderF(~I),'female')) + nansum(strcmp(jour.(refJour(r)).genderF(~I),'male'));
    end

    ax1 = axes('position',[0.65 0.15 0.20 0.70]); hold on; grid off; 
        set(ax1,'Color','none','linewidth',1.0,'fontsize',10,'Box','off','xaxislocation','bottom','xdir','reverse'); hold on;
        set(ax1,'Xtick',[],'Xticklabel',[],'YScale','linear');
        ylabel([{'Number studies'}]); 
        ylim([1 350]); 	
        xlim([0 K+1]); 

    for i = 1:length(refJour)
        b1 = bar(ax1,[i-barWidth/2],[nStudies(i,1)],'BarWidth',barWidth,'LineWidth',1); 
            b1(1).FaceColor = CT1(2,:); b1(1).EdgeColor = CT4(end-1,:); 
        b2 = bar(ax1,[i+barWidth/2],[nStudies(i,2)],'stacked','BarWidth',barWidth,'LineWidth',1); 
            b2(1).FaceColor = CT1(4,:); b2(1).EdgeColor = CT4(end-1,:); 
    end

    set(gca,'view',[90 -90])

% add legend
l1 = legend(ax1,[b1 b2],[join(string(timeBins(1,:)),'-')],[join(string(timeBins(2,:)),'-')],'Box','off','Orientation','Vertical'); 
    l1.Position = [ 0.0775    0.8775    0.2299    0.0521]; % [0.2963    0.2366    0.2314    0.0462]; 
    l1.ItemTokenSize = [12,8]; % default is [30,18] 
    
% add some descriptors
ax = axes('position', [0.46 0.87 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[{'Community'},{'representation'}],'Fontsize',9,'Fontweight','normal','color',CT1(3,:),'HorizontalAlignment','left','Rotation',30);
ax = axes('position', [0.55 0.86 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,['Gender parity'],'Fontsize',9,'Fontweight','normal','color',CT4(end-2,:),'HorizontalAlignment','left','Rotation',30); 
ax = axes('position', [0.23 0.865 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,'Higher Impact Factor','Fontsize',10,'Fontweight','normal','color',CT4(end-2,:),'HorizontalAlignment','center','FontAngle','italic');
ax = axes('position', [0.23 0.135 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,'Lower Impact Factor','Fontsize',10,'Fontweight','normal','color',CT4(end-2,:),'HorizontalAlignment','center','FontAngle','italic');

cd(homeFold); 

% Fig3.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), Apr 2023

% this file generates Figure 3 of Koffman et al. (in revision at Nature Geoscience)
% exploring the shifting nature of women coauthor networks in ice core science across gender and career stage

clear

homeFold = cd; 
cd ../; workFold = cd;
addpath(genpath('Dependencies/')); 

%% User-defined input parameters

mc = 1000; % number of monte carlo gender simulations (suggested mc >= 1000)
timeBins = [1979 1999; % choose two time periods only, min/max must be >= 1979 and <= 2011, respectively
            2000 2011]; % input the time bins over which to compare -- make sure the size(timeBins,1) = 2,

%% Load the abstract data

if isfile("abstracts.mat")
    load('abstracts.mat')
else
    error('Missing file "abstracts.mat". Please first run saveAbstracts.m')
end

%% Load persons data

if isfile("persons.mat")
    load('persons.mat')
else
    error('Missing file "persons.mat". Please first run savePersons.m')
end

%% Refine dataset

% define each author's first year of occurrence and then remove
% all authors with < minNumPapers, and without at least one occurrence in
% at least one junior and one senior scientist designation 

% define each persons career start
persons.yearStart = nan(size(persons.year)); 
for i = 1:size(persons.year)
    persons.yearStart(i) = nanmin(persons.year{i});
end

% define the minimum number of papers per author
minNumPapers = 2; 

% remove all non-identified genders and all persons with <minNumPapers studies
I = strcmp(persons.gender,"") | persons.nStudies < minNumPapers;
persons.year(I) = [];
persons.cite(I) = [];
persons.coauthors(I) = [];
persons.coauthorC(I) = [];
persons.journals(I) = [];
persons.genderR(I) = [];
persons.genderF(I) = [];
persons.genderN(I) = [];
persons.nations(I) = [];
persons.firstAuthor(I) = [];
persons.nStudies(I) = [];
persons.probability(I) = [];
persons.persons(I) = [];
persons.gender(I) = [];
persons.genderMC(I,:) = [];
persons.yearStart(I,:) = [];

%% determine where in the career stage each author is! 

% go through all persons and, for each study, define that person as one of the following:
% 1 = junior man
% 2 = senior man
% 3 = junior woman
% 4 = senior woman

if mc > size(persons.genderMC,2); mc = size(persons.genderMC,2); end % check that mc isn't too large!
nYears = 10; % number of years assumed between junior and senior author
cStage = cell(size(simpleNames)); 
h = waitbar(0, 'Assigning career stage ... please wait...');
for i = 1:size(simpleNames,1)
    if ~isempty(simpleNames{i})
    cStage{i} = cell(size(simpleNames{i})); 
    for j = 1:size(simpleNames{i},1)
        if ~isempty(simpleNames{i}{j})
        cStage{i}{j} = nan(size(simpleNames{i}{j},1),mc); 
        for k = 1:size(simpleNames{i}{j},1)
            I = strcmp(simpleNames{i}{j}(k),persons.persons);
            if sum(I) > 0 % if there's a name match...
            for m = 1:mc
                if ((years(i) - persons.yearStart(I)) < nYears) && strcmp(persons.genderMC(I,m),"male")
                    cStage{i}{j}(k,m) = 1;
                elseif ((years(i) - persons.yearStart(I)) >= nYears) && strcmp(persons.genderMC(I,m),"male")
                    cStage{i}{j}(k,m) = 2;
                elseif ((years(i) - persons.yearStart(I)) < nYears) && strcmp(persons.genderMC(I,m),"female")
                    cStage{i}{j}(k,m) = 3;
                elseif ((years(i) - persons.yearStart(I)) >= nYears) && strcmp(persons.genderMC(I,m),"female")
                    cStage{i}{j}(k,m) = 4;
                end
            end
            end
        end
        end
    end
    end
    waitbar(i / size(cStage,1))
end
close(h)

% sort by first author -- one structure
jm = cell(size(timeBins,1),m); % junior man
sm = cell(size(timeBins,1),m); % senior man
jw = cell(size(timeBins,1),m); % junior woman
sw = cell(size(timeBins,1),m); % senior woman
for g = 1:size(timeBins,1)
    yearsC = find(years >= nanmin(timeBins(g,:)) & years <= nanmax(timeBins(g,:)));
    h = waitbar(0, 'Assigning coauthor networks ... please wait...');
    for i = min(yearsC):max(yearsC)
        if ~isempty(cStage{i})
        for j = 1:size(cStage{i},1)
            if length(cStage{i}{j}) > 1 % make sure more than 1 author
            for m = 1:mc
                if cStage{i}{j}(1,m) == 1 
                    jm{g,m} = [jm{g,m}; cStage{i}{j}(2:end,m)];
                elseif cStage{i}{j}(1,m) == 2 
                    sm{g,m} = [sm{g,m}; cStage{i}{j}(2:end,m)];
                elseif cStage{i}{j}(1,m) == 3 
                    jw{g,m} = [jw{g,m}; cStage{i}{j}(2:end,m)];
                elseif cStage{i}{j}(1,m) == 4 
                    sw{g,m} = [sw{g,m}; cStage{i}{j}(2:end,m)];
                end
            end
            end
        end
        end
        waitbar((i - min(yearsC) + 1) / (max(yearsC) - min(yearsC) + 1));
    end
    close(h)
end

%% calculate each of the author proportions and coauthor proportions

% calculate the proportion of each coauthor type for each period
allAuth = cell(size(timeBins,1),1); % senior woman
for g = 1:size(timeBins,1)
    m1 = 1; 
    yearsC = find(years >= nanmin(timeBins(g,:)) & years <= nanmax(timeBins(g,:)));
    for i = min(yearsC):max(yearsC)
        if ~isempty(cStage{i})
        for j = 1:size(cStage{i},1)
            m2 = size(cStage{i}{j},1); 
            if m2 > 1 % make sure more than 1 author
                allAuth{g}(m1:(m2+m1-1),:) = cStage{i}{j};
                % allAuth{g} = [allAuth{g}; cStage{i}{j}];
                m1 = size(allAuth{g},1) + 1;
            end
        end
        end
    end
    I = isnan(allAuth{g}(:,1)); 
    allAuth{g}(I,:) = []; 
end

authProp = nan(4,size(allAuth,1),mc);
h = waitbar(0, 'Calculating Monte Carlo proportions ... please wait...');
for m = 1:mc
    for i = 1:4
        for j = 1:size(allAuth,1)
            authProp(i,j,m) = nansum(allAuth{j}(:,m) == i) ./ size(allAuth{j},1);
        end
    end
    waitbar(m / mc)
end
close(h);
authProp = cat(2,authProp,authProp(:,2,:) - authProp(:,1,:));

%% Define coauthor proportion by gender and career stage

% reminder:
% 1 = junior man (jm)
% 2 = senior man (sm)
% 3 = junior woman (jw)
% 4 = senior woman (sw)

labels = ["jm","sm","jw","sw"]; 
% create a 3D matrix with the proportions 
coauthProp = cell(size(timeBins,1),1); % senior woman
for i = 1:length(coauthProp)
    coauthProp{i} = nan(length(labels),length(labels),mc);
end
for k = 1:size(timeBins,1)
    for m = 1:mc
        for j = 1:length(labels)
            coauthProp{k}(1,j,m) = sum(jm{k,m} == j) ./ sum(~isnan(jm{k,m})); % row1
            coauthProp{k}(2,j,m) = sum(sm{k,m} == j) ./ sum(~isnan(sm{k,m})); % row2
            coauthProp{k}(3,j,m) = sum(jw{k,m} == j) ./ sum(~isnan(jw{k,m})); % row3
            coauthProp{k}(4,j,m) = sum(sw{k,m} == j) ./ sum(~isnan(sw{k,m})); % row4
        end
    end
end
coauthPropDiff = coauthProp{end} - coauthProp{1};

%% Plot 

cd(workFold)

fig1 = figure; hold on;
    set(fig1,'units','centimeters','position',[2.5,0.5,11,20]);
    ax = gca; ax.Visible = 'off';
    set(fig1,'PaperPositionMode','auto');         
    set(fig1,'PaperOrientation','landscape');

% assign colors
cd Figures/cbrewer/
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,5); 
    CT2 = cbrewer('seq','YlOrBr' ,8); CT2(end,:) = []; CT2(end,:) = []; CT2(end,:) = []; CT2(end,:) = []; 
    CT3 = cbrewer('seq','Blues' ,5); 
    CT4 = cbrewer('seq','Greys' ,5);  
warning('on','all')
cd ../../

axisWidth = 0.60; 
K = size(authProp,2); 
barWidth = 1/6; 
fs = 10; % fontsize
upp = 75; % upper uncertainty range
low = 25; % lower uncertainty range

% Proportion of all authors:

ax0 = axes('Position',[0.15 0.80 axisWidth 0.125]); hold on; 
    set(ax0,'Color','none','linewidth',1.0,'fontsize',10,'Box','on','xaxislocation','bottom','yaxislocation','left',...
        'xtick',[1:K]','Xticklabel',[join(string(timeBins),'-'); "Difference"]); hold on; % ,'xcolor','none'
    ylabel([{'Proportion'},{'all authors'}],'fontsize',10); 
    ylim([-0.2 0.60]); 	
    xlim([0+3*barWidth K+1-3*barWidth]);
    % dividers
    plotVals = nanmean([ [1:(K-1)]',[2:(K)]'],2);
    for i = 1:K-1
        plot([plotVals(i),plotVals(i)],[-0.2 0.60],'--','linewidth',0.5,'color',CT4(end-1,:)); 
    end

    for i = 1:K
        b1 = bar(ax0,[i-3*barWidth/2],[median(authProp(1,i,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b1.FaceColor = CT2(2,:); b1.EdgeColor = CT4(end,:);
        p1 = plot(ax0,[i-3*barWidth/2, i-3*barWidth/2],[prctile(authProp(1,i,:),low), prctile(authProp(1,i,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b2 = bar(ax0,[i-1*barWidth/2],[median(authProp(2,i,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b2.FaceColor = CT2(end-1,:); b2.EdgeColor = CT4(end,:);
        p2 = plot(ax0,[i-1*barWidth/2, i-1*barWidth/2],[prctile(authProp(2,i,:),low), prctile(authProp(2,i,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b3 = bar(ax0,[i+1*barWidth/2],[median(authProp(3,i,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b3.FaceColor = CT1(2,:); b3.EdgeColor = CT4(end,:);
        p3 = plot(ax0,[i+1*barWidth/2, i+1*barWidth/2],[prctile(authProp(3,i,:),low), prctile(authProp(3,i,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b4 = bar(ax0,[i+3*barWidth/2],[median(authProp(4,i,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b4.FaceColor = CT1(end-1,:); b4.EdgeColor = CT4(end,:);
        p4 = plot(ax0,[i+3*barWidth/2, i+3*barWidth/2],[prctile(authProp(4,i,:),low), prctile(authProp(4,i,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
    end

    l1 = legend([b1 b2 b3 b4],'Junior man',...
        'Senior man','Junior woman','Senior woman','Box','off','Orientation','Vertical','Fontsize',8); 
        l1.Position = [0.7517    0.8421    0.2265    0.0873]; 
        l1.ItemTokenSize = [8,8]; % default is [30,18]

% redefine K
K = size(coauthPropDiff,2); 
xticklabels = ["Junior man-led"; "Senior man-led"; "Junior woman-led"; "Senior woman-led"];

% Proportion coauthor pre-2000

ax1 = axes('Position',[0.15 0.60 axisWidth*4/3 0.125]); hold on; 
    set(ax1,'Color','none','linewidth',1.0,'fontsize',10,'Box','on','xaxislocation','bottom','yaxislocation','left',...
        'xtick',[1:K],'xticklabel',[]); hold on; % ,'xcolor','none'
    ylabel([{'Proportion'},{'coauthors'}],'fontsize',10); 
    ylim([0 0.60]); 	
    xlim([0+3*barWidth K+1-3*barWidth]);
    % dividers
    plotVals = nanmean([ [1:(K-1)]',[2:(K)]'],2);
    for i = 1:K-1
        plot([plotVals(i),plotVals(i)],[-0.2 0.60],'--','linewidth',0.5,'color',CT4(end-1,:)); 
    end
    title(join(string(timeBins(1,:)),'-'),'Fontsize',10,'FontWeight','normal')

    for i = 1:K
        b1 = bar(ax1,[i-3*barWidth/2],[median(coauthProp{1}(i,1,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b1.FaceColor = CT2(2,:); b1.EdgeColor = CT4(end,:);
        p1 = plot(ax1,[i-3*barWidth/2, i-3*barWidth/2],[prctile(coauthProp{1}(i,1,:),low), prctile(coauthProp{1}(i,1,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b2 = bar(ax1,[i-1*barWidth/2],[median(coauthProp{1}(i,2,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b2.FaceColor = CT2(end-1,:); b2.EdgeColor = CT4(end,:);
        p2 = plot(ax1,[i-1*barWidth/2, i-1*barWidth/2],[prctile(coauthProp{1}(i,2,:),low), prctile(coauthProp{1}(i,2,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b3 = bar(ax1,[i+1*barWidth/2],[median(coauthProp{1}(i,3,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b3.FaceColor = CT1(2,:); b3.EdgeColor = CT4(end,:);
        p3 = plot(ax1,[i+1*barWidth/2, i+1*barWidth/2],[prctile(coauthProp{1}(i,3,:),low), prctile(coauthProp{1}(i,3,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b4 = bar(ax1,[i+3*barWidth/2],[median(coauthProp{1}(i,4,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b4.FaceColor = CT1(end-1,:); b4.EdgeColor = CT4(end,:);
        p4 = plot(ax1,[i+3*barWidth/2, i+3*barWidth/2],[prctile(coauthProp{1}(i,4,:),low), prctile(coauthProp{1}(i,4,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
    end

% Proportion coauthor pre-2000

ax2 = axes('Position',[0.15 0.44 axisWidth*4/3 0.125]); hold on; 
    set(ax2,'Color','none','linewidth',1.0,'fontsize',10,'Box','on','xaxislocation','bottom','yaxislocation','left',...
        'xtick',[1:K],'xticklabel',[]); hold on; % ,'xcolor','none'
    ylabel([{'Proportion'},{'coauthors'}],'fontsize',10); 
    ylim([0 0.60]); 	
    xlim([0+3*barWidth K+1-3*barWidth]);
    % dividers
    plotVals = nanmean([ [1:(K-1)]',[2:(K)]'],2);
    for i = 1:K-1
        plot([plotVals(i),plotVals(i)],[-0.2 0.60],'--','linewidth',0.5,'color',CT4(end-1,:)); 
    end
    title(join(string(timeBins(2,:)),'-'),'Fontsize',10,'FontWeight','normal')

    for i = 1:K
        b1 = bar(ax2,[i-3*barWidth/2],[median(coauthProp{2}(i,1,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b1.FaceColor = CT2(2,:); b1.EdgeColor = CT4(end,:);
        p1 = plot(ax2,[i-3*barWidth/2, i-3*barWidth/2],[prctile(coauthProp{2}(i,1,:),low), prctile(coauthProp{2}(i,1,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b2 = bar(ax2,[i-1*barWidth/2],[median(coauthProp{2}(i,2,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b2.FaceColor = CT2(end-1,:); b2.EdgeColor = CT4(end,:);
        p2 = plot(ax2,[i-1*barWidth/2, i-1*barWidth/2],[prctile(coauthProp{2}(i,2,:),low), prctile(coauthProp{2}(i,2,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b3 = bar(ax2,[i+1*barWidth/2],[median(coauthProp{2}(i,3,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b3.FaceColor = CT1(2,:); b3.EdgeColor = CT4(end,:);
        p3 = plot(ax2,[i+1*barWidth/2, i+1*barWidth/2],[prctile(coauthProp{2}(i,3,:),low), prctile(coauthProp{2}(i,3,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b4 = bar(ax2,[i+3*barWidth/2],[median(coauthProp{2}(i,4,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b4.FaceColor = CT1(end-1,:); b4.EdgeColor = CT4(end,:);
        p4 = plot(ax2,[i+3*barWidth/2, i+3*barWidth/2],[prctile(coauthProp{2}(i,4,:),low), prctile(coauthProp{2}(i,4,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
    end

% Proportion coauthor difference

ax3 = axes('Position',[0.15 0.28 axisWidth*4/3 0.125]); hold on; 
    set(ax3,'Color','none','linewidth',1.0,'fontsize',10,'Box','on','xaxislocation','bottom','yaxislocation','left',...
        'xtick',[1:K],'xticklabel',[]); hold on; % ,'xcolor','none'
    ylabel([{'\DeltaProportion'},{'coauthors'}],'fontsize',10); 
    ylim([-0.15 0.15]); set(ax3,'ytick',-0.2:0.1:0.2)
    xlim([0+3*barWidth K+1-3*barWidth]);
    % labels
    adj = -0.245; 
    text(1,adj,[{"Junior"},{"man-led"},{"studies"}],'color',CT4(end,:),'HorizontalAlignment','Center','Fontsize',10); 
    text(2,adj,[{"Senior"},{"man-led"},{"studies"}],'color',CT4(end,:),'HorizontalAlignment','Center','Fontsize',10); 
    text(3,adj,[{"Junior"},{"woman-led"},{"studies"}],'color',CT4(end,:),'HorizontalAlignment','Center','Fontsize',10); 
    text(4,adj,[{"Senior"},{"woman-led"},{"studies"}],'color',CT4(end,:),'HorizontalAlignment','Center','Fontsize',10); 
    % dividers
    plotVals = nanmean([ [1:(K-1)]',[2:(K)]'],2);
    for i = 1:K-1
        plot([plotVals(i),plotVals(i)],[-0.2 0.60],'--','linewidth',0.5,'color',CT4(end-1,:)); 
    end
    title([strcat('Difference (post minus pre-',num2str(timeBins(end,1)),')')],'Fontsize',10,'FontWeight','normal')

    for i = 1:K
        b1 = bar(ax3,[i-3*barWidth/2],[median(coauthPropDiff(i,1,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b1.FaceColor = CT2(2,:); b1.EdgeColor = CT4(end,:);
        p1 = plot(ax3,[i-3*barWidth/2, i-3*barWidth/2],[prctile(coauthPropDiff(i,1,:),low), prctile(coauthPropDiff(i,1,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b2 = bar(ax3,[i-1*barWidth/2],[median(coauthPropDiff(i,2,:))],'BarWidth',barWidth,'LineWidth',0.5); 
            b2.FaceColor = CT2(end-1,:); b2.EdgeColor = CT4(end,:);
        p2 = plot(ax3,[i-1*barWidth/2, i-1*barWidth/2],[prctile(coauthPropDiff(i,2,:),low), prctile(coauthPropDiff(i,2,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b3 = bar(ax3,[i+1*barWidth/2],[median(coauthPropDiff(i,3,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b3.FaceColor = CT1(2,:); b3.EdgeColor = CT4(end,:);
        p3 = plot(ax3,[i+1*barWidth/2, i+1*barWidth/2],[prctile(coauthPropDiff(i,3,:),low), prctile(coauthPropDiff(i,3,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
        b4 = bar(ax3,[i+3*barWidth/2],[median(coauthPropDiff(i,4,:))],'stacked','BarWidth',barWidth,'LineWidth',0.5); 
            b4.FaceColor = CT1(end-1,:); b4.EdgeColor = CT4(end,:);
        p4 = plot(ax3,[i+3*barWidth/2, i+3*barWidth/2],[prctile(coauthPropDiff(i,4,:),low), prctile(coauthPropDiff(i,4,:),upp)],'LineWidth',1,'Color',CT4(end,:));  
    end

ax = axes('position', [0.20 0.8175 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[strcat("{\itn} = ",num2str(size(allAuth{1},1)))],'Fontsize',8,'Fontweight','normal','color',CT4(end-1,:),'HorizontalAlignment','left','FontAngle','normal');
ax = axes('position', [0.40 0.8175 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[strcat("{\itn} = ",num2str(size(allAuth{2},1)))],'Fontsize',8,'Fontweight','normal','color',CT4(end-1,:),'HorizontalAlignment','left','FontAngle','normal');

ax = axes('position', [0.89 0.745 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,['(b)'],'Fontsize',14,'Fontweight','bold','color',CT4(end,:),'HorizontalAlignment','left','FontAngle','normal');
ax = axes('position', [0.70 0.945 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[{'(a)'}],'Fontsize',14,'Fontweight','bold','color',CT4(end,:),'HorizontalAlignment','left','FontAngle','normal');

cd(homeFold)

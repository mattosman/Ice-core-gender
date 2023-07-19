% EDF1.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), Apr 2023

% this file generates Extended Data Figure 1 of Koffman et al. (in revision at Nature Geoscience)
% showing patterns in citation rate, coauthorship, and internationality by gender

clear
homeFold = cd; 
cd ../; workFold = cd;
addpath(genpath('Dependencies/')); 

%% User-defined input parameters

timeBins = [1969, 2011; 
            2012, 2021]; % input the time bins over which to compare -- make sure the size(timeBins,1) = 2 or 3

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

%% Determine metrics by each study where first author gender is denoted

study.type = [];
study.year = [];
study.nations = [];
study.coauthors = [];
study.cite = [];
study.genderR = []; % gender ratio
study.genderN = []; % number of females 
study.genderF = []; % gender of first author

for i = 1:length(years)
    for j = 1:length(gender{i})
        if ~isempty(gender{i}{j})
            if (gender{i}{j}(1) == "") == 0
                % country
                ca = unique(country.country{i}{j}); 
                if length(ca) == 1 % if the coauthor countries are the same as the first author' country.. 
                    study.type = [study.type; "domestic"]; 
                elseif length(ca) > 1
                    study.type = [study.type; "international"]; 
                else
                    study.type = [study.type; ""]; 
                end
                % year
                study.year = [study.year; years(i)];
                % cite rate 
                study.cite = [study.cite; log(citations{i}(j) ./ timeElapsed{i}(j))];
                % number of nationalities
                if isempty(ca) % if the coauthor countries are the same as the first author country.. 
                    study.nations = [study.nations; 0]; 
                else
                    study.nations = [study.nations; length(ca)]; 
                end            
                % number of coauthors
                study.coauthors = [study.coauthors; length(country.country{i}{j}(2:end))]; 
                % gender ratio
                if sum(~strcmp(gender{i}{j},"")) >= 2 && sum(~strcmp(gender{i}{j},""))/length(gender{i}{j}) >= 1.0 % at least more than one gender listed and at least half genders reported
                    study.genderR = [study.genderR; sum(strcmp(gender{i}{j}(2:end),"female")) ./...
                        (sum(strcmp(gender{i}{j}(2:end),"female")) + sum(strcmp(gender{i}{j}(2:end),"male")))]; 
                else
                    study.genderR = [study.genderR; NaN]; 
                end
                % number of females
                study.genderN = [study.genderN; sum(strcmp(gender{i}{j}(2:end),"female"))]; 
                % gender of first author
                study.genderF = [study.genderF; gender{i}{j}(1)]; 
            end
        end
    end
end

% preallocate
male.type = cell(size(timeBins,1),1); 
male.nations = cell(size(timeBins,1),1); 
male.coauthors = cell(size(timeBins,1),1); 
male.cite = cell(size(timeBins,1),1); 
male.genderR = cell(size(timeBins,1),1); 
male.genderN = cell(size(timeBins,1),1); 
female = male; % copy preallocated cells
% partition by man vs. woman led
for i = 1:size(timeBins,1)
    Im = study.year >= timeBins(i,1) & study.year <= timeBins(i,2) & strcmp(study.genderF,"male");
    If = study.year >= timeBins(i,1) & study.year <= timeBins(i,2) & strcmp(study.genderF,"female");
    male.type{i} = study.type(Im); female.type{i} = study.type(If); 
    male.nations{i} = study.nations(Im); female.nations{i} = study.nations(If); 
    male.cite{i} = study.cite(Im); female.cite{i} = study.cite(If); 
        I = isinf(male.cite{i}) | isnan(male.cite{i}); male.cite{i}(I) = []; % remove nan or inf val's
        I = isinf(female.cite{i}) | isnan(female.cite{i}); female.cite{i}(I) = []; % remove nan or inf val's
    male.coauthors{i} = study.coauthors(Im); female.coauthors{i} = study.coauthors(If); 
    male.genderR{i} = study.genderR(Im); female.genderR{i} = study.genderR(If); 
    male.genderN{i} = study.genderN(Im); female.genderN{i} = study.genderN(If); 
end

%% Plot

cd(workFold) 

fig1 = figure; hold on;
    if size(timeBins,1) == 3
    set(fig1,'units','centimeters','position',[0.5,0.5,13,15]);
    elseif size(timeBins,1) == 2
    set(fig1,'units','centimeters','position',[0.5,0.5,10.5,15]);
    end
    ax = gca; ax.Visible = 'off';
    set(fig1,'PaperPositionMode','auto');         
    set(fig1,'PaperOrientation','landscape');

% assign colors
cd Figures/cbrewer/
warning('off','all')  
    CT1 = cbrewer('seq','YlOrBr' ,9); CT1(end,:) = []; CT1(end,:) = []; CT1(end,:) = []; CT1(1,:) = []; 
    CT2 = cbrewer('seq','Reds' ,5);
    CT3 = cbrewer('seq','Purples' ,5);    
    CT4 = cbrewer('seq','Greys' ,5);  
    CT5 = cbrewer('seq','Blues' ,5);  
    gold = [0.975, 0.65, 0.125]; 
warning('on','all')
cd ../../

K = size(timeBins,1); 
adj = 0.025; 
barWidth = 0.40; 
widthVal = 0.20; 
box1Alpha = 0.50; 

ax0 = axes('position',[0.15 0.10 0.675 0.75]); hold on; grid off; 
    set(ax0,'Color','none','linewidth',1.0,'fontsize',10,'Box','on','xtick',[],'ytick',[]); hold on;
    plotVals = nanmean([ [1:(K-1)]',[2:(K)]'],2); 
    ylim([0 1]); 	
    xlim([0+0.5 K+1-0.5]);
    for i = 1:size(timeBins,1)-1
        plot([plotVals(i),plotVals(i)],[0 1],'--','linewidth',0.5,'color',CT4(end-1,:)); 
    end

% Proportion studies international 
ax1 = axes('position',[0.15 0.10 0.675 0.10]); hold on; grid off; 
    set(ax1,'Color','none','linewidth',1.0,'fontsize',10,'Box','off','xaxislocation','bottom','yaxislocation','left'); hold on;
    set(ax1,'Xtick',[1:K],'Xticklabel',join(string(timeBins),'-'),'XTickLabelRotation',0);
    ylabel([{'Number studies'}]); 
    xlabel(ax1,[{'Year range'}]); 
    if size(timeBins,1) == 3
    ylim([0 900]); set(ax1,'Ytick',[0:300:900]); 
    elseif size(timeBins,1) == 2
    ylim([0 1500]); set(ax1,'Ytick',[0:500:1500]); 
    end
    xlim([0+0.5 K+1-0.5]);

    for i = 1:size(timeBins,1)
        % male
        b1 = bar(ax1,[i-barWidth/2],[length(male.type{i})],'stacked','BarWidth',barWidth/4,'LineWidth',1); 
            b1(1).FaceColor = CT4(end-3,:); b1(1).EdgeColor = CT4(end-1,:); 
        % female
        b2 = bar(ax1,[i+barWidth/2],[length(female.type{i})],'stacked','BarWidth',barWidth/4,'LineWidth',1); 
            b2(1).FaceColor = CT4(end-3,:); b2(1).EdgeColor = CT4(end-1,:); 
    end

% number nations per study
ax2 = axes('position',[0.15 0.22 0.675 0.20]); hold on; grid off; 
    set(ax2,'Color','none','linewidth',1.0,'fontsize',10,'Box','off','xaxislocation','bottom','yaxislocation','right','xcolor','none'); hold on;
    set(ax2,'Xtick',[1:K],'Xticklabel',join(string(timeBins),'-'),'XTickLabelRotation',0);
    ylabel([{'Number nations'},{'per study'}]); 
    ylim([0 6]); set(ax2,'Ytick',[0:2:6]); 
    xlim([0+0.5 K+1-0.5]);

    for i = 1:size(timeBins,1)
        % men
        xLoc = i-barWidth/2; cDat = male.nations{i}; cDat(isnan(cDat)) = [];
            % violin plot
            h = figure(101); cd Figures/Violinplot/; % procedure: run violin plot to get the plotting data:
            violins = violinplot(cDat,[0],'Width',widthVal,'Bandwidth',0.35); % ,'Bandwidth',10);
            vertices = violins.ViolinPlot.Vertices; 
                upper = [vertices(1:(length(vertices)/2),1)-1, vertices(1:(length(vertices)/2),2)]; 
                lower = [vertices((length(vertices)/2+1):end,1)-1, vertices((length(vertices)/2+1):end,2)]; 
            close(h); cd ../../
            box1a = fill(ax2,[lower(:,1) + xLoc; upper(:,1) + xLoc] + adj, [lower(:,2); upper(:,2)],CT1(2,:)); 
                box1a.EdgeColor = 'k'; box1a.LineWidth = 0.5; box1a.FaceAlpha = box1Alpha; box1a.EdgeAlpha = box1Alpha + 0.2; 
            % 95 percentile range
            l = plot(ax2,[xLoc xLoc] + adj,[prctile(cDat,2.5), prctile(cDat,97.5)],'Linewidth',1,'Color',0.3.*[1 1 1]);
            % Interquartile range
            box2a = fill(ax2,[xLoc-widthVal/4, xLoc-widthVal/4, xLoc+widthVal/4, xLoc+widthVal/4] + adj,...
                [prctile(cDat,25), prctile(cDat,75), prctile(cDat,75), prctile(cDat,25)],CT1(3,:)); 
                box2a.EdgeColor = 0.3.*[1 1 1]; box2a.FaceAlpha = 1.0; 
            % mean
            l = plot(ax2,[xLoc-widthVal/4, xLoc+widthVal/4] + adj,[nanmean(cDat), nanmean(cDat)],'Linewidth',1.5,'Color',CT1(end,:));
        % women
        xLoc = i+barWidth/2; cDat = female.nations{i}; cDat(isnan(cDat)) = [];
            % violin plot
            h = figure(101); cd Figures/Violinplot/; % procedure: run violin plot to get the plotting data:
            violins = violinplot(cDat,[0],'Width',widthVal,'Bandwidth',0.35); % ,'Bandwidth',10);
            vertices = violins.ViolinPlot.Vertices; 
                upper = [vertices(1:(length(vertices)/2),1)-1, vertices(1:(length(vertices)/2),2)]; 
                lower = [vertices((length(vertices)/2+1):end,1)-1, vertices((length(vertices)/2+1):end,2)]; 
            close(h); cd ../../
            box1b = fill(ax2,[lower(:,1) + xLoc; upper(:,1) + xLoc] - adj, [lower(:,2); upper(:,2)],CT2(2,:)); 
                box1b.EdgeColor = 'k'; box1b.LineWidth = 0.5; box1b.FaceAlpha = box1Alpha; box1b.EdgeAlpha = box1Alpha + 0.2; 
            % 95 percentile range
            l = plot(ax2,[xLoc xLoc] - adj,[prctile(cDat,2.5), prctile(cDat,97.5)],'Linewidth',1,'Color',0.3.*[1 1 1]);
            % Interquartile range
            box2b = fill(ax2,[xLoc-widthVal/4, xLoc-widthVal/4, xLoc+widthVal/4, xLoc+widthVal/4] - adj,...
                [prctile(cDat,25), prctile(cDat,75), prctile(cDat,75), prctile(cDat,25)],CT2(3,:)); 
                box2b.EdgeColor = 0.3.*[1 1 1]; box2b.FaceAlpha = 1.0; 
            % mean
            l = plot(ax2,[xLoc-widthVal/4, xLoc+widthVal/4] - adj,[nanmean(cDat), nanmean(cDat)],'Linewidth',1.5,'Color',CT2(end,:));
    end

% citation rate

% number nations per study
ax2 = axes('position',[0.15 0.44 0.675 0.20]); hold on; grid off; 
    set(ax2,'Color','none','linewidth',1.0,'fontsize',10,'Box','off','xaxislocation','bottom','yaxislocation','left','xcolor','none'); hold on;
    set(ax2,'Xtick',[1:K],'Xticklabel',join(string(timeBins),'-'),'XTickLabelRotation',0);
    ylabel([{'Annual citation'},{'rate per study'},{'(log-transformed)'}]); 
    ylim([-3 5]); set(ax2,'Ytick',[-2:2:4]); 
    xlim([0+0.5 K+1-0.5]);

    for i = 1:size(timeBins,1)
        % men
        xLoc = i-barWidth/2; cDat = male.cite{i}; cDat(isnan(cDat)) = [];
            % violin plot
            h = figure(101); cd Figures/Violinplot/; % procedure: run violin plot to get the plotting data:
            violins = violinplot(cDat,[0],'Width',widthVal); % ,'Bandwidth',10);
            vertices = violins.ViolinPlot.Vertices; 
                upper = [vertices(1:(length(vertices)/2),1)-1, vertices(1:(length(vertices)/2),2)]; 
                lower = [vertices((length(vertices)/2+1):end,1)-1, vertices((length(vertices)/2+1):end,2)]; 
            close(h); cd ../../
            box1a = fill(ax2,[lower(:,1) + xLoc; upper(:,1) + xLoc] + adj, [lower(:,2); upper(:,2)],CT1(2,:)); 
                box1a.EdgeColor = 'k'; box1a.LineWidth = 0.5; box1a.FaceAlpha = box1Alpha; box1a.EdgeAlpha = box1Alpha + 0.2; 
            % 95 percentile range
            l = plot(ax2,[xLoc xLoc] + adj,[prctile(cDat,2.5), prctile(cDat,97.5)],'Linewidth',1,'Color',0.3.*[1 1 1]);
            % Interquartile range
            box2a = fill(ax2,[xLoc-widthVal/4, xLoc-widthVal/4, xLoc+widthVal/4, xLoc+widthVal/4] + adj,...
                [prctile(cDat,25), prctile(cDat,75), prctile(cDat,75), prctile(cDat,25)],CT1(3,:)); 
                box2a.EdgeColor = 0.3.*[1 1 1]; box2a.FaceAlpha = 1.0; 
            % mean
            l = plot(ax2,[xLoc-widthVal/4, xLoc+widthVal/4] + adj,[nanmean(cDat), nanmean(cDat)],'Linewidth',1.5,'Color',CT1(end,:));
        % women
        xLoc = i+barWidth/2; cDat = female.cite{i}; cDat(isnan(cDat)) = [];
            % violin plot
            h = figure(101); cd Figures/Violinplot/; % procedure: run violin plot to get the plotting data:
            violins = violinplot(cDat,[0],'Width',widthVal); % ,'Bandwidth',10);
            vertices = violins.ViolinPlot.Vertices; 
                upper = [vertices(1:(length(vertices)/2),1)-1, vertices(1:(length(vertices)/2),2)]; 
                lower = [vertices((length(vertices)/2+1):end,1)-1, vertices((length(vertices)/2+1):end,2)]; 
            close(h); cd ../../
            box1b = fill(ax2,[lower(:,1) + xLoc; upper(:,1) + xLoc] - adj, [lower(:,2); upper(:,2)],CT2(2,:)); 
                box1b.EdgeColor = 'k'; box1b.LineWidth = 0.5; box1b.FaceAlpha = box1Alpha; box1b.EdgeAlpha = box1Alpha + 0.2; 
            % 95 percentile range
            l = plot(ax2,[xLoc xLoc] - adj,[prctile(cDat,2.5), prctile(cDat,97.5)],'Linewidth',1,'Color',0.3.*[1 1 1]);
            % Interquartile range
            box2b = fill(ax2,[xLoc-widthVal/4, xLoc-widthVal/4, xLoc+widthVal/4, xLoc+widthVal/4] - adj,...
                [prctile(cDat,25), prctile(cDat,75), prctile(cDat,75), prctile(cDat,25)],CT2(3,:)); 
                box2b.EdgeColor = 0.3.*[1 1 1]; box2b.FaceAlpha = 1.0; 
            % mean
            l = plot(ax2,[xLoc-widthVal/4, xLoc+widthVal/4] - adj,[nanmean(cDat), nanmean(cDat)],'Linewidth',1.5,'Color',CT2(end,:));
    end

% number of coauthors by gender

ax3 = axes('position',[0.15 0.65 0.675 0.20]); hold on; grid off; 
    set(ax3,'Color','none','linewidth',1.0,'fontsize',10,'Box','off','xaxislocation','bottom','yaxislocation','right','xcolor','none'); hold on;
    set(ax3,'Xtick',[1:K],'Xticklabel',join(string(timeBins),'-'),'XTickLabelRotation',0);
    ylabel([{'Number woman'},{'coauthors per study'}]); 
    % ylim([0 20]); set(ax2,'Ytick',[0:4:20]); 
    ylim([0 6.5]); set(ax3,'Ytick',[0:2:6]); 
    xlim([0+0.5 K+1-0.5]);

    for i = 1:size(timeBins,1)
        % men
        xLoc = i-barWidth/2; cDat = male.genderN{i}; cDat(isnan(cDat)) = [];
            % violin plot
            h = figure(101); cd Figures/Violinplot/; % procedure: run violin plot to get the plotting data:
            violins = violinplot(cDat,[0],'Width',widthVal,'Bandwidth',0.325); % ,'Bandwidth',10);
            vertices = violins.ViolinPlot.Vertices; 
                upper = [vertices(1:(length(vertices)/2),1)-1, vertices(1:(length(vertices)/2),2)]; 
                lower = [vertices((length(vertices)/2+1):end,1)-1, vertices((length(vertices)/2+1):end,2)]; 
            close(h); cd ../../
            box1a = fill(ax3,[lower(:,1) + xLoc; upper(:,1) + xLoc] + adj, [lower(:,2); upper(:,2)],CT1(2,:)); 
                box1a.EdgeColor = 'k'; box1a.LineWidth = 0.5; box1a.FaceAlpha = box1Alpha; box1a.EdgeAlpha = box1Alpha + 0.2; 
            % 95 percentile range
            l = plot(ax3,[xLoc xLoc] + adj,[prctile(cDat,2.5), prctile(cDat,97.5)],'Linewidth',1,'Color',0.3.*[1 1 1]);
            % Interquartile range
            box2a = fill(ax3,[xLoc-widthVal/4, xLoc-widthVal/4, xLoc+widthVal/4, xLoc+widthVal/4] + adj,...
                [prctile(cDat,25), prctile(cDat,75), prctile(cDat,75), prctile(cDat,25)],CT1(3,:)); 
                box2a.EdgeColor = 0.3.*[1 1 1]; box2a.FaceAlpha = 1.0; 
            % mean
            l = plot(ax3,[xLoc-widthVal/4, xLoc+widthVal/4] + adj,[nanmean(cDat), nanmean(cDat)],'Linewidth',1.5,'Color',CT1(end,:));
        % women
        xLoc = i+barWidth/2; cDat = female.genderN{i}; cDat(isnan(cDat)) = [];
            % violin plot
            h = figure(101); cd Figures/Violinplot/; % procedure: run violin plot to get the plotting data:
            violins = violinplot(cDat,[0],'Width',widthVal,'Bandwidth',0.325); % ,'Bandwidth',10);
            vertices = violins.ViolinPlot.Vertices; 
                upper = [vertices(1:(length(vertices)/2),1)-1, vertices(1:(length(vertices)/2),2)]; 
                lower = [vertices((length(vertices)/2+1):end,1)-1, vertices((length(vertices)/2+1):end,2)]; 
            close(h); cd ../../
            box1b = fill(ax3,[lower(:,1) + xLoc; upper(:,1) + xLoc] - adj, [lower(:,2); upper(:,2)],CT2(2,:)); 
                box1b.EdgeColor = 'k'; box1b.LineWidth = 0.5; box1b.FaceAlpha = box1Alpha; box1b.EdgeAlpha = box1Alpha + 0.2; 
            % 95 percentile range
            l = plot(ax3,[xLoc xLoc] - adj,[prctile(cDat,2.5), prctile(cDat,97.5)],'Linewidth',1,'Color',0.3.*[1 1 1]);
            % Interquartile range
            box2b = fill(ax3,[xLoc-widthVal/4, xLoc-widthVal/4, xLoc+widthVal/4, xLoc+widthVal/4] - adj,...
                [prctile(cDat,25), prctile(cDat,75), prctile(cDat,75), prctile(cDat,25)],CT2(3,:)); 
                box2b.EdgeColor = 0.3.*[1 1 1]; box2b.FaceAlpha = 1.0; 
            % mean
            l = plot(ax3,[xLoc-widthVal/4, xLoc+widthVal/4] - adj,[nanmean(cDat), nanmean(cDat)],'Linewidth',1.5,'Color',CT2(end,:));
    end

ax = axes('position', [0.16 0.895 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,['Man-led'],'Fontsize',10,'Fontweight','normal','color',CT1(end-1,:),'HorizontalAlignment','left','FontAngle','normal');
ax = axes('position', [0.16 0.87 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[{'Woman-led'}],'Fontsize',10,'Fontweight','normal','color',CT2(end-1,:),'HorizontalAlignment','left','FontAngle','normal');

cd(homeFold) 

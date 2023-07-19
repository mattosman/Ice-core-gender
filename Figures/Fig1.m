% Fig1.m
% written by Matt Osman (mattosman@arizona.edu / mo549@cam.ac.uk), Apr 2023

% this file generates Figure 1 of Koffman et al. (in revision at Nature Geoscience)
% highlighting trends in women authorship in ice core science over the past half-century

clear
homeFold = cd; 
cd ../; workFold = cd;
addpath(genpath('Dependencies/')); 

%% User-defined input parameters:

mc = 1000; % number of monte carlo gender simulations (suggested mc >= 1000)

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

%% calculate gender ratio

if isfile("Fig1_mc.mat")
    load Fig1_mc.mat
    fileExists = true;
else
    fileExists = false;
    grA = []; 
end

% compute + save output only if mc doesn't exist big or mc is specified as larger than what already exists
if ~fileExists || size(grA,2) < mc

numStudies = zeros(size(years)); 
grF = nan(size(years,1),mc); 
grC = nan(size(years,1),mc); 
grA = nan(size(years,1),mc); 
maleFemale = ["male","female"];
h = waitbar(0, 'Creating random gender permutations...');
for k = 1:mc
    genderF = cell(size(gender)); % gender first author
    genderC = cell(size(gender)); % gender of all coauthors
    genderA = cell(size(gender)); % gender of all authors
    for i = 1:length(gender)
        m = 1; % roaming index for number of studies
        if ~isempty(abstracts{i})
        for j = 1:size(gender{i},1)
            if ~isempty(gender{i}{j})
                % first author
                genderFirst = gender{i}{j}(1);
                if ~strcmp(genderFirst,"")
                    r = rand;
                    if r > prob{i}{j}(1); genderFirst = maleFemale(~strcmp(genderFirst,maleFemale)); end
                end
                genderF{i} = [genderF{i}; genderFirst]; 
                % coauthor
                if size(gender{i}{j},1) > 1
                genderCoauthor = gender{i}{j}(2:end);
                for c = 1:numel(genderCoauthor)
                    if ~strcmp(genderCoauthor(c),"")
                        r = rand;
                        % if random value is greater than prob value, then assign opposite gender as genderize
                        if r > prob{i}{j}(c+1); genderCoauthor(c) = maleFemale(~strcmp(genderCoauthor(c),maleFemale)); end 
                    end
                end
                genderC{i} = [genderC{i}; genderCoauthor]; 
                end
                % all author
                genderA{i} = [genderA{i}; [genderFirst;genderCoauthor]]; 
                if sum(~strcmpi(gender{i}{j}(1),"")) > 0
                numStudies(i) = m; 
                m = m+1;
                end
            end
        end
        I = ~strcmp(genderF{i},""); 
        grF(i,k) = sum(strcmp(genderF{i}(I),"female")) / ...
                (sum(strcmp(genderF{i}(I),"female")) + sum(strcmp(genderF{i}(I),"male"))); % .*100;
        I = ~strcmp(genderC{i},""); 
        grC(i,k) = sum(strcmp(genderC{i}(I),"female")) / ...
                (sum(strcmp(genderC{i}(I),"female")) + sum(strcmp(genderC{i}(I),"male"))); % .*100;
        I = ~strcmp(genderA{i},""); 
        grA(i,k) = sum(strcmp(genderA{i}(I),"female")) / ...
                (sum(strcmp(genderA{i}(I),"female")) + sum(strcmp(genderA{i}(I),"male"))); % .*100);
        end
    end
    waitbar(k / mc)
end
close(h)
    
% save output only if adequately big
save('Dependencies/Data/Fig1_mc.mat','grA','grC','grF','mc','numStudies','maleFemale','-v7.3'); 

end

%% %%%%%%%%%%%%% Plot it up %%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(workFold); 

% set figure parameters
fig1 = figure; hold on;
    set(fig1,'units','centimeters','position',[2.5,0.5,20,13.5]);
    ax = gca; ax.Visible = 'off';
    set(fig1,'PaperPositionMode','auto');         
    set(fig1,'PaperOrientation','landscape');

    topWidth = 0.60;
    btmHt = 0.275;
    totHt = 0.5; 
    spaceBetween = 0.125; % space between bottom plots
    bottomWidth = (topWidth - spaceBetween) / 2; % 0.025 == space between plots
    startWidth = 0.20; 
    startWidth(2) = startWidth(1) + bottomWidth + spaceBetween; 
    panel1Lim = [0 0.55];
    panel2Lim = [-0.35 0.30];
   
% assign colors
cd Figures/cbrewer/
warning('off','all')
    CT1 = cbrewer('seq','Reds' ,6); 
    CT2 = cbrewer('seq','Blues' ,5);
    CT3 = cbrewer('seq','Purples' ,5);    
    CT4 = cbrewer('seq','Greys' ,5);  
    gold = [0.975, 0.65, 0.125]; 
warning('on','all')
cd ../../
markerSize = 1.2; 

% background axis
ax = axes('position',[startWidth(1) 0.25 topWidth totHt]); hold on; grid off; 
	set(ax,'Color','none','linewidth',1.0,'fontsize',9,'Box','on','xtick',[],'ytick',[]); hold on;

% gender proportion axis
ax = axes('position',[startWidth(1) 0.25 topWidth btmHt]); hold on; grid off; 
	set(ax,'Color','none','linewidth',1.0,'fontsize',9,'Box','off','xaxislocation','bottom'); hold on;
    xlim([min(years) max(years)+1]); 
    ylim([panel1Lim]); 
    set(gca,'ytick',[0:0.1:0.4])
    ylabel('Proportion women','fontsize',10); % xlabel('Year','fontsize',11); 
    % set(gca,'YTick',log([1,10,100]),'YTickLabel',[{'1'},{'10'},{'100'}],'yminortick','on'); 

% genderize uncertainties
for i = 1:length(years)
    if ~isnan(grF(i,1)) && numStudies(i) > 0
    l1 = plot([years(i) years(i)],[prctile(grF(i,:),2.5), prctile(grF(i,:),97.5)], 'Linewidth', 0.5, 'Color', 0.3.*[1 1 1] ); %  CT1(end-1,:)
    end
    if ~isnan(grA(i,1)) && numStudies(i) > 0
    l2 = plot([years(i) years(i)],[prctile(grA(i,:),2.5), prctile(grA(i,:),97.5)], 'Linewidth', 0.5, 'Color', 0.3.*[1 1 1] ); % CT2(end-1,:)
    end
end
 
% scatter
I = ~isnan(grF(:,1)) & numStudies > 0; 
s1 = scatter(years(I),nanmean(grF(I,:),2), markerSize .* numStudies(I),'filled','Linewidth',0.5); 
    s1.MarkerFaceColor = CT1(end-2,:); s1.MarkerEdgeColor = CT4(end,:); 
I = ~isnan(grA(:,1)) & numStudies > 0;
s2 = scatter(years(I),nanmean(grA(I,:),2), markerSize .* numStudies(I),'filled','Linewidth',0.5); 
    s2.MarkerFaceColor = CT1(end-0,:); s2.MarkerEdgeColor = CT4(end,:); 

% legend
offAxis = 1000; 
p1 = scatter(ax, offAxis, offAxis, 30, CT1(end-0,:), 'filled');          
	p1.Marker = 'o'; p1.MarkerEdgeColor = CT4(end,:); p1.MarkerFaceAlpha = 1.0; p1.LineWidth = 0.5; 
p2 = scatter(ax, offAxis, offAxis, 30, CT1(end-2,:), 'filled');          
	p2.Marker = 'o'; p2.MarkerEdgeColor = CT4(end,:); p2.MarkerFaceAlpha = 1.0; p2.LineWidth = 0.5; 
p3 = scatter(ax, offAxis, offAxis, 30, CT1(end-4,:), 'filled'); p2.Marker = 'o'; % 'square';         
	p3.Marker = 'o'; p3.MarkerEdgeColor = CT4(end,:); p3.MarkerFaceAlpha = 1.0; p3.LineWidth = 0.5; 
[leg,icons] = legend([p1 p2 p3], 'All author', 'First author', 'First author minus coauthor',...
    'Orientation','vertical','Fontsize',9,'Box','off','NumColumns',2);
    leg.Position = [0.39    0.765    0.1896    0.0789]; 

% gender excess proportion axis
ax = axes('position',[startWidth(1) (0.25 + (totHt - diff(panel2Lim)*btmHt/diff(panel1Lim))) topWidth (diff(panel2Lim)*btmHt/diff(panel1Lim))]); hold on; grid off; 
	set(ax,'Color','none','linewidth',1.0,'fontsize',9,'Box','off','xaxislocation','bottom','yaxislocation','right','xcolor','none'); hold on;
    xlim([min(years) max(years)+1]); 
    ylim(panel2Lim); 
    set(gca,'ytick',[-0.2:0.1:0.3])
    ylabel('\DeltaProportion women','fontsize',10); % xlabel('Year','fontsize',11); 

% plot the genderize uncertainties
for i = 1:length(years)
    if ~isnan(grF(i,1)) && ~isnan(grC(i,1)) && numStudies(i) > 0
    l3 = plot([years(i) years(i)],[prctile(sort(grF(i,:),'ascend') - sort(grC(i,:),'ascend'),2.5), ...
        prctile(sort(grF(i,:),'ascend') - sort(grC(i,:),'ascend'),97.5)], 'Linewidth', 0.5, 'Color', 0.3.*[1 1 1] ); % gold
    end
end
% line across zero
l4 = plot([years(1) years(end)],[0 0], 'Linewidth', 1, 'Color', 0.3.*[1 1 1],'Linestyle','-.'); 
 
% scatter
I = ~isnan(grF(:,1)) & numStudies > 0; 
s3 = scatter(years(I), nanmean(grF(I,:)-grC(I,:),2) , markerSize .* numStudies(I),'filled','Linewidth',0.5); % same as nanmean(grF(I,:),2) - nanmean(grA(I,:),2); 
    s3.MarkerFaceColor = CT1(end-4,:); s3.MarkerEdgeColor = CT4(end,:); 

% legend
offAxis = 1000; 
p1 = scatter(ax, offAxis, offAxis, 1, CT4(end-2,:), 'filled');          
	p1.Marker = 'o'; p1.MarkerEdgeColor = CT4(end,:); p1.MarkerFaceAlpha = 1.0; p1.LineWidth = 0.5; 
p2 = scatter(ax, offAxis, offAxis, 10, CT4(end-2,:), 'filled');          
	p2.Marker = 'o'; p2.MarkerEdgeColor = CT4(end,:); p2.MarkerFaceAlpha = 1.0; p2.LineWidth = 0.5; 
p3 = scatter(ax, offAxis, offAxis, 100, CT4(end-2,:), 'filled'); p2.Marker = 'o'; % 'square';         
	p3.Marker = 'o'; p3.MarkerEdgeColor = CT4(end,:); p3.MarkerFaceAlpha = 1.0; p3.LineWidth = 0.5; 
[leg,icons] = legend([p1 p2 p3], '1 publication', '10 publications', '100 publications',...
    'Orientation','vertical','Fontsize',9,'Box','off','NumColumns',2);
    M = findobj(icons,'type','patch');
        set(M(1),'MarkerSize',sqrt(1));
        set(M(2),'MarkerSize',sqrt(10));
        set(M(3),'MarkerSize',sqrt(100));
    leg.Position = [0.1775    0.765    0.1896    0.0789]; 

ax = axes('position', [0.1525 0.595 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[{'Women'},{' '},{'-represented'},{'as first author'}],'Fontsize',9,'Fontweight','normal',...
        'color',CT4(end-2,:),'HorizontalAlignment','center','FontAngle','italic','Rotation',90);
ax = axes('position', [0.1525 0.595 0.02 0.02]); hold on; box off; set(gca, 'visible', 'off')
    t1 = text(0,0,[{' '},{'under | over'},{' '},{' '}],'Fontsize',9,'Fontweight','normal',...
        'color',CT4(end-2,:),'HorizontalAlignment','center','FontAngle','normal','Fontweight','bold','Rotation',90);
annotation(gcf,'arrow',[0.145 0.145],[0.5225 0.4875],'Color',CT4(end-2,:).*[1 1 1]);
annotation(gcf,'arrow',[0.145 0.145],[0.6705 0.7055],'Color',CT4(end-2,:).*[1 1 1]);

cd(homeFold); 

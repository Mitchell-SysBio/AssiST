function [wellCenters, wellNames, h,h1,phenotypes,phenotypesReclass, reClassWell, reClass, drugMGC, drugMIC, drugs, pClass]=inferPlatePhenotypes(img, curString, platemapFileName, net, tfFlip, tfRotate, numNoGrowth,numRestricted,numFull, pClassThreshold)
%% classifyPlatePhynotypesMICMGC.m — Derive MIC/MGC phenotypes from per-well predictions
% Summary:
%   Converts per-well class probabilities/predictions into drug response phenotypes,
%   MIC/MGC indices, and visualization layers for a plate.
%
% Requirements:
%   Image Processing Toolbox (for visualization).
%
% Dependencies:
%   inputDialog2.m (optional reclass).
%
% Inputs:
%   img         double  Plate image
%   curString   string  File Name
%   drugs       vector(string)      drug names 
%   drugMap     matrix(string)      Drug names for each well
%   concMap     matrix(double)      Drug concentrations for each well 
%   typeMap     matrix(string)      Drug type for each well
%   net         DAG/Series          Trained network for classify()
%   tfFlip/ tfRotate                rotate/flip flags
%   numNoGrowth/numRestricted/numFull  number corresponding to each class for each growth type
%   pClassThreshold    number corresponding to confidence threshold for probability of class inference
%
% Outputs:
%  drugMGCMaxConc, drugMIC matrix      MGC & MIC call for each drug
%  phenotypes              matrix      Encoded per-well phenotypes.
%  phenotypesReclass       matrix      Post-manual reclassification phenotypes.
%  pClass                  matrix      posterior probabilities of each class per image
% 
% Output Images:
%  PhenotypeImages: phenotype color overlaid on plate after neural network classification
%  and user reclassification
%  QAImages: MIC and MGC overlaid on PhenotypeImages 
%  MetricImages: MIC and MGC overlaid on to plate with drug concentrations 
%
%
% Notes:
%   - Always coerce images to [50 50 1] single before classify to match the net.
%   - Safer manual reclassification: validate well IDs and class ranges; ignore bad tokens.
%   - Direction inference seeded from first non-control drug; deterministic fallback if none.
%% house keeping
n_NN_input_layer = 50; % size of input layer for neural network classifier (should match the train_NN.m)
allPhenotypes = net.Layers(end).Classes; % a list of all phenotypes the classifier knows
numClasses = numel(allPhenotypes);

% No Growth = black to red tones
temp = hot(numel(numNoGrowth)+3);
myColorMap =[zeros(1,3); temp(1:end-4,:)];

% Restricted Growth = green tones
temp = flipud(summer(numel(numRestricted)+2));
myColorMap = [myColorMap; temp(3:end,:)];

% Full Growth = blue tones
temp = sky(numel(numFull)+2); 
myColorMap = [myColorMap; temp(3:end,:)];

myAlphaMap = 0.3;

idRows = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
idColumns = {'1','2','3','4','5','6','7','8','9','10','11','12'};
wellNames = cell(96,1);
wellCenters = nan(96,2);
MICcolor = [137,81,41]./255;
MGCcolor = [0 0 0];


%% Display images

% adjust image (rotate and flip)
if(tfFlip), img = fliplr(img); end
if(tfRotate), img = imrotate(img,90); end

img2 = imadjust(img);
resizeFactor = 1000/min(size(img2)); % adjust the image so the short axis is 1000 pixels
img2 = imresize(img2,resizeFactor);

% Display the image
h = figure('color','w');
movegui(h,"north")
imshow(img2); hold on;
title(sprintf("%s - Mark the middle of well A1, THEN mark middle of well H12",curString), 'Interpreter','none');


% Ask user to select the upper left and bottom right wells
[x, y] = ginput(2); % Get 2 points from the user

% Define the grid dimensions for a 96-well plate (8 rows, 12 columns)
rows = 8;
cols = 12;

% Calculate the distances between wells
rowDist = (y(2) - y(1)) / (rows - 1);
colDist = (x(2) - x(1)) / (cols - 1);

% Determine well Positions 
radius = round(0.325*mean([rowDist,colDist]));
from_i_to_rc = nan(96,2);
wellNames = cell(96,1);

for row = 0:rows-1
    for col = 0:cols-1
        i = row*12+col+1;
        wellNames{i} = [idRows{row+1} idColumns{col+1}];
        wellCenters(i,1) = x(1) + col * colDist;
        wellCenters(i,2) = y(1) + row * rowDist;
        wellCenters(i,3) = row+1;
        wellCenters(i,4) = col+1;

        from_i_to_rc(i,[1 2]) = [row+1,col+1]; % record conversion from i index to row and col index

        rect = [wellCenters(i,1)-radius, wellCenters(i,2)-radius, 2*radius, 2*radius];
    end
end

%% classifing individual wells
l = radius*2;
mask =strel('disk',radius,0);
phenotypes = categorical([]); % place holder for phenotypes
pClass = nan(96,numClasses);
for i = 1:96
    rect = [wellCenters(i,1)-radius, wellCenters(i,2)-radius, l, l];
    curWellImage= imcrop(img2, rect);
    curWellImage(~mask.Neighborhood) = 0;
    if size(curWellImage,3)==3 % change from RGB to gray
        curWellImage = rgb2gray(curWellImage); 
    end
    img = im2single(imresize(curWellImage,[50 50]));
    img = reshape(img,[n_NN_input_layer n_NN_input_layer 1]); % resize image to match proportions of NN classifier
    [phenotypes(i),pClass(i,:)] = classify(net, img);
end

radius = radius*1.25;
for i=1:length(allPhenotypes)
    curPhenotype = allPhenotypes(i);
    tfFitPhenotype = (phenotypes == curPhenotype); % Logical indexing
    for j=1:length(tfFitPhenotype)
        if(tfFitPhenotype(j))
            rect = [wellCenters(j,1)-radius, wellCenters(j,2)-radius, 2*radius, 2*radius];
            rectangle('Position', rect,'FaceColor', myColorMap(i,:),'FaceAlpha', myAlphaMap,'Curvature', [1, 1]);
            x = [rect(1), rect(1), rect(1) + rect(3), rect(1) + rect(3)];
            y = [rect(2), rect(2) + rect(4), rect(2) + rect(4), rect(2)];
           
            % Determine if classification is below the threshold
            if isnan(pClassThreshold)
                text(x(1),y(1),num2str(i),'BackgroundColor','w','Color','k');
            elseif max(pClass(j,:))<pClassThreshold
                text(x(1),y(1),num2str(i),'BackgroundColor','w','Color','r');
            else
                text(x(1),y(1),num2str(i),'BackgroundColor','w','Color','k');
            end 
        end
    end
end
%% reclassify based on user input
[reClassWell, reClass] = inputDialog2(myColorMap, net.Layers(end).Classes);
phenotypesReclass = phenotypes;
if ~isempty(reClass)
    newClass = reClass;
    wellChanged = reClassWell;
    for reClassed = 1:length(newClass)
        idx = strcmp(wellNames, wellChanged{reClassed});
        phenotypesReclass(idx) = allPhenotypes(newClass(reClassed));
     end
    % changing colormaps to match reclassed
    rectAll = findall(h,'Type', 'Rectangle');
    set(rectAll, 'Visible', 'off'); % removing patch color
    for i=1:length(allPhenotypes)
        curPhenotype = allPhenotypes(i);
        tfFitPhenotype = (phenotypesReclass == curPhenotype); % Logical indexing
        for j=1:length(tfFitPhenotype)
            if(tfFitPhenotype(j))
                rect = [wellCenters(j,1)-radius, wellCenters(j,2)-radius, 2*radius, 2*radius];
                rectangle('Position', rect,'FaceColor', myColorMap(i,:),'FaceAlpha', myAlphaMap,'Curvature', [1, 1]);

                x = [rect(1), rect(1), rect(1) + rect(3), rect(1) + rect(3)];
                y = [rect(2), rect(2) + rect(4), rect(2) + rect(4), rect(2)];
                
                text(x(1),y(1),num2str(i),'BackgroundColor','w','Color','k','EdgeColor','k');
            end
        end
    end
end

%% Saving Phenotype Image
saveas(h,['./PhenotypeImages/' curString ' - PT.png']);

%% Stopping if MIC and MGC classification is unneeded
if strcmp(platemapFileName, 'none')
    drugMGC = [];
    drugMIC = []; 
    drugs = [];
    h1=[];
    return % exit out of function 
end 

%% Loading information for classifying MIC and MGC
% load plate map data and identify plate images for analysis
drugMap = categorical(readcell(platemapFileName,'Range','B3:M10'));
concMap = readmatrix(platemapFileName,'Range','B14:M21');
typeMap = readcell(platemapFileName,'Range','B25:M32');

% prepare place holder for drug sensitivity data
drugs = unique(drugMap); nDrugs = numel(drugs);

%% Classifying MIC and MGC
% generate clickmap
% RULES: full growth was detected = 6_cloud, 7_pellet; Restricted growth was detected = 4_cloud_small, 5_pellet_small
% For Restricted growth MGC will be between lower concentration and 1 step below 
% If two Restricted growth then chose lower concentration
% No growth:  1_bubbles, 2_empty, 3_pellet_tiny
phenotype96 = reshape(phenotypesReclass, 12,8)';
clickMap = cellfun(@(x) str2double(x(1)),string(phenotype96)); % making clickmap by making phenotypes into matrix in plate format
drugMGC = nan(nDrugs,1); wellIdxMGC = nan(nDrugs,1);
drugMIC = nan(nDrugs,1); wellIdxMIC = nan(nDrugs,1);
MICplots = cell(nDrugs,2); MGCplots = cell(nDrugs,2); 
for iDrug = 1:nDrugs
    curDrug = drugs(iDrug);
    inxDrugUS = find(drugMap==curDrug);
    curConc = concMap(inxDrugUS);
    stattype = typeMap{inxDrugUS(1)};
    [curConc, sortidx] = sort(curConc,'descend'); % making sure it's sorted from high to low 
    inxDrug = inxDrugUS(sortidx); % reordering inxDrug to make sure it's following high to low of concentration
    % determining which direction is from high to low concentration
    direction = "descend"; % intitializing just in case control is first
    if sum(curConc) ~= 0 % not positive control 
        [R,C] = ind2sub([8,12],inxDrug); 
        % --- must be one row or one column
        isRow = (max(R)-min(R))==0 & (max(C)-min(C))>0;
        isCol = (max(C)-min(C))==0 & (max(R)-min(R))>0;
    
        % --- order along the line (left→right for a row, top→bottom for a column)
        if isRow
            [x, ord] = sort(C);                 % left (small C) -> right (large C)
            dir_words = ["right","left"]; % "left → right","right → left"
        else
            [x, ord] = sort(R);                 % top (small R) -> bottom (large R)
            dir_words = ["descend","ascend"];%["top → bottom","bottom → top"]
        end
        y  = curConc(ord);
        d  = diff(y);                           % stepwise differences
        inc = all(d > 0);
        dec = all(d < 0);
        if inc
            direction = dir_words(1);
        elseif dec
            direction = dir_words(2);
        end 
    end 
    curDrugDilutionFactor = curConc(1)/curConc(2); % calculating the fold conc steps
    curClickValues = clickMap(inxDrug); % order is high to low concentration of drug
    
    % making curClickValues into no growth, Restricted growth, full growth
    for ng = numNoGrowth
        % making no growth = 0 
        idx = find(curClickValues == ng);
        curClickValues(idx) = 0; 
    end 
    for mg = numRestricted
        % making Restricted growth = 1
        idx = find(curClickValues == mg);
        curClickValues(idx) = 1; 
    end 
    for fg = numFull
        % making full growth = 2
        idx = find(curClickValues == fg);
        curClickValues(idx) = 2; 
    end 
    % curClickValues(curClickValues==6) = 2; % full growth was detected = 6_cloud
    % curClickValues(curClickValues==7) = 2; % full growth was detected = 7_pellet
    % curClickValues(curClickValues==4) = 1; % Restricted growth was detected = 4_cloud_small
    % curClickValues(curClickValues==5) = 1; % Restricted growth was detected = 5_pellet_small
    growthInd = find(curClickValues>0); % finding indexes of all that have growth

    if ~isempty(growthInd)
        gIndex = 1;
        maxClickValue = curClickValues(growthInd(gIndex));

        % handing skip well case
        if ~all(curClickValues(growthInd(gIndex):end)) % if there is a no growth well in lower concentrations
            if length(growthInd) == 1 % skip well and all lower concentrations is no growth
                % case #3 - no growth reported - will report as next lowest concentration using dilution factor
                maxClickValue = 0;
            else % skip well and there is growth in lower concentrations
                gIndex = 2;
                maxClickValue = curClickValues(growthInd(gIndex)); % next lowest concentration with growth is used
            end
        end

        % assigning concentration values
        % MGC: concentration at which there is maximum growth
        % MIC: minimum concentration at which there is no growth
        if maxClickValue == 2 % case #1 - full growth was detected = 7_cloud, 8_pellet
            curMaxConc = curConc(growthInd(gIndex));
            idxMGC = inxDrug(growthInd(gIndex));% getting index of well with MGC concentration
            % getting index of well with MIC concentration 
            if growthInd(gIndex)-1 == 0
                % MIC is higher than tested concentration
                idxMIC = 0; 
                curMICmaxConc = max(curConc)*curDrugDilutionFactor; % MIC is saved as 1 step up
            else
                idxMIC = inxDrug(growthInd(gIndex)-1);
                curMICmaxConc = curConc(growthInd(gIndex)-1);
            end 
        elseif maxClickValue == 1 % case #2 = Restricted growth was detected = 5_cloud_small, 6_pellet_small
            mgIDX = curClickValues(growthInd)==1; % getting Restricted growth idx 
            firstZ = find(mgIDX==0, 1); % geting first full growth well 
            if length(growthInd(gIndex:end))>1 && ~isempty(firstZ)
                % getting MGC concentration 
                % if there are multiple Restricted growth, then choose lowest
                % concentration and do average of that concentration and 1
                % step down
                curMaxConcRestricted = curConc(growthInd(firstZ-1)); % getting lowest conc
                curMaxConcRestrictedStepDown = curMaxConcRestricted/curDrugDilutionFactor;
                curMaxConc = mean([curMaxConcRestricted curMaxConcRestrictedStepDown]);
                % getting indexes of well with MGC concentration
                idxMGC = inxDrug(growthInd(firstZ-1))+0.5;

                % getting MIC concentration 
                % If bacteriostatic drug then use 1 concentration above
                % full growth (ignoring tiny and small pellets) 
                if strcmp(stattype, 'bacteriostatic') 
                    if growthInd(gIndex)-1 == 0
                        idxMIC = 0; % MIC is higher than tested concentration
                        curMICmaxConc = max(curConc)*curDrugDilutionFactor;
                    else
                        if sum(curClickValues(growthInd)==2) == 0 
                            % no full growths then MIC is lowest tested
                            % concentration 
                            curMICmaxConc = min(curConc);
                            % getting index of well with MIC concentration
                            idxMIC = inxDrug(min(curConc) == curConc);
                        else 
                            % find full growth well at highest conc and get concentration 1
                            % higher
                            fidx = min(find(curClickValues(growthInd)==2)); 
                            idxMIC = inxDrug(growthInd(fidx)-1);
                            curMICmaxConc = curConc(growthInd(fidx)-1);
                        end 
                    end
                else 
                    % Bactericidal so Restricted growth counts 
                    if growthInd(gIndex)-1 == 0
                        idxMIC = 0; % MIC is higher than tested concentration
                        curMICmaxConc = max(curConc)*curDrugDilutionFactor;
                    else
                        idxMIC = inxDrug(growthInd(gIndex)-1);
                        curMICmaxConc = curConc(growthInd(gIndex)-1);
                    end
                end 
            else % only 1 Restricted growth
                curMaxConcRestricted = curConc(growthInd(gIndex));
                curMaxConcRestrictedStepDown = curMaxConcRestricted/curDrugDilutionFactor;
                curMaxConc = mean([curMaxConcRestricted curMaxConcRestrictedStepDown]);
                % getting index of well with MGC concentration
                idxMGC = inxDrug(growthInd(gIndex))+0.5;

                % getting MIC concentration 
                % If bacteriostatic drug then use 1 concentration above
                % full growth (ignoring tiny and small pellets) 
                if strcmp(stattype, 'bacteriostatic') 
                    if growthInd(gIndex)-1 == 0
                        idxMIC = 0; % MIC is higher than tested concentration
                        curMICmaxConc = max(curConc)*curDrugDilutionFactor;
                    else
                        if sum(curClickValues(growthInd)==2) == 0 
                            % no full growths then MIC is lowest tested
                            % concentration 
                            curMICmaxConc = min(curConc);
                            % getting index of well with MIC concentration
                            idxMIC = inxDrug(min(curConc) == curConc);
                        else 
                            % find full growth well  at highest conc and get concentration 1
                            % higher
                            fidx = min(find(curClickValues(growthInd)==2)); 
                            idxMIC = inxDrug(growthInd(fidx)-1);
                            curMICmaxConc = curConc(growthInd(fidx)-1);
                        end 
                    end
                else 
                    % Bactericidal so Restricted growth counts 
                    if growthInd(gIndex)-1 == 0
                        idxMIC = 0; % MIC is higher than tested concentration
                        curMICmaxConc = max(curConc)*curDrugDilutionFactor;
                    else
                        idxMIC = inxDrug(growthInd(gIndex)-1);
                        curMICmaxConc = curConc(growthInd(gIndex)-1);
                    end
                end 

            end
        else % skip well and no other growth
            % case #3 - no growth reported - will report as next lowest concentration using dilution factor
            curMaxConc = min(curConc)/curDrugDilutionFactor;
            % getting index of well with MGC concentration
            idxMGC =inxDrug(min(curConc) == curConc)+0.5;

            curMICmaxConc = min(curConc);
            % getting index of well with MIC concentration 
            idxMIC = inxDrug(min(curConc) == curConc);
        end

    else % case #3 - no growth reported - will report as next lowest concentration using dilution factor
        curMaxConc = min(curConc)/curDrugDilutionFactor;
        % getting index of well with MGC concentration
        idxMGC =inxDrug(min(curConc) == curConc)+0.5;

        curMICmaxConc = min(curConc);
        % getting index of well with MIC concentration
        idxMIC = inxDrug(min(curConc) == curConc);
    end

    % MGC: storing concentration with maximum growth
    drugMGC(iDrug) = curMaxConc;
    wellIdxMGC(iDrug) = idxMGC;

    % MIC: storing 1st concentration with no growth
    drugMIC(iDrug)= curMICmaxConc;
    wellIdxMIC(iDrug)= idxMIC;
   
    % Labeling MIC circle 
    if idxMIC ==0
        % if MIC is higher than test concentrations
        [r,c] = ind2sub([8,12], inxDrug(1)); 
        aidxMIC = find(wellCenters(:,3)==r & wellCenters(:,4)==c); % finding right index for wellCenters
            if strcmp(direction, "descend")
                x = linspace(wellCenters(aidxMIC,1)-radius, wellCenters(aidxMIC,1)+radius,10);
                y = zeros(1,10)+ wellCenters(aidxMIC,2)+radius;
            elseif strcmp(direction, "ascend")
                x = linspace(wellCenters(aidxMIC,1)-radius, wellCenters(aidxMIC,1)+radius,10);
                y = zeros(1,10)+ wellCenters(aidxMIC,2)-radius;
            elseif strcmp(direction, "left")
                x = zeros(1,10)+ wellCenters(aidxMIC,1)-radius;
                y = linspace(wellCenters(aidxMIC,2)-radius, wellCenters(aidxMIC,2)+radius,10); 
            elseif strcmp(direction, "right")
                x = zeros(1,10)+ wellCenters(aidxMIC,1)+radius;
                y = linspace(wellCenters(aidxMIC,2)-radius, wellCenters(aidxMIC,2)+radius,10); 
            end
            plot(x,y, 'Color', MICcolor,'LineWidth',4)
            MICplots{iDrug,2} = [x;y];
            MICplots{iDrug,1} = 'line';
    else
        [r,c] = ind2sub([8,12], idxMIC); 
        aidxMIC = find(wellCenters(:,3)==r & wellCenters(:,4)==c); % finding right index for wellCenters
        rect = [wellCenters(aidxMIC,1)-radius, wellCenters(aidxMIC,2)-radius, 2*radius, 2*radius];
        MICplots{iDrug,2} = rect; 
        MICplots{iDrug,1} = 'circle';
        rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', MICcolor, 'LineWidth', 4,'LineStyle','-');
    end
    
    % Label MGC circle
    if mod(idxMGC,1)==0.5
        % If Restricted growth 
         [r,c] = ind2sub([8,12], idxMGC-0.5); 
        aidxMGC = find(wellCenters(:,3)==r & wellCenters(:,4)==c); % finding right index for wellCenters
            if strcmp(direction, "descend")
                x = linspace(wellCenters(aidxMGC,1)-radius, wellCenters(aidxMGC,1)+radius,10);
                y = zeros(1,10)+ wellCenters(aidxMGC,2)-radius;
             elseif strcmp(direction, "ascend")
                x = linspace(wellCenters(aidxMGC,1)-radius, wellCenters(aidxMGC,1)+radius,10);
                y = zeros(1,10)+ wellCenters(aidxMGC,2)+radius; 
            elseif strcmp(direction, "left")
                x = zeros(1,10)+ wellCenters(aidxMGC,1)+radius;
                y = linspace(wellCenters(aidxMGC,2)-radius, wellCenters(aidxMGC,2)+radius,10);
            elseif strcmp(direction, "right")
                x = zeros(1,10)+ wellCenters(aidxMGC,1)-radius;
                y = linspace(wellCenters(aidxMGC,2)-radius, wellCenters(aidxMGC,2)+radius,10); 
            end
            plot(x, y, 'Color', MGCcolor,'LineWidth',4)
            MGCplots{iDrug,2} = [x;y];
            MGCplots{iDrug,1} = 'line';

    else
        [r,c] = ind2sub([8,12], idxMGC); 
        aidxMGC = find(wellCenters(:,3)==r & wellCenters(:,4)==c); % finding right index for wellCenters
        rect = [wellCenters(aidxMGC,1)-radius, wellCenters(aidxMGC,2)-radius, 2*radius, 2*radius];
        MGCplots{iDrug,2} = rect; 
        MGCplots{iDrug,1} = 'circle';
        rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', MGCcolor, 'LineWidth', 4,'LineStyle','-');
    end 
end
title(sprintf("%s: MIC = Brown; MGC = Black",curString));
%% saving image
saveas(h,['./QAimages/' curString ' - QA.png']);
%% Plotting MIC and MGC onto drug concentrations as image
rowSign = ["A", "B", "C", "D", "E","F","G","H"];
colSign = ["01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"];
% Display the image
h1 = figure('color','w');
movegui(h1,"north")
imshow(ones(size(img2)))
hold on;
title(sprintf("%s: MIC = Brown; MGC = Black",curString));
for c1=1:12
    text(wellCenters(c1,1), min(wellCenters(:,2))-1.5*radius, colSign(c1), 'HorizontalAlignment','center')
end 
uniWC = unique(wellCenters(:,2));
for r1 = 1:8
    text(min(wellCenters(:,1))-1.5*radius,uniWC(r1), rowSign(r1), 'HorizontalAlignment','right')
end 

for r = 1:8
for c = 1:12
    wc = find(wellCenters(:,3)==r & wellCenters(:,4)==c); % finding right index for wellCenters
    rect = [wellCenters(wc,1)-radius, wellCenters(wc,2)-radius, 2*radius, 2*radius];
    rectangle('Position', rect,'FaceColor', [ 1 1 1], 'Curvature', [1, 1]);
    text(wellCenters(wc,1), wellCenters(wc,2),num2str(concMap(r,c)), 'HorizontalAlignment', 'center')  
end 
end 

for i = 1:length(MICplots)
 % Label MIC 
    if strcmp(MICplots{i,1},"line")
        plot(MICplots{i,2}(1,:),MICplots{i,2}(2,:), 'Color', MICcolor,'LineWidth',4)
    else
        rectangle('Position', MICplots{i,2},'Curvature', [1, 1], 'EdgeColor', MICcolor, 'LineWidth', 4,'LineStyle','-');
    end

    % Label MGC 
    if strcmp(MGCplots{i,1},"line")
        plot(MGCplots{i,2}(1,:),MGCplots{i,2}(2,:), 'Color', MGCcolor,'LineWidth',4)
    else
        rectangle('Position', MGCplots{i,2},'Curvature', [1, 1], 'EdgeColor', MGCcolor, 'LineWidth', 4,'LineStyle','-');
    end
end 
%% saving image
saveas(h1,['./MetricImages/' curString ' - M.png']);
end 

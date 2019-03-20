%% Add your names in a comment here at the beginning of the code!

%Jocelyn and Leafia aka Shark and Leaf

% Instructions: Follow through this code step by step, while also referring
% to the overall instructions and questions from the lab assignment sheet.

%% 1. Read in the monthly gridded CO2 data from the .csv file
% The data file is included in your repository as “LDEO_GriddedCO2_month_flux_2006c.csv”
% Your task is to write code to read this in to MATLAB
% Hint: you can again use the function “readtable”, and use your first data lab code as an example.
%<--
CO2data = readtable('LDEO_GriddedCO2_month_flux_2006c.csv');

%% 2a. Create new 3-dimensional arrays to hold reshaped data
%Find each unique longitude, latitude, and month value that will define
%your 3-dimensional grid
longrid = unique(CO2data.LON); %finds all unique longitude values
latgrid = unique(CO2data.LAT); %<-- following the same approach, find all unique latitude values
monthgrid = unique(CO2data.MONTH); %<-- following the same approach, find all unique months

%Create empty 3-dimensional arrays of NaN values to hold your reshaped data
    %You can make these for any variables you want to extract - for this
    %lab you will need PCO2_SW (seawater pCO2) and SST (sea surface
    %temperature)
PCO2_SW = NaN*zeros(length(longrid),length(latgrid),length(monthgrid));%<--
SST = NaN*zeros(length(longrid),length(latgrid),length(monthgrid));%<--

%% 2b. Pull out the seawater pCO2 (PCO2_SW) and sea surface temperature (SST)
%data and reshape it into your new 3-dimensional arrays

for i= 1:length(CO2data.LON)
           latindex = find(latgrid==CO2data.LAT(i));
           lonindex = find(longrid==CO2data.LON(i));
           PCO2_SW(lonindex,latindex,CO2data.MONTH(i))=CO2data.PCO2_SW(i);
           SST(lonindex,latindex,CO2data.MONTH(i))=CO2data.SST(i);
end%<--


%% 3a. Make a quick plot to check that your reshaped data looks reasonable
%Use the imagesc plotting function, which will show a different color for
%each grid cell in your map. Since you can't plot all months at once, you
%will have to pick one at a time to check - i.e. this example is just for
%January

imagesc(SST(:,:,1))

%% 3b. Now pretty global maps of one month of each of SST and pCO2 data.
%I have provided example code for plotting January sea surface temperature
%(though you may need to make modifications based on differences in how you
%set up or named your variables above).

figure(1); clf
worldmap world
contourfm(latgrid, longrid, SST(:,:,1)','linecolor','none');
c = colorbar('southoutside'); 
c.Label.String = '\it Sea Surface Temperature [^oC]';
geoshow('landareas.shp','FaceColor','black');
%title('January Sea Surface Temperature')

%%
%Check that you can make a similar type of global map for another month
%and/or for pCO2 using this approach. Check the documentation and see
%whether you can modify features of this map such as the contouring
%interval, color of the contour lines, labels, etc.

figure(2); clf
worldmap world
contourfm(latgrid, longrid, SST(:,:,2)','linecolor','none');
c = colorbar('southoutside');
c.Label.String = '\it Sea Surface Temperature [^oC]';
geoshow('landareas.shp','FaceColor','black');
title('February Sea Surface Temperature');

%%
figure(3); clf
worldmap world
contourfm(latgrid, longrid, PCO2_SW(:,:,1)','linecolor','none');
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('January pCO_2');

%% 4. Calculate and plot a global map of annual mean pCO2

PCO2_mean = mean(PCO2_SW,3); %Calculate mean pCO2 along the 3rd dimension (over the months)

figure(4); clf
worldmap world
contourfm(latgrid, longrid, PCO2_mean','linecolor','none');
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('Annual Mean Seawater pCO_2');

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2

%reference year: 2000
%365ppm 
%data obtained from Full Record Graph of the Keeling Curve from scripps
%institution of oceanography, UC San Diego --- https://scripps.ucsd.edu/programs/keelingcurve/
%selected data source because we noticed that it was highly cited in online
%material targeted towards general audiences

refyear = 365; %CO2 ppm of reference year
difference= PCO2_mean - refyear; %difference between mean PCO2 and the reference CO2

figure(5); clf
worldmap world
contourfm(latgrid, longrid, difference','linecolor','none');
cmocean('balance', 'pivot', 0);
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('Difference between Annual Mean pCO2 and 2000 Atmospheric pCO2');

%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle
%<--
SST_mean = mean(SST,3); %Calculate the mean SST along the 3rd dimension (over all the months)
PCO2_BP = PCO2_SW.*exp(0.0423.*(repmat(SST_mean,1,1,12)-SST));
PCO2_T = PCO2_mean.*exp(0.0423.*(SST-repmat(SST_mean,1,1,12)));

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

BATSlat=32; %32.833;   
BATSlon=297.5; %295.833;
Rosslat=-76; %-76.83;   
Rosslon=177.5; %176 is actual lon
Papalat=48; %50 is the actual lat but since lat is in steps of 4 we had to choose the closest cell further from shore;       
Papalon=212.5;%215; is the actual lon but its in steps of 5 starting at 2.5

ind_BATS = find(CO2data.LAT==BATSlat&CO2data.LON==BATSlon);
ind_Ross = find(CO2data.LAT==Rosslat&CO2data.LON==Rosslon);
ind_Papa = find(CO2data.LAT==Papalat&CO2data.LON==Papalon);

PCO2_extracted = NaN*zeros(12,3);
SST_extracted = NaN*zeros(12,3);
PCO2_BP_extracted = NaN*zeros(12,3);
PCO2_T_extracted = NaN*zeros(12,3);

for i=1:length(ind_BATS)
    PCO2_extracted(i,1) = CO2data.PCO2_SW(ind_BATS(i));
    SST_extracted(i,1) = CO2data.SST(ind_BATS(i));
end

for i=1:length(ind_Ross)
    PCO2_extracted(i,2) = CO2data.PCO2_SW(ind_Ross(i));
    SST_extracted(i,2) = CO2data.SST(ind_Ross(i));
end

for i=1:length(ind_Papa)
    PCO2_extracted(i,3) = CO2data.PCO2_SW(ind_Papa(i));
    SST_extracted(i,3) = CO2data.SST(ind_Papa(i));
end

MeanPCO2_BATS = mean(PCO2_extracted(:,1));
MeanPCO2_Ross = mean(PCO2_extracted(:,2));
MeanPCO2_Papa = mean(PCO2_extracted(:,3));

MeanSST_BATS = mean(SST_extracted(:,1));
MeanSST_Ross = mean(SST_extracted(:,2));
MeanSST_Papa = mean(SST_extracted(:,3));

for i=1:length(ind_BATS)
    PCO2_BP_extracted(i,1) = PCO2_extracted(i,1)*exp(0.0423*(MeanSST_BATS-SST_extracted(i,1)));
    PCO2_T_extracted(i,1) = MeanPCO2_BATS*exp(0.0423*(SST_extracted(i,1)-MeanSST_BATS));
end

for i=1:length(ind_Ross)
    PCO2_BP_extracted(i,2) = PCO2_extracted(i,2)*exp(0.0423*(MeanSST_Ross-SST_extracted(i,2)));
    PCO2_T_extracted(i,2) = MeanPCO2_Ross*exp(0.0423*(SST_extracted(i,2)-MeanSST_Ross));
end

for i=1:length(ind_Papa)
    PCO2_BP_extracted(i,3) = PCO2_extracted(i,3)*exp(0.0423*(MeanSST_Papa-SST_extracted(i,3)));
    PCO2_T_extracted(i,3) = MeanPCO2_Papa*exp(0.0423*(SST_extracted(i,3)-MeanSST_Papa));
end
%% Create Figure for BATS Site
figure(6); clf
subplot(2,1,1);
hold on
title('Observed pCO_2 and Sea Surface Temperature at BATS Site')
yyaxis left
plot(monthgrid, PCO2_extracted(:,1), '-o', 'LineWidth', 1)
ylabel('Observed pCO_2 [µatm]')
yyaxis right
plot(monthgrid, SST_extracted(:,1), '-or', 'LineWidth', 1)
ylabel('Sea Surface Temperature [^oC]')
xlabel('Months')
legend({'Observed pCO_2', 'Sea Surface Temperature'}, 'Location', 'northwest')
legend('boxoff')

subplot(2,1,2);
hold on
plot(monthgrid,PCO2_extracted(:,1), '-ok', 'LineWidth', 1)
plot(monthgrid,PCO2_BP_extracted(:,1), '-o', 'LineWidth', 1) 
plot(monthgrid,PCO2_T_extracted(:,1), '-or', 'LineWidth', 1)
xlabel('Months')
ylabel('pCO_2 [µatm]')
legend({'Observed pCO_2','Seasonal Biophysical Effect', 'Seasonal Temperature Effect'}, 'Location', 'northwest') %label stuff + change the legend's location
legend('boxoff') %remove the box from the legend
title('Seasonal Variability in pCO_2 at BATS Site')
%% Create Figure for Ross Sea
figure(7); clf
subplot(2,1,1);
hold on
title('Observed pCO_2 and Sea Surface Temperature at Ross Sea')
yyaxis left
plot(monthgrid, PCO2_extracted(:,2), '-o', 'LineWidth', 1)
ylabel('Observed pCO_2 [µatm]')
yyaxis right
plot(monthgrid, SST_extracted(:,2), '-or', 'LineWidth', 1)
ylabel('Sea Surface Temperature [^oC]')
xlabel('Months')
legend({'Observed pCO_2', 'Sea Surface Temperature'}, 'Location', 'southeast')
legend('boxoff')

subplot(2,1,2);
hold on
plot(monthgrid,PCO2_extracted(:,2), '-ok', 'LineWidth', 1)
plot(monthgrid,PCO2_BP_extracted(:,2), '-o', 'LineWidth', 1) 
plot(monthgrid,PCO2_T_extracted(:,2), '-or', 'LineWidth', 1)
xlabel('Months')
ylabel('pCO_2 [µatm]')
legend({'Observed pCO_2','Seasonal Biophysical Effect', 'Seasonal Temperature Effect'}, 'Location', 'southeast') %label stuff + change the legend's location
legend('boxoff') %remove the box from the legend
title('Seasonal Variability in pCO_2 at Ross Sea')
%% Create Figure for Weather Station Papa
figure(8); clf
subplot(2,1,1);
hold on
title('Observed pCO_2 and Sea Surface Temperature at Weather Station Papa')
yyaxis left
plot(monthgrid, PCO2_extracted(:,3), '-o', 'LineWidth', 1)
ylabel('Observed pCO_2 [µatm]')
yyaxis right
plot(monthgrid, SST_extracted(:,3), '-or', 'LineWidth', 1)
ylabel('Sea Surface Temperature [^oC]')
xlabel('Months')
legend({'Observed pCO_2', 'Sea Surface Temperature'}, 'Location', 'northwest')
legend('boxoff')

subplot(2,1,2);
hold on
plot(monthgrid,PCO2_extracted(:,3), '-ok', 'LineWidth', 1)
plot(monthgrid,PCO2_BP_extracted(:,3), '-o', 'LineWidth', 1) 
plot(monthgrid,PCO2_T_extracted(:,3), '-or', 'LineWidth', 1)
xlabel('Months')
ylabel('pCO_2 [µatm]')
legend({'Observed pCO_2','Seasonal Biophysical Effect', 'Seasonal Temperature Effect'}, 'Location', 'northwest') %label stuff + change the legend's location
legend('boxoff') %remove the box from the legend
title('Seasonal Variability in pCO_2 at Weather Station Papa')

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above Reproduce Figure 6! ***extension option**
SeasonID = SST - repmat(SST_mean,1,1,12);%calculate the difference between the observed SST and the annual mean SST for each grid cell
diff_PCO2 = NaN*zeros(size(SeasonID,1),size(SeasonID,2)); %create an empty 72x40 array to store the seasonal differences in pCO2
for i = 1:size(SeasonID,1)
    for j = 1:size(SeasonID,2)
        [maxPCO2_val,maxPCO2_idx] = max(PCO2_SW(i,j,:),[],3);%calculate the maximum PCO2 value and its index along the 3rd dimension
        [minPCO2_val,minPCO2_idx] = min(PCO2_SW(i,j,:),[],3);%calculate the minimum PCO2 value and its index along the 3rd dimension
        if SeasonID(i,j,maxPCO2_idx) > 0 %when the SST is greater than the mean, diff_PCO2 = maxPCO2 - minPCO2
            diff_PCO2(i,j) = maxPCO2_val - minPCO2_val;
        elseif SeasonID(i,j,maxPCO2_idx) < 0 %when the SST is less than the mean, diff_PCO2 = minPCO2 - maxPCO2
            diff_PCO2(i,j) = minPCO2_val - maxPCO2_val;
        end
    end
end

figure(9); clf %but actually Figure 6 in Takahashi et al. (2002) lol
hold on
worldmap world
contourfm(latgrid, longrid, diff_PCO2','linecolor','none');
scatterm(BATSlat, BATSlon, 50, 'oy', 'filled');
scatterm(Rosslat, Rosslon, 50, 'oy', 'filled');
scatterm(Papalat, Papalon, 50, 'oy', 'filled');
cmocean('balance', 'pivot', 0)
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('Seasonal Mean Monthly Amplitude of pCO_2 in Seawater');
%% Reproduce Figure 7
diff_BP = max(PCO2_BP,[], 3)-min(PCO2_BP,[], 3); %This is Equation 3 in Takahashi et al. (2002); this is the equation used in Fig.7 of that paper!
figure(10); clf %but actually Figure 7 in Takahashi et al. (2002) lol
hold on
worldmap world
contourfm(latgrid, longrid, diff_BP','linecolor','none');
scatterm(BATSlat, BATSlon, 50, 'om', 'filled');
scatterm(Rosslat, Rosslon, 50, 'om', 'filled');
scatterm(Papalat, Papalon, 50, 'om', 'filled');
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('Seasonal Biophysical Drawdown of Seawater pCO_2');
%% Reproduce Figure 8
diff_T = max(PCO2_T,[], 3)-min(PCO2_T,[], 3); %This is Equation 4 in Takahashi et al. (2002); this is the equation used in Fig.8 of that paper!
figure(11); clf %but actually Figure 8 in Takahashi et al. (2002) lol
hold on
worldmap world
contourfm(latgrid, longrid, diff_T','linecolor','none');
scatterm(BATSlat, BATSlon, 50, 'om', 'filled');
scatterm(Rosslat, Rosslon, 50, 'om', 'filled');
scatterm(Papalat, Papalon, 50, 'om', 'filled');
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('Seasonal Temperature Effect on Seawater pCO_2');
%% Reproduce Figure 9
TminusB = diff_T - diff_BP;
figure(12); clf %but actually Figure 9 in Takahashi et al. (2002) lol
hold on
worldmap world
contourfm(latgrid, longrid, TminusB','linecolor','none');
scatterm(BATSlat, BATSlon, 50, 'oy', 'filled');
scatterm(Rosslat, Rosslon, 50, 'oy', 'filled');
scatterm(Papalat, Papalon, 50, 'oy', 'filled');
cmocean('balance', 'pivot', 0)
c = colorbar('southoutside');
c.Label.String = '\it pCO_2 [µatm]';
geoshow('landareas.shp','FaceColor','black');
title('Difference Between Seasonal Temperature Effect and Biophysical Effect on pCO_2');



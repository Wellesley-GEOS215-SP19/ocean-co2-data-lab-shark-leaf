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
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January Sea Surface Temperature (^oC)')

%Check that you can make a similar type of global map for another month
%and/or for pCO2 using this approach. Check the documentation and see
%whether you can modify features of this map such as the contouring
%interval, color of the contour lines, labels, etc.

%<--
figure(2); clf
worldmap world
contourfm(latgrid, longrid, SST(:,:,2)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('February Sea Surface Temperature (^oC)')

figure(3); clf
worldmap world
contourfm(latgrid, longrid, PCO2_SW(:,:,1)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January pCO2 (µatm)')

%% 4. Calculate and plot a global map of annual mean pCO2
%<--

PCO2_mean = mean(PCO2_SW,3);

figure(4); clf
worldmap world
contourfm(latgrid, longrid, PCO2_mean','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Annual Mean pCO2 (µatm)')

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2
%<--
%reference year: 2000
%365ppm 
%data obtained from Full Record Graph of the Keeling Curve from scripps
%institution of oceanography, UC San Diego --- https://scripps.ucsd.edu/programs/keelingcurve/
%selected data source because we noticed that it was highly cited in online
%material targeted towards general audiences

refyear = 365;
difference= PCO2_mean - refyear;

figure(5); clf
worldmap world
contourfm(latgrid, longrid, difference','linecolor','none');
cmocean('balance', 'pivot', 0)
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Difference between Annual Mean pCO2 and 2000 Atmospheric pCO2 (µatm)')


%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle
%<--
SST_mean = mean(SST,3); %b/c on the third dimension
PCO2_BP = PCO2_SW.*exp(0.0423.*(repmat(SST_mean,1,1,12)-SST));
PCO2_T = PCO2_mean.*exp(0.0423.*(SST-repmat(SST_mean,1,1,12)));

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

%<--
BATSlat=32; %32.833;   
BATSlon=297.5; %295.833;
Rosslat=-76; %-76.83;   
Rosslon=177.5; %176 is actual lon
Papalat=48; %50 is the actual lat but since lat is in steps of 4 we had to choose the closest cell further from shore;       
Papalon=212.5 ;%215; is the actual lon but its in steps of 5 starting at 2.5

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

figure(6); clf
hold on
plot(monthgrid,PCO2_BP_extracted(:,1))
plot(monthgrid,PCO2_T_extracted(:,1))
plot(monthgrid,PCO2_extracted(:,1))
title('BATS')

figure(7); clf
hold on
plot(monthgrid,PCO2_BP_extracted(:,1))
plot(monthgrid,PCO2_T_extracted(:,1))
plot(monthgrid,PCO2_extracted(:,1))
title('BATS SST')

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above

%<--Fig 7 was PCO2_BP
%Fig 8 PCO2_T
%fig 9 is PCO2_T - PCO2_BP

%copy plot and use scatterm on map
%use contourfm


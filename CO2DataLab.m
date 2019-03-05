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
figure(6); clf
worldmap world
%contourfm(latgrid, longrid, difference','linecolor','none');
scatterm(32.833,295.833,100)
scatterm(-76.83,
cmocean('balance', 'pivot', 0)
colorbar
geoshow('landareas.shp','FaceColor','black')
title('location')

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above

%<--

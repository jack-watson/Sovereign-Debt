
function countryHeatmap(vals, stateNames, cmap)
% Adapted from MATLAB Central user DGM
% Country shapefile sourced from: 
% https://www.naturalearthdata.com/downloads/10m-cultural-vectors/

% country names
% stateNames = {'Belgium', 'Bulgaria', 'Czech Republic', 'Denmark', 'Germany','Estonia','Ireland',...
%     'Greece','Spain','France','Croatia','Italy', 'Cyprus','Latvia', 'Lithuania', 'Luxembourg',...
%     'Hungary', 'Malta', 'Netherlands','Austria', 'Poland', 'Portugal','Romania','Slovenia', 'Slovakia',...
%     'Finland', 'Sweden'};

% cmap = jet(256); % specify the colormap
% Random vector of values (one for each country)

if nargin < 3 % default colormap is "hot" with the gradient reversed
    cmap = flipud(hot(256));
end

numRegions = length(stateNames);
coloridx = round(rescale(vals,1,size(cmap,1))); % rescale the data to use as an index into cmap

% Map
figure('Color','w')
xh = worldmap('world');
setm(xh,'MapProjection','miller')
tightmap
geoshow ('landareas.shp', 'FaceColor', 'white');
for i = 1:1:numRegions
    bordersm(stateNames{i}, 'facecolor', cmap(coloridx(i), :))
end
colorbar
colormap(cmap) % specify the colormap
clim([min(vals), max(vals)]) % set caxis to match datarange

end
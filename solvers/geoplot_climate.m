function geoplot_climate(CurrentSubs,Beta,options)
%% Chu Hong 13 April 2021
% Plot weighted grid on map

arguments 
    CurrentSubs (2,1) double = [31;138];  % given location by indices, NOT (lon,lat), default : Darkar
    Beta (:,1) double = randn(10512,1) % weight to plot in each location, length == total points 
    options.LatitudeVector (:,1) double = (90:-2.5:-90)' % default : grid of length 73 with step 2.5
    options.LongitudeVector (:,1) double = (2.5*[0,1:143])' % default : grid of length 144 with step 2.5
    options.count_long_to_lat (1,1) double = 1;
end

%% Apprearance 
SmallPointSize = 4; ZoomPointSize = 30;
MainPointColor = [0,0,1]; 
ZoomBoxColor = [0,0.75,0];

%% Compute parameters
LatitudeVector = options.LatitudeVector;
LongitudeVector = options.LongitudeVector;

NumLat = length(LatitudeVector);
NumLong = length(LongitudeVector); 
NumPoints = NumLat*NumLong;

if options.count_long_to_lat % if data was read longitude by longitude 
    CurrentIndex = sub2ind(flip([NumLat,NumLong]),CurrentSubs(2),CurrentSubs(1));
    BetaFull = [Beta(1:CurrentIndex-1);0;Beta(CurrentIndex+1:end)];
    BetaFullScaled = BetaFull./max(BetaFull);
    BetaFullScaled = reshape(reshape(BetaFullScaled,NumLong,NumLat)',[],1); 
else % if data was read latitude by latitude
    CurrentIndex = sub2ind([NumLat,NumLong],CurrentSubs(1),CurrentSubs(2));
    BetaFull = [Beta(1:CurrentIndex-1);0;Beta(CurrentIndex+1:end)];
    BetaFullScaled = BetaFull./max(BetaFull);
end

%% Plot
% Make the current point to be center
CurrentLat = LatitudeVector(CurrentSubs(1)); CurrentLong = LongitudeVector(CurrentSubs(2));
LatLimit = [-90,90]; LongLimit = [-180,180];
LongitudeVector = LongitudeVector - 360*(LongitudeVector>CurrentLong+180);
LongitudeVector = LongitudeVector + 360*(LongitudeVector<CurrentLong-180);
LongLimit = [min(LongitudeVector),max(LongitudeVector)];
CurrentLat = LatitudeVector(CurrentSubs(1)); CurrentLong = LongitudeVector(CurrentSubs(2));
% Start to plot
fig = figure;
fig.Position = [100 100 1e3 0.95e3];
% big map
gx = geoaxes;
geobasemap('streets-light'); hold on 

geolimits(gx,LatLimit,LongLimit);
hold on   

ax = axes; axis([0 1 0 1]);ax.Visible = 'off'; 
hold on
[~,HeavyIdx] = maxk(BetaFullScaled,20);    % plot 20 largest entries                              
% select most heavy locations to plot (otherwise, plot super slow)
cmap = colormap(gx,flipud(hot));
for jj = HeavyIdx'
    [idx1,idx2] = ind2sub([NumLat,NumLong],jj);
    neig_lat = LatitudeVector(idx1); neig_long = LongitudeVector(idx2);
    colorcode = cmap(max(1,round((BetaFullScaled(jj))*256)),:);
    point = geoplot(gx,min(90,max(-90,neig_lat+[-1.25,1.25])),neig_long*ones(1,2),...
        'LineWidth',SmallPointSize,'Color',colorcode);
    point.Color(4) = BetaFullScaled(jj);
    hold on
end
geoplot(gx,CurrentLat+[-1.25,1.25],CurrentLong*ones(1,2),...
        'LineWidth',SmallPointSize,'Color',MainPointColor);
hold on
cb = colorbar(gx,'southoutside');
cb.Position = [0.25 0.2 0.55 0.02];
%     cdata = cb.Face.Texture.CData; %% trying to get transperent colorbar
%     but fail:))))))
%     cdata(end,:) = uint8(linspace(0,255,256));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;
hold off 

% zoom map
geoplot(gx,[-10,10,10,-10,-10]+CurrentLat,[-10,-10,10,10,-10]+CurrentLong,'Color',ZoomBoxColor);
ax2 = axes('position',[.55 .55 .35 .35]); box on % put box around new pair of axes
ax2.Visible = 'off'; ax2.XColor = 'red';
for jj = HeavyIdx'
    [idx1,idx2] = ind2sub([NumLat,NumLong],jj);
    neig_lat = LatitudeVector(idx1); neig_long = LongitudeVector(idx2);
    colorcode = cmap(max(1,round((BetaFullScaled(jj))*256)),:);
    point = geoplot(min(90,max(-90,neig_lat+[-1.25,1.25])),neig_long*ones(1,2),...
        'LineWidth',ZoomPointSize,'Color',colorcode);
    point.Color(4) = BetaFullScaled(jj);
    hold on
end
gx2 = point.Parent; 
geoplot(gx2,CurrentLat+[-1.25,1.25],CurrentLong*ones(1,2),...
        'LineWidth',ZoomPointSize,'Color',MainPointColor);
hold on
gx2.MapCenter = [CurrentLat,CurrentLong];
geolimits(gx2,[-10,10]+CurrentLat,[-10,10]+CurrentLong);
gx2.AxisColor = ZoomBoxColor;
gx2.LatitudeLabel.String = '';gx2.LongitudeLabel.String = '';
hold off   
end


function show_sonoma(kmlFile, opts)
    if nargin < 2, opts = struct; end
    if ~isfield(opts,'N'),      opts.N = 2000;    end
    if ~isfield(opts,'method'), opts.method = 'makima'; end

    if nargin < 1 || isempty(kmlFile)   %kml file declaration
        [f,p] = uigetfile({'*.kml','KML files (*.kml)'}, 'Choose a KML polyline');
        if isequal(f,0), disp('Canceled.'); return; end
        kmlFile = fullfile(p,f);
    end
    assert(exist(kmlFile,'file')==2, 'File not found: %s', kmlFile);

  
    [lat,lon,alt] = kml_read_longest_linestring(kmlFile);
    assert(numel(lat)>=2,'Could not find a valid LineString in the KML.');

    % ENU anchor 
    if ~isfield(opts,'lat0'), opts.lat0 = lat(1); end
    if ~isfield(opts,'lon0'), opts.lon0 = lon(1); end
    if ~isfield(opts,'h0'),   opts.h0   = (isempty(alt)||all(~isfinite(alt))) * 0 + ...
                                        (~isempty(alt))*median(alt,'omitnan'); end

    [x_raw,y_raw,z_raw] = geodetic2enu_local(lat,lon,alt,opts.lat0,opts.lon0,opts.h0);

    [x,y,z] = resample_polyline_strict(x_raw,y_raw,z_raw,opts.N,opts.method); %ENU 

    L = sum(hypot(diff(x),diff(y)));
    ttl = sprintf('%s — %d pts, length %.1f m', strip_filename(kmlFile), numel(x), L);

    % ----- plot -----
    fig = figure('Name','KML Track Viewer','Color','w');
    t = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

    nexttile(t,1);  % top-down ENU
    plot(x,y,'-','LineWidth',2); axis equal; grid on
    xlabel('E (m)'); ylabel('N (m)');
    title({'Track (top-down ENU)', ttl});

    nexttile(t,2);  % 3D
    plot3(x,y,z,'-','LineWidth',2); grid on; axis equal
    xlabel('E (m)'); ylabel('N (m)'); zlabel('Up (m)');
    title('Track (3D)'); view(3); rotate3d on
end

% ==================== helpers ====================

function [lat,lon,alt] = kml_read_longest_linestring(fname)
% Parse KML and return the longest <coordinates> block as a polyline.
    txt = fileread(fname);
    blocks = regexp(txt,'<coordinates[^>]*>(.*?)</coordinates>','tokens');
    best = struct('lat',[],'lon',[],'alt',[],'score',-inf);

    for b = 1:numel(blocks)
        s = strtrim(blocks{b}{1});
        if isempty(s), continue; end
        rows = regexp(s,'\s+','split');  % rows like "lon,lat[,alt]"
        llon = zeros(0,1); llat = zeros(0,1); lalt = zeros(0,1);
        for k = 1:numel(rows)
            if isempty(rows{k}), continue; end
            p = regexp(rows{k},',','split');
            if numel(p) >= 2
                llon(end+1,1) = str2double(p{1}); %#ok<AGROW>
                llat(end+1,1) = str2double(p{2});
                if numel(p) >= 3, lalt(end+1,1) = str2double(p{3}); else, lalt(end+1,1) = 0; end
            end
        end
        ok = isfinite(llon)&isfinite(llat)&isfinite(lalt);
        llon=llon(ok); llat=llat(ok); lalt=lalt(ok);
        if numel(llon) < 2, continue; end

        % if polygon is closed, drop the duplicate last point
        if hypot(llon(end)-llon(1), llat(end)-llat(1)) < 1e-12
            llon(end)=[]; llat(end)=[]; lalt(end)=[];
        end

        score = sum(hypot(diff(llon),diff(llat))); % proxy for path length
        if score > best.score
            best = struct('lat',llat,'lon',llon,'alt',lalt,'score',score);
        end
    end
    lat = best.lat; lon = best.lon; alt = best.alt;
end

function [x,y,z] = geodetic2enu_local(lat,lon,alt,lat0,lon0,h0)
% Use Mapping Toolbox if present, otherwise a small-area approximation.
    try
        wgs84 = wgs84Ellipsoid('meters');
        [x,y,z] = geodetic2enu(lat,lon,alt,lat0,lon0,h0,wgs84);
    catch
        Re = 6378137;                   % Earth radius (m)
        x = (lon - lon0) * pi/180 .* (Re*cosd(lat0));
        y = (lat - lat0) * pi/180 .* Re;
        if isempty(alt), alt = zeros(size(lat)); end
        z = alt - h0;
    end
end

function [xo,yo,zo] = resample_polyline_strict(x,y,z,N,method)
    if nargin < 5 || isempty(method), method = 'makima'; end

    x = x(:); y = y(:);
    if nargin < 3 || isempty(z), z = zeros(size(x)); else, z = z(:); end

    ok = isfinite(x)&isfinite(y)&isfinite(z);
    x=x(ok); y=y(ok); z=z(ok);

    if numel(x) >= 2
        d = hypot(diff(x),diff(y));
        keep = [true; d > 1e-12];   % drop consecutive duplicates
        x=x(keep); y=y(keep); z=z(keep);
    end
    assert(numel(x)>=2,'Polyline has <2 valid points after cleaning.');

    s = [0; cumsum(hypot(diff(x),diff(y)))];   % arc length
    [s, ia] = unique(s,'stable');              % keep vectors in sync
    x = x(ia); y = y(ia); z = z(ia);
    assert(numel(s)==numel(x),'Internal length mismatch after unique().');

    so = linspace(0,s(end),N).';
    xo = interp1(s,x,so,method);
    yo = interp1(s,y,so,method);
    zo = interp1(s,z,so,method);
end

function s = strip_filename(p)
    [~,s,ext] = fileparts(p); s = [s ext];
end

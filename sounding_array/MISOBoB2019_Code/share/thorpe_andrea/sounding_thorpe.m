% Try calculating thorpe scale for the soundings. Using data from the
% ascent which has been processed because the drops don't contain height
clear
%% User-defined fields

savefig = 1;

plev = 50;
mixing_eff = 0.25;
z2 = 5:5:(26*1000); % for regrid. 5-m spacing up to 26km

datdir = '/Volumes/cruise/CURRENT/share/data/radiosonde/Ride/';
figdir = '/Users/andrea/Documents/miso-bob/soundings/plots/turb/';
%% Load data and calculate reordered theta_v profile for each

% List files
files  = dir([datdir 'raw_data_second*.txt']);
files  = {files.name};
nfiles = length(files);

for ifile = 1:nfiles
    
    disp(['Working on ' num2str(ifile) '/' num2str(nfiles)])
    
    % Read using readtable
    dat = readtable([datdir files{ifile}]);

    % Calculate virtual potential temperature, and read in height
    % theta_v = theta * (1 + 0.61r); r = water vapor mixing ratio
    theta = cellfun(@str2num,dat.PotTemp(2:end));
    q     = cellfun(@str2num,dat.SpHum(2:end))/1000;
    z     = cellfun(@str2num,dat.HeightMSL(2:end));

    % Calculate thorpe scales
    [theta_v(ifile,:), epsi(ifile,:), K(ifile,:)] = ...
        eddyDiffusivityThorpe(theta, q, z, z2, plev, mixing_eff);
    
    % Get time from file name
    ymd_zz = files{ifile}(17:29);
    time = datenum(ymd_zz, 'yyyymmdd_HHMM');
    doy(ifile) = time - datenum(year(datetime(time,'convertfrom','datenum')),...
        1,1) + 1;
    
    % Grab rh and regrid to z2. First have to delete rows
    % that have duplicate values at a single level.
    [~, ui] = unique(z);
    row2del = setdiff(1:length(z), ui);
    rh_tmp = cellfun(@str2num,dat.RH(2:end));
    z(row2del)      = [];
    rh_tmp(row2del) = [];
    rh(ifile,:) = interp1(z, rh_tmp, z2);
    
    % From rh and temperature, calculate rh over ice for temperatures that
    % are below freezing
    T = cellfun(@str2num,dat.Temp(2:end));
    T(row2del) = [];
    T = interp1(z, T, z2); % interpolate to fixed grid
    vpr_pres = ClausiusClapeyron(T) .* (rh(ifile,:)/100);
    rh_ice(ifile,:) = vpr_pres ./ ps_ice(T-273.15);
    rh_ice(ifile,T > 273.15) = 0; % set rh to 0 above freezing

end

disp('Done')
%% Plot
% Filled contour of epsilon (TKE dissipation rate) (or eddy diffusivity)
% with time, with line contours of virtual potential temperature overlaid

% Contour levels
colmap = parula(12);
clev_thetav = 295:5:600;

figure('color','white','Position',[254 106 1034 382],'colormap',colmap); hold on
pp = pcolor(doy, z2/1000, log10(epsi)');
pp.EdgeColor = 'none';
contour(doy, z2/1000, theta_v', clev_thetav, 'linecolor','k',...
    'showtext','on','labelspacing',550)

set(gca,'FontSize',16)
xlabel('Year Day')
ylabel('Height (km)')
ylim([0 18])
caxis([-4 -1.0])
ch = colorbar;
ch.Label.String = 'log_{10} \epsilon';

if savefig
    figname = [figdir 'thorpe_epsilon.jpg'];
    export_fig(figname, '-r200')
end

%% Plot with clouds
% Filled contour of epsilon (TKE dissipation rate) (or eddy diffusivity)
% with filled contour showing relative humidities above a certain threshold
% Note: Made a scatter plot comparing the rh over ice to the value for
% epsilon. There doesn't appear to be any relationship between the two
% based only on the scatterplot

% Contour levels
colmap = parula(12);
clev_thetav = 295:5:600;
colmap_cloud = flipud(greys(100));
colmap_cloud(75:end,:) = [];

figure('color','white','Position',[254 106 1034 382]); hold on
% Contour clouds
contour(doy, z2/1000, rh',90:1:95,...
    'linecolor',[1 1 1] * 0.8,'linewidth',4);
contour(doy, z2/1000, rh',96:1:100,...
    'linecolor',[1 1 1] * 0.6,'linewidth',4);
% Contour ICE clouds
contour(doy, z2/1000, 100*rh_ice',125:5:149,...
    'linecolor',[1 1 1] * 0.8,'linewidth',4);
contour(doy, z2/1000, 100*rh_ice',150:5:250,...
    'linecolor',[1 1 1] * 0.6,'linewidth',4);

% Contour epsilon
pp = pcolor(doy, z2/1000, log10(epsi)');
pp.EdgeColor = 'none';
contour(doy, z2/1000, theta_v', clev_thetav, 'linecolor','k',...
    'showtext','on','labelspacing',550)

set(gca,'FontSize',16)
xlabel('Year Day')
ylabel('Height (km)')
ylim([0 18])
caxis([-4 -1.0])
ch = colorbar;
ch.Label.String = 'log_{10} \epsilon';

if savefig
    figname = [figdir 'thorpe_epsilon_wclouds.jpg'];
    export_fig(figname, '-r200')
end

%% Plot some information: frequency and mean value

figure('color','white')
gap = [0.03 0.03];
mh  = [0.17 0.025];
mw  = [0.1 0.04];
[ax axpos] = tight_subplot(1,2,gap,mh,mw);

% Frequency turbulent
subplot(1,2,1)
    freqturb = sum(isfinite(epsi),1)/size(epsi,1);
    plot(freqturb, z2/1000, 'linewidth', 2)
    ylim([0 18])
    set(gca,'FontSize',16, 'position',axpos{1},'xtick',0:0.2:0.8)
    xlabel({'Frequency of Significant','Turbulence Detected'})
    ylabel('Height (km)')
    grid on

% Mean value of epsilon
subplot(1,2,2)
    plot(nanmean(epsi,1), z2/1000, 'linewidth', 2)
    ylim([0 18])
    set(gca,'FontSize',16, 'position',axpos{2},'xscale','log',...
        'xtick',[1e-5 1e-4 1e-3 1e-2],'yticklabel',[])
    xlabel('Mean \epsilon')
    grid on
    
if savefig
    figname = [figdir 'thorpe_epsilon_stats.jpg'];
    export_fig(figname, '-r200')
end
%% Scatter plot for all radiosondes
% figure('position',[500 500 790 300], 'color','white'); hold on
% colmap = (jet(nfiles));
% ptsize = 80;
% 
% for ifile = 1:nfiles
%     subplot(1,2,1); hold on
%         scatter(epsi(ifile,:), z2/1000, ptsize, colmap(ifile,:), '.')
%         
%     subplot(1,2,2); hold on
%         scatter(K(ifile,:), z2/1000, ptsize, colmap(ifile,:), '.')
% end
% 
% subplot(1,2,1)
%     set(gca,'fontsize',16,'xscale','log')
%     ylabel('Height (km)')
%     xlabel('Epsilon')
%     xlim([0 0.2])
%     grid on
% subplot(1,2,2)
%     xtx = [1e-1 1e0 1e1];
%     set(gca,'fontsize',16,'xscale','log','xtick',xtx)
%     ylabel('Height (km)')
%     xlabel('Eddy Diffusivity')
%     grid on
%% FUNCTION
% Compute the eddy diffusivity coefficients for turbulent sections of
% soundings

function [theta_v, epsi, K] = eddyDiffusivityThorpe(T, q, z, z2, p, mixing_eff)
    % Input variables
    %   T  - potential temperature (K)
    %   q  - specific humidity (kg/kg)
    %   z  - height coordinate
    %   z2 - regrid to this vertical grid
    %   p  - confidence level for separating inversions from overturns
    %   mixing_eff - mixing efficiency to calculate eddy diffusivity
    %                from epsilon
    % Output variables
    %   theta_v - interpolated and sorted virtual potential temperature
    %   epsi    - TKE dissipation rate
    %   K       - eddy diffusivity
    
    g = 9.81; % gravity
    
    % Use the wilson2010_lookup function to load the table of tnr values
    % that need to be surpassed in order to detect a physical overturn from
    % a noise-induced inversion. The function takes a confidence level
    % (/100) as input. Note: this can be passed as input into the function
    % instead of calling within the function if this is being applied to
    % many profiles to increase speed (readtable function takes ~0.05s)
    tnrmin = wilson2010_lookup(p);
    
    % Compute virtual potential temperature from temperature and 
    % humidity information
    r        = q ./ (1-q); % mixing ratio needed for theta_v
    theta_v  = T .* (1 + (0.61 * r));
    
    % Regrid theta_v to 5-m vertical resolution. First have to delete rows
    % that have duplicate values at a single level.
    [~, ui] = unique(z);
    row2del = setdiff(1:length(z), ui);
    z(row2del)       = [];
    theta_v(row2del) = [];
    
    theta_v = interp1(z, theta_v, z2);
    
    % Estimate the noise variance by removing a linear trend from short
    % segments of data and calculating the mean of the squared differences.
    % The average from all of the sequences is an estimate of twice the
    % noise variance. Use a segment of length 5 obs
    nlevs = length(z2);
    for ilev = 1:nlevs-4
        seqvar(ilev) = nanmean(nandetrend(...
            z2(ilev:ilev+4),theta_v(ilev:ilev+4)).^2);
    end
    
    sigma_n = sqrt(0.5 * nanmean(seqvar));
    
    % Compute an intermediate profile of theta_v by sweeping the profile
    % from bottom to top and keeping the value for a succeeding point the
    % same as the preceeding point if the difference between the two does
    % not exceed measurement noise (sigma_n). Do this twice, the second
    % time sweeping from the top and then average the two profiles. 
    theta_v_int1 = theta_v;
    theta_v_int2 = theta_v;
    for ilev = 2:nlevs
        tdiff = abs(theta_v(ilev) - theta_v_int1(ilev-1));
        
        if tdiff < sigma_n
            theta_v_int1(ilev) = theta_v_int1(ilev-1);
        end
    end
    for ilev = nlevs-1:-1:1
        tdiff = abs(theta_v(ilev+1) - theta_v_int2(ilev));
        
        if tdiff < sigma_n
            theta_v_int2(ilev) = theta_v_int2(ilev+1);
        end
    end
    theta_v_int = (theta_v_int1 + theta_v_int2)/2;
    
    % Compute displacements from a sorted profile using the intermediate
    % profile. Displacements - D
    [theta_v_sort inds D] = thorpeSort(theta_v_int, z2);
    
    % Compute the brunt-vaisala frequency using the sorted profile
    % N = sqrt(g/theta * dtheta/dz)
    dtheta_dz = forwDiff_v(theta_v_sort, z2(2) - z2(1));
    N2 = (g.* dtheta_dz)./ theta_v_sort; % ignore unstable areas
    
    % Locate overturns using the displacements profile as places where: the
    % sum of the displacements over the length is zero; the displacements
    % are negative in the lower half and positive in the upper half
    % --> Inversions are places where there are consecutive non-zero values
    cD = cumsum(D);
    A = cD ~= 0;
    n = double(diff([~A(1);A(:)]) == 1);
    v = accumarray(cumsum(n).*A(:)+1,1);
    n(n == 1) = v(2:end);
    ii_start = find(n ~= 0); % Index of the first point of the inversion
    ii_end = ii_start + n(ii_start)-1; % Index of the last point
    % --> Inversions also have the requirement that the displacements are
    % generally negative at the bottom & positive at the top
    if any(cD > 0)
        error(['One of the inversions violates condition 2.. ' ...
            'recode so this inversion is deleted'])
    end
    
    % Use the sorted profile to determine the local trend-to-noise
    % ratio (TNR), which is the local trend divided by the standard
    % deviation of the noise
    %          (X(i+1) - X(i-1)) / (2 * sigma_n)
    tnr = NaN([nlevs 1]);
    for ilev = 2:nlevs-1
        tnr(ilev) = (theta_v_sort(ilev+1) - theta_v_sort(ilev-1)) ...
            / (2 * sigma_n);
    end
    
    % Separate physical overturns from inversions due to noise. Choose a
    % confidence level, and then from a look-up-table of the range for n
    % independent identically distributed noise that are normally
    % distributed between 0 and 1, determine the minimum tnr that needs to
    % be surpassed for the overturn. Use the average tnr over the inversion
    % to determine the tnr. 
    % Loop through each inversion
    
    epsi = NaN([nlevs 1]);
    K = NaN([nlevs 1]);
   
    for ii = 1:length(ii_start)
        % Determine whether or not the inversion can be classified as an
        % overturn
        tnr_inv = nanmean(tnr(ii_start(ii):ii_end(ii)));
        n_inv   = ii_end(ii) - ii_start(ii) + 1;
        
        % If we can classify the inversion as an overturn
        if tnr_inv > tnrmin(n_inv)
            % Calculate the thorpe length
            thorpe_L = rms(D(ii_start(ii):ii_end(ii)));

            % Calculate epsilon
            % Ck * L_T^2 * N^3
            epsi(ii_start(ii):ii_end(ii)) = ...
                1.0 * thorpe_L^2 * ...
                nanmean(N2(ii_start(ii):ii_end(ii)))^(3/2);
            
            % Calculate the eddy diffusivity
            K(ii_start(ii):ii_end(ii)) = mixing_eff .* ...
                epsi(ii_start(ii):ii_end(ii)) ./ ...
                nanmean(N2(ii_start(ii):ii_end(ii)));
        end 
    end
end
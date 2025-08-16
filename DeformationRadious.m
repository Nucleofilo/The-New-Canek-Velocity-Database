clear; clc; close all;

% Load your CTD file (replace with correct path if needed)
load('D:\Papers\Paper_III\Data_Methods\Raw_per_Campaign\Canek\Canek_CTD_historico');
lom = lon;
lam = lat;

[yucx, yucy] = muadro( [-86.6128 -84.9265], [21.4965 21.8998]);
[flox, floy] = muadro( [-82.0871 -80.7244], [22.9682 24.6864]);

inyuc = inpolygon(lom,lam, yucx, yucy);
inflo = inpolygon(lom,lam, flox, floy);
into = inyuc | inflo;


lom = lom(into);
lam = lam(into);
sp = sp(into);
p = p(into);
te = te(into);
fe = fe(into);


[sa, in_ocean] = gsw_SA_from_SP(sp(:),p(:),lom(:),lam(:));
ct = gsw_CT_from_t(sa,te(:),p(:));
sigt = gsw_sigma0(sa,ct);

plot(lom, lam, 'o')

%%
% Initialize output variables
Ro = []; lat = []; lon = []; depth = [];

% Unique cast identifiers
lances = unique(fe); 
ncasts = length(lances); 
aa = 1;

for i = 1 : ncasts
    index = (fe == lances(i)); % index of current cast

    % Extract variables
    sdu = sa(index); 
    tdu = ct(index);
    pdu = p(index);
    ladu = lam(index);
    londu = lom(index);

    % Remove NaNs
    bad = isnan(sdu + tdu + pdu + ladu);
    sdu(bad) = []; tdu(bad) = []; pdu(bad) = []; ladu(bad) = [];

    % Only process casts with enough vertical range
    if min(pdu) <= 20 && max(pdu) >= 650 
        % Compute buoyancy frequency
        [N2, p_mid] = gsw_Nsquared(sdu, tdu, pdu, mean(ladu));
        N2(N2 < 0) = 0;
        N = sqrt(N2); % units of s^-1

        % Convert pressure to depth (positive downward)
        z_mid = -gsw_z_from_p(p_mid, mean(ladu)); % meters

        % Interpolation grid in depth (positive downward)
        dz = 10; 
        z_target = (0:dz:max(pdu))';

        % Interpolate N(z) onto uniform z grid
        try
        N_interp = interp1(z_mid, N, z_target, 'linear', 'extrap');

        % Compute local Coriolis parameter
        f = coriolis(mean(ladu)); % in s^-1

        % Apply WKB formula (Chelton Eq. A.20), m = 1
        Ro_m = (1 / (pi * f)) * trapz(z_target, N_interp); % in meters

        % Save results
        Ro(aa)   = Ro_m * 1e-3; % convert to km
        lat(aa)  = mean(ladu);
        lon(aa)  = mean(londu);
        depth(aa) = max(pdu);

        inyu(aa) = inpolygon(lon(aa),lat(aa), yucx, yucy);
        catch

        end
        aa = aa + 1;
    end
end
Ro(Ro < 10) = NaN; 
% Optional: scatterplot
figure;
scatter(depth(inyu), Ro(inyu), 30, 'filled'); hold on;
scatter(depth(~inyu), Ro(~inyu), 30, 'filled');
grid on;
ll = legend('Yucatan Channel', 'Florida Strait');
xlabel('Depth (m)');
ylabel('\lambda_1 (km)');
title('First Baroclinic Rossby Radius (WKB Approx.)');
%%

close all

figurax;set(gcf, 'pos', [-1918 258 1414 738])
axes('pos', [0.5552 0.2290 0.3498 0.6960])
[xv, yv] = muadro([-86.6296 -84.8879], [21.4141 22.0294]);
inyu = inpolygon(lon, lat, xv, yv);

plot(depth(inyu), Ro(inyu), 'ok', 'markerfacecolor',  rgb('RoyalBlue'))

b = polyfit(depth(inyu), Ro(inyu), 3);
yp = polyval(b, 600:2100);
hold on;
plot(600:2100, yp, 'b', 'linewidth', 3)

rmyu = mean(yp)
[min( Ro(inyu)) max( Ro(inyu))]


%

hold on;
[xv, yv] = muadro([-82.1774 -80.4355], [22.9898 24.5000]);
infl = inpolygon(lon, lat, xv, yv);

plot(depth(infl),Ro(infl), 'sk', 'markerfacecolor', rgb('OrangeRed'), 'markersize', 9)

b = polyfit(depth(infl), Ro(infl), 2);
yp = polyval(b, 750:1600);
hold on;
plot(750:1600, yp, 'r', 'linewidth', 3)

rmflo = mean(yp)
[min( Ro(infl)) max( Ro(infl))]

title( '\textbf{First Baroclinic Rossby Radius of Deformation}' );

xlabel('Profile depth (m)')
ylabel('Rossby radius of deformation (km)')


ax0(1) = plot(NaN,NaN, 'sk', 'markerfacecolor', rgb('OrangeRed'), 'markersize', 9);
ax0(2) = plot(NaN,NaN, 'ok', 'markerfacecolor',  rgb('RoyalBlue'));


ll = legend(ax0, {'Florida', 'Yucatan'}, 'Location', 'Best');
ll.FontSize = 13;
ll.Position = [0.5584 0.8546 0.0803 0.0625];
% text(627.2189, 97.7688, '(Chelton et al., 1998)', 'FontSize', 13);
grid on;

Fs = findall(gcf,'-property','FontSize');
for k = 1 : length(Fs)
    if strcmp(get(Fs(k), 'type'), 'axes')
        if Fs(k).FontSize < 12
            Fs(k).FontSize = 12;
        end
    end
end

set(findall(gcf,'-property','Interpreter'),'Interpreter', 'latex');
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex');

text(2078, 98.5564, '\textbf{(I)}','Interpreter', 'latex', 'fontsize', 17);

%%
axes('pos', [0.753889674681754 0.234540527223454 0.147807705924052 0.19308951974927]);
box on;

plot(lon(infl), lat(infl), '.', 'color', 'r', 'markersize', 7);
hold on;

plot(lon(inyu), lat(inyu), '.', 'color',  'b', 'markersize', 7);
hold on;

axis equal tight
axis([-87.3, -80.7, 20.5, 25])
text(-87.1071, 24.2678, '\textbf{(II)}','Interpreter', 'latex', 'fontsize', 17);
set(gca, 'YtickLabels', {}, 'XtickLabels', {})

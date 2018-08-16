% ShieldingWrapperScript.m is a wrapper script for calculating topographic and slope
% shielding of cosmogenic nuclide production for the simplified catchment examples
% presented in DiBiase (2018)
%  
% The wrapper script reproduces the model as described in the main text, and generates
% the model results for Figures 3, 4, and 5 of DiBiase (2018). The figures generated for
% the paper were run at a resolution of 500, which takes many hours to run on a typical desktop.
% 
%
% Required functions:
%    SkylineFunction.m     Function for calculating skyline shielding integration limit
%                            as a function of azimuth
%
%    IntensityFunction.m   Function for calculating the normalized radiation flux as a
%                            function of slope position and depth below the surface
% 
% Reference:
%    DiBiase, R.A., 2018. Short Communication: Increasing vertical attenuation length
%        of cosmogenic nuclide production on steep slopes negates topographic shielding
%        corrections for catchment erosion rates. Earth Surf. Dynam. Discuss., in review,
%        https://doi.org/10.5194/esurf-2018-48
%
%
% Last updated by Roman DiBiase rad22@psu.edu 8/15/2018
%
%

tic

% define model parameters
resolution = 100;                           % resolution of model (500 is value used for figures in paper, but slow to run)

alpha_deg = 0:10:80;                        % vector of hillslope angles, degrees
alpha = deg2rad(alpha_deg);                 % conversion of hillslope angle to radians

Lh = 500;                                   % horizontal hillslope length (m)
dx = Lh/resolution;                         % horizontal grid spacing (m)
x = (0:dx:Lh);                              % vector of mapview hillslope distances (m)

lambda = 160;                               % nominal mass attenuation length (g cm^-2)
lambda_uni = 1.3*lambda;                    % mass attenuation length for unidirectional incoming radiation (g cm^-2)
m = 2.3;                                    % exponent of intensity function (typically assumed to be 2.3)
density = 2.7;                              % rock density (g cm^-3)
dz = lambda/(resolution*density)/100;       % vertical grid spacing (m)
z = (0:dz:(40*lambda/density/100))';        % vector of vertical depths (m)

UnshieldedProductionProfile = zeros(size(z));           % production on horizontal unshielded surface

% initialize variables for "interior" catchment calculation
S_interior = zeros(length(z),length(x));                % temporary array of shielding factors for interior catchment case
MeanInteriorShieldingFactor = zeros(length(alpha),1);   % vector of mean shielding factors as a function of hillslope angle
depth5pc = zeros(length(alpha),length(x));              % vector of depth at which production is 5% of local surface production
C_eff_interior = zeros(length(alpha),length(x));        % array containing effective shielding as a function of distance for each interior case
S_interior_grid = cell(length(alpha_deg),1);            % cell array containing shielding factor data for all interior cases

% initialize variables for "exterior" catchment calculation
S_exterior = zeros(length(z),length(x));                % temporary array of shielding factors for exterior catchment case
MeanExteriorShieldingFactor = zeros(length(alpha),1);   % vector of mean shielding factors as a function of catchment dip angle
C_eff_exterior = zeros(length(alpha),length(x));        % array containing effective shielding as a function of distance for each exterior case
S_exterior_grid = cell(length(alpha_deg),1);            % cell array containing shielding factor data for all exterior cases

% Figure sizing and preparation
Fig3PosX = [0.065 0.375 0.685 0.065 0.375 0.685 0.065 0.375 0.685];
Fig3PosY = [0.595 0.595 0.595 0.44 0.44 0.44 0.285 0.285 0.285];
Fig3 = figure(3);
set(Fig3,'Units','inches','Position',[1 2 8.5 11])
Fig4 = figure(4);
set(Fig4,'Units','inches','Position',[2 2 8.5 11])
Fig5 = figure(5);
set(Fig5,'Units','inches','Position',[3 2 8.5 11])


% Calculate production profile as a function of depth for unshielded case
for j = 1:length(z)
    % Define function handle for intensity function
    f_intensity = @(phi,theta)IntensityFunction(z(j),density,lambda_uni,0,phi,theta,m);
    
    % Solve Equation 5 of DiBiase (2018)
    UnshieldedProductionProfile(j) = 2.*integral2(f_intensity,0,pi,0,pi/2);
end

% Calculate proxy of surface concentration for unshielded sample (denominator of Equation 13 in DiBiase (2018)
C_unshielded = sum(UnshieldedProductionProfile);

% Main loop through hillslope angles/catchment dip, across all hillslope positions, and at all depths
for n = 1:length(alpha)
    toc
    for i = 2:length(x)
        disp(['alpha = ' num2str(alpha_deg(n)) ' deg; x = ' num2str(x(i)) ' m/' num2str(Lh) ' m']) % display progress        
        for j = 1:length(z)
            
            % Define function handles for intensity and skyline functions for interior case
            f_intensity = @(phi,theta)IntensityFunction(z(j),density,lambda_uni,alpha(n),phi,theta,m);
            f_sky_interior = @(phi)SkylineFunction( x(i), z(j), Lh, alpha(n), phi, 1);
            
            % Solve Equation 5 of DiBiase (2018)
            S_interior(j,i) = 2.*integral2(f_intensity,0,pi,f_sky_interior,pi/2);
            
            % Define function handle for skyline function for exterior case
            f_sky_exterior = @(phi)SkylineFunction( x(i), z(j), Lh, alpha(n), phi, 0);
            
            % Solve Equation 5 of DiBiase (2018)
            S_exterior(j,i) = 2.*integral2(f_intensity,0,pi,f_sky_exterior,pi/2);
        end
        
        % Determine depth at which production is 5% of surface as proxy for effective attenuation length (Equation 9 of DiBiase (2018)
        depth5pc(n,i) = z((abs(S_interior(:,i)./S_interior(1,i)-0.05)==min(abs(S_interior(:,i)./S_interior(1,i)-0.05))));
    end
    
    % Assume no shielding of ridgeline sample
    S_interior(:,1)=UnshieldedProductionProfile;
    
    % Save key data for plotting (interior case)
    S_interior_grid{n} = S_interior;
    C_eff_interior(n,:) = sum(S_interior_grid{n})./C_unshielded;        % Equation 13 in DiBiase (2018)
    depth5pc(n,1) = z((abs(S_interior(:,1)./S_interior(1,1)-0.05)==min(abs(S_interior(:,1)./S_interior(1,1)-0.05))));
    MeanInteriorShieldingFactor(n) = mean(C_eff_interior(n,2:end));     % Equation 14 in DiBiase (2018)
    
    % Save key data for plotting (exterior case)
    S_exterior(:,1)=UnshieldedProductionProfile;
    S_exterior_grid{n} = S_exterior;    
    C_eff_exterior(n,:) = sum(S_exterior_grid{n})./C_unshielded;        % Equation 13 in DiBiase (2018)
    MeanExteriorShieldingFactor(n) = mean(C_eff_exterior(n,2:end));     % Equation 14 in DiBiase (2018)
    
    % Assign line colors for Figure 4        
    switch alpha_deg(n)
        case 0
            LineColor = [0 0 0];
        case 10
            LineColor = [0 0 0.72];
        case 20
            LineColor = [0 0.1 1];
        case 30
            LineColor = [0 0.45 1];
        case 40
            LineColor = [0 0.82 0.6];
        case 50
            LineColor = [0 0.66 0.07];
        case 60
            LineColor = [1 0.9 0];
        case 70
            LineColor = [1 0.62 0];
        case 80
            LineColor = [1 0 0];
    end
    
    
    % Generate Figure 3 of DiBiase (2018)
    figure(3)
    hold on
    subplot('Position', [Fig3PosX(n) Fig3PosY(n) 0.28 0.115])
    imagesc([0 1],[0 max(z).*density./lambda.*100],log(S_interior_grid{n}),[-3 0])
    xlim([0 1])
    ylim([0 5])    
    xticks(0:.2:1)
    yticks(0:1:5)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','FontSize',7)
    title(['\alpha = ' num2str(alpha_deg(n)) '\circ'],'FontName','Times','FontSize',9)
    
    % Generate Figure 4 of DiBiase (2018)
	figure(4)
    subplot('Position', [0.32 0.61 0.35 0.19])
    hold on
	plot(x(2:end)/Lh,S_interior_grid{n}(1,2:end),'Color',LineColor)
    xlim([0 1])
    ylim([0 1])
    box on    
    xticks(0:.2:1)
    yticks(0:.2:1)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','FontSize',7)
    ylabel('Surface skyline shielding, S_0','FontName','Times','FontSize',8)
    
    subplot('Position', [0.32 0.38 0.35 0.19])
    hold on
    plot(x(2:end)./Lh,depth5pc(n,2:end)./(3.*lambda./density./100),'Color',LineColor)
    xlim([0 1])
    ylim([0 4])
    box on    
    xticks(0:.2:1)
    yticks(0:1:4)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','FontSize',7)
    ylabel('Normalized effective attenuation length, \Lambda/\Lambda_e_f_f','FontName','Times','FontSize',8)
	
	subplot('Position', [0.32 0.15 0.35 0.19])
	hold on
    plot(x(2:end)/Lh,C_eff_interior(n,2:end),'Color',LineColor)
    xlim([0 1])
    ylim([0 2])
    box on
    xticks(0:.2:1)
    yticks(0:.5:2)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','FontSize',7)
    xlabel('Normalized distance from ridgeline, x/L_{h}','FontName','Times','FontSize',8)
    ylabel('Total effective shielding factor, C_e_f_f','FontName','Times','FontSize',8)
end


% Clean up and print Figure 3 to pdf
figure(3)
subplot('Position', [Fig3PosX(1) Fig3PosY(1) 0.28 0.115])
ylabel('Normalized vertical depth, \rhoz/\Lambda','FontName','Times','FontSize',8)
subplot('Position', [Fig3PosX(4) Fig3PosY(4) 0.28 0.115])
ylabel('Normalized vertical depth, \rhoz/\Lambda','FontName','Times','FontSize',8)
subplot('Position', [Fig3PosX(7) Fig3PosY(7) 0.28 0.115])
ylabel('Normalized vertical depth, \rhoz/\Lambda','FontName','Times','FontSize',8)
xlabel('Normalized distance from ridgeline, x/L_{h}','FontName','Times','FontSize',8)
subplot('Position', [Fig3PosX(8) Fig3PosY(8) 0.28 0.115])
xlabel('Normalized distance from ridgeline, x/L_{h}','FontName','Times','FontSize',8)
subplot('Position', [Fig3PosX(9) Fig3PosY(9) 0.28 0.115])
xlabel('Normalized distance from ridgeline, x/L_{h}','FontName','Times','FontSize',8)
print('Figure3','-dpdf','-fillpage')

% Print Figure 4 to pdf
figure(4)
print('Figure4','-dpdf','-fillpage')

% Generate Figure 5 and print to pdf
figure(5)
subplot('Position', [0.12 0.7 0.33 0.12])
plot(alpha_deg,MeanInteriorShieldingFactor,'k.','MarkerSize',10)
box on
grid on
xlim([0 90])
ylim([0.8 1.2])
xticks(0:10:90)
yticks(0.8:.1:1.2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',7)
xlabel('Mean hillslope angle,  \alpha (degrees)','FontName','Times','FontSize',8)
ylabel('Mean shielding factor','FontName','Times','FontSize',8)
title('Interior catchments','FontName','Times','FontSize',9)

subplot('Position', [0.57 0.7 0.33 0.12])
plot(alpha_deg,MeanExteriorShieldingFactor,'k.','MarkerSize',10)
box on
grid on
xlim([0 90])
ylim([0.5 2.0])
xticks(0:10:90)
yticks(0.5:.5:2.0)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',7)
xlabel('Dip of plane fit to ridgelines,  \beta (degrees)','FontName','Times','FontSize',8)
ylabel('Mean shielding factor','FontName','Times','FontSize',8)
title('Exterior catchments','FontName','Times','FontSize',9)
toc
print('Figure5','-dpdf','-fillpage')

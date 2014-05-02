%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Communication Simulator        %%%%%%
%%%%%%%   Created by: Zachary Tschirhart %%%%%%
%%%%%%%   For Team PROPHACY Mission      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function communication_sim
sc_altitude  = 370; %km
earth_radius = 6378.1; %km
sc_radius    = sc_altitude + earth_radius; %km
planck_c     = 6.626e-34;
c            = 299792458; %m/s  
%[altitude, temp, pres, rho, c, g, mu, nu, k, n, n_sum] = ...
%    atmosphere(); 

%  Wavelength chart
%  http://www.lexellaser.com/techinfo_wavelengths.htm
%    Far-infrared   1.54e-6 meters
%    Near-infrared  7e-7 to 1.35e-6 meters
%    Visible        4.16e-7 to 6.94e-7 meters
%    Fiber optic comm operates around the 7.5e-7 to 1.4e-6 meters region
%    Atmosphere attenuation significantly increases around 1.45e-6
%     According to Laser experts, the best wavelength to use for
%     space communication, would be 1060, since Silicon's
%     transmisitivity peaks at this wavelength and any water
%     absorption can be neglected

%lambda = 4.16e-7:1e-9:1.54e-6;  % Full range, can change this
lambda = 1.06e-6;
Pt     = 0.1:0.01:14;            % .1 to 10 watts power transmitted
Dt     = (0.001:0.001:0.045);    % Diameter of the laser appature in meters
Dr     = Dt;  %Using the same apature on both TX and RX
range  = 300:100:(max_range(sc_radius,earth_radius).*1000); % distances 
                                                      % between sc
                                                      % meters
beam_waist = 0.0001; % beam waist in meters
TX_RX_system_loss = 1;%10^(-57.8/20);  % Unitless, converted from dB
                                    % More info: 
                                    %Free-Space Laser
                                    %Communications: 
                                    %Principles and Advances
                                    % By Arun K. Majumdar
nu = 100; % photons/bit is the typical sensitivity of a silicon receiver

i = 1;
for r = range
  [Pr(i),laser_area_at_RX] = OFSL(Pt(end), TX_RX_system_loss, ...
                             TX_RX_system_loss, beam_waist, ...
                             Dr(end), r, lambda);
  diameter = 2*sqrt(laser_area_at_RX/pi);
  max_datarate(i) = datarate(Pr(i), nu, lambda, c, planck_c);
  theta(i) = pointing_req(Dr(end), 2*sqrt(laser_area_at_RX/pi), r);
  i = i + 1;
end

  [Pr,laser_area_at_RX] = OFSL(Pt(end), TX_RX_system_loss, ...
                             TX_RX_system_loss, beam_waist, ...
                             Dr(end), 500000, lambda)
  diameter = 2*sqrt(laser_area_at_RX/pi)
  max_datarate = datarate(Pr, nu, lambda, c, planck_c)
  theta = pointing_req(Dr(end), 2*sqrt(laser_area_at_RX/pi), 500000)


  %figure('Name','Range vs. pointing requirements');
  %plot(range, theta);

  %figure('Name', 'Range vs. Power received');
  %plot(range, Pr);

  %figure('Name','Range vs. max datarate');
  %plot(range, max_datarate);

end


% Wrapper function for the atmo function, so we can see outputs
% easily
% Output:
%   altitude:   Total Reporting Altitudes[0<=alt<=1000 km][km]
%   temp:       Temperature array[0<=alt<=1000 km][K]
%   pres:       Pressure array[0<=alt<=1000 km][Pa]
%   rho:        Density array[0<=alt<=1000 km][kg/m^3]
%   c:          Speed of sound array[0<=alt<=86 km][m/s]
%   g:          Gravity array[0<=alt<=1000 km][m/s^2]
%   mu:         Dynamic Viscosity array[0<=alt<=86 km][N*s/m^2]
%   nu:         Kinematic Viscosity array[0<=alt<=86 km][m^2/s]
%   k:          Coefficient of Thermal Conductivity
%                 array[0<=alt<=86 km][W/(m*K)]
%   n:          Number Density of individual gases
%                 (N2 O O2 Ar He H)[86km<=alt<=1000km][1/m^3]
%   n_sum:      Number Density of total gases
%                 [86km<=alt<=1000km][1/m^3]
function [altitude, temp, pres, rho, c, g, mu, nu, k, n, n_sum] = ...
      atmosphere()
  [altitude, Z_L, Z_U, temp, pres, rho, c, g, mu, nu, k, n, n_sum] = atmo();
end


% Optical free space link equation Assuming flat beam power width
%   For more complex and accurate solution, add Gaussian diffusion.
% Input:
%   Pt     = Power transmitted (W)
%   Lt/r   = Transmitter/Receiver Loss 
%   Dt/r   = Transmitter/Receiver Aperture Diameter (m)
%   R      = Range (m)
%   lambda = Wavelength (m)
% Output:
%   Pr     = Power Received
function [Pr,laser_area_at_RX] = OFSL(Pt, Lt, Lr, Dt, Dr, R, lambda)
  apature_loss     = 1 - exp(-2);
  beam_div_angle   = lambda / (2*pi*Dt);
  power_density_TX = Pt/(pi*(Dt/2)^2);
  laser_area_at_RX = pi*(R*beam_div_angle)^2;
  power_density_RX = Pt/laser_area_at_RX;
  ideal_Pr         = power_density_RX*(pi*(Dr/2)^2);
  Pr               = Lr*Lt*ideal_Pr*apature_loss^2;
end


% Pointing accuracy 
% TODO: Gaussian equation and power losses
% Input:
%   Dr     = Receiver Aperture Diameter (m)
%   beam_d = Diameter of beam at Receiver (m)
%   R      = Range (m)
% Output:
%   theta  = +/- degrees pointing requirement

function [theta] = pointing_req(Dr, beam_d, R)
  theta = sind(abs(beam_d - Dr)/R);
end



%Maximum data rate that can be acheived with power received
% Input: 
%   Pr       = Power Received (W)
%   nu       = Receiver Sensitivity (photons/bit)
%   lambda   = Wavelength (m)
%   c        = Speed of light constant (m)
%   planck_c = Planck's constant
% Output:
%   dr = Max datarate (bits/sec)
function [dr] = datarate(Pr, nu, lambda, c, planck_c)
  v = c/lambda;
  dr = Pr/(planck_c*v*nu);
end

%Gaussian beam divergance and beam width
% Input: 
%   lambda = Wavelength (m)
%   w0     = Beam Waist (m)
%   R      = Range (m)
% Output:
%   w      = Spot size radius (m)
function [w] = beam_width(lambda, w0, R)
  zr = (pi*w0^2)/lambda;
  w = w0*sqrt(1 + (R/zr)^2);
end

% Function to find the absolute maximum communication range of
% satillites given they need to be within LOS
% Input:
%   sc_radius               = Radius of orbit from center of gravity (km)
%   earth_radius            = Radius of Earth (km)
% Output:
%   max_distance_between_sc = Distance between spacecraft (km)
function [max_distance_between_sc] = max_range(sc_radius, ...
                                               earth_radius)
  max_distance_between_sc = 2*sqrt(sc_radius^2 - earth_radius^2);
end

% Function to propogate altitude values over a beam traveled 
% Input:
%   distance_between_sc = Distance between spacecraft (km)
%   sc_radius           = Radius of orbit from center of gravity (km)
%   earth_radius        = Radius of Earth (km)
%   division            = Deltax of signal propogation 
% Output:
%   altitudes           = List of altitudes from center of earth (km)
%   valid               = Valid signal (> earth radius for now) 
function [altitudes, valid] = list_altitudes(distance_between_sc, ...
                                             sc_radius, ...
                                             earth_radius, ...
                                             division)
  angle_between_crafts = 2*asin(distance_between_sc/(2*sc_radius));
  lowest_radius = sc_radius * cos(angle_between_crafts / 2);

  if(lowest_radius <= earth_radius) valid = 0;
  else valid = 1;
  end
  
  i = 1;
  for x = -distance_between_sc/2:division:distance_between_sc/2
    altitude(i) = sqrt(x^2 + lowest_altitude^2) - earth_radius;
    i = i + 1;
  end
end


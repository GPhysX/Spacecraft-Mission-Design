%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Communication Simulator        %%%%%%
%%%%%%%   Created by: Zachary Tschirhart %%%%%%
%%%%%%%   For Team PROPHACY Mission      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function communication_sim
sc_altitude  = 370; %km
earth_radius = 6378.1; %km
sc_radius    = sc_altitude + earth_radius; %km
%[altitude, temp, pres, rho, c, g, mu, nu, k, n, n_sum] = ...
%    atmosphere(); 

%  Wavelength chart
%  http://www.lexellaser.com/techinfo_wavelengths.htm
%    Far-infrared   1.54e-6 meters
%    Near-infrared  7e-7 to 1.35e-6 meters
%    Visible        4.16e-7 to 6.94e-7 meters
%    Fiberoptics operate around the 7.5e-7 to 1.4e-6 meters region
%    Atmosphere attenuation significantly increases around 1.45e-6 meters
lambda = 4.16e-7:1e-9:1.54e-6;  % Full range, can change this
Pt     = 0.1:0.01:5;            % .1 to 5 watts power transmitted
At     = pi.*((0.001:0.001:0.045).^2); % Area of the laser
                                       % appature in meters
Ar     = At;  %Using the same apature on both TX and RX
range  = 1:1:(max_range(sc_radius,earth_radius).*1000); % distances 
                                                      % between sc meters
TX_RX_system_loss = 10^(-57.8/20);  % Unitless, converted from dB
                                    % More info: 
                                    %Free-Space Laser
                                    %Communications: 
                                    %Principles and Advances
                                    % By Arun K. Majumdar
normal_range = 1; %meters
Pr = OFSL(Pt(end), TX_RX_system_loss, TX_RX_system_loss, At(end), ...
          Ar(end), normal_range, lambda(1))


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


% Optical free space link equation
% Input:
%   Pt     = Power transmitted (W)
%   Lt/r   = Transmitter/Receiver Loss 
%   At/r   = Transmitter/Receiver Aperture Area (m)
%   R      = Range (m)
%   lambda = Wavelength (m)
% Output:
%   Pr     = Power Received

function [Pr] = OFSL(Pt, Lt, Lr, At, Ar, R, lambda)
%%Apature loss was added without verification, need to verify %%
  apature_loss = 1 - exp(-2)
  RX_gain      = (4*pi*Ar)/lambda^2;
  range_loss   = (lambda/(4*pi*R))^2;
  TX_gain      = (4*pi*At)/lambda^2;
  Pr           = Lr*RX_gain*range_loss*TX_gain*Lt*Pt* ...
                   apature_loss^2;
end

%Maximum data rate that can be acheived with power received
% Input: 
%   Pr = Power Received (W)
%   nu = Receiver Sensitivity (photons/bit)
% Output:
%   dr = Max datarate (bits/sec)
function [dr] = datarate(Pr, nu)
  dr = Pr/nu;
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
  %  angle_between_crafts = 2*acos(earth_radius/sc_radius);
  %  max_distance_between_sc =
  %  sin(angle_between_crafts/2)*(2*sc_radius);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Communication Simulator %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function communication_sim
sc_altitude  = 370; %km
earth_radius = 6378.1; %km
sc_radius    = sc_altitude + earth_radius; %km

[Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo()

%Altitude of weather/clouds varies from about 6 to 25 km.

% TODO: Take into account pointing accuracy loss
% TODO: Calculate power loss though atmosphere

%% We need to figure out the altitude for each portion of the
%% integration so we can find rho




end




%% The distance between the spacecraft determine how much the
%% signal is traveling through the atmosphere.
%% In order to get the altitude of the middle of the beam (Which is
%% the lowest point in the atmosphere that the beam will traval), a
%% simple geometric relation is made.

function [altitudes, valid] = list_altitudes(distance_between_sc, ...
                                             sc_radius, division)
  angle_between_crafts = 2*asin(distance_between_sc/(2*sc_radius));
  lowest_altitude = sc_radius * cos(angle_between_crafts / 2);
  for i = 1:division:
    
  end
end


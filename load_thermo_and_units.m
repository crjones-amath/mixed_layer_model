% load_thermo_and_units.m:  script to load thermodynamic path, constants,
% and common unit conversions
global g stefan r_earth Omega Rd Rv Cp Cpv Cw L delta psurf pref
addpath ../thermo/
here = pwd;

thermo_constants;               % load common thermo constants

gkg         = 1e-3;             % conversion for q (g/kg) -> (kg/kg)
kJkg        = 1e3;              % conversion for MSE kJ/kg -> J/kg
hPa         = 100;              % hPa -> Pa;
days_to_sec = 24*3600;          % days -> sec
mm_per_day  = 1e3.*days_to_sec; % (m/s) -> (mm/day)
mm_per_sec  = 1e3;              % (m/s) -> (mm/s)
g_per_m2    = 1e3;              % (kg/m2) -> (g/m2)
cm3         = 1e-6;             % (cm^3)  -> (m^3);
%---------------------------
% DYCOMS-specified constants
%---------------------------
Cp    = 1015;
L     = 2.47e6;
psurf = 1017.8.*hPa;

gam         = g./Cp;            % dry adiabatic lapse rate
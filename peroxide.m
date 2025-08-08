%% Front Matter

% Don't regen the database, but everything else should get cleared
% Otherwise weird state can persist in the mixture causing bad things
clearvars -except DB

import combustiontoolbox.databases.NasaDatabase;
import combustiontoolbox.rocket.*;
import combustiontoolbox.core.*;

DB = NasaDatabase();

sys = ChemicalSystem(DB);
solver = RocketSolver('problemType', 'ROCKET_IAC');

%% Design Conditions
p_amb = 0.2391; % 35,000 ft, standard day
T_amb = 218.9; %K
a_amb = 297; % m/s
M_cruise = 0.85;
L_D = 19;
eta_f = 0.95;
eta_shaft = 0.97;
dry_mass = 190000; % kg
prop_mass = 400000; % kg
tank_fraction = 0.1; % FIXME: Find this value
g = 9.81; % m/s^2
R_air = 287.05; %J/kgK
thrust_margin = 1.2;

% Rocket Turbine Design Parameters
pi_t = 1/30;
eta_t = 0.9;

%% Define mixture
mix = Mixture(sys);

%The problem we're going to solve here is part A of HW8 just as a test
fuel_mass_frac = 8; % percent
peroxide_strength = 70; % percent
chamber_pressure = 75; % bar
fuel_temp = 300; % K
fuel = {'RP_1'};
oxidizer = {'H2O2bLb', 'H2ObLb'};
fprintf("Propellant is %.1f%% RP-1 / %.1f%% High Test Peroxide\n", ...
    fuel_mass_frac, 100 - fuel_mass_frac);
fprintf("\t(Peroxide is %.1f%% strength)\n", peroxide_strength);

fuel_molar_weight = DB.getProperty(fuel, 'W');
oxidizer_molar_weight_vec = DB.getProperty(oxidizer, 'W');

oxidizer_molar_weight = (oxidizer_molar_weight_vec * [peroxide_strength
    (100 - peroxide_strength)]) / 100;

fuel_moles = fuel_mass_frac / fuel_molar_weight;
oxidizer_moles = [peroxide_strength (100-peroxide_strength)] * ...
    (100 - fuel_mass_frac) / (100 * oxidizer_molar_weight);

set(mix, fuel, 'fuel', fuel_moles);
set(mix, oxidizer, 'oxidizer', oxidizer_moles);
mix.config.compositionUnits = 'mass fraction';
print(mix);

%% Solve for the gas generator exit conditions
solver.FLAG_SUBSONIC = true;
mix.setTemperature(fuel_temp); %K
mix.setPressure(chamber_pressure); %bar
mix.areaRatio = 2; % Need to put something to keep the toolbox happy

fprintf('Input Conditions: %.1f K / %.1f bar\n', ...
    mix.T, mix.p, mix.areaRatio);

[inj_mix, chamber_mix, throat_mix, ~] = solver.solve(mix);

fprintf('Chamber Conditions are %.1f K / %.1f bar\n', ...
    chamber_mix.T, chamber_mix.p);

%% Turbine Analysis

tau_t = 1 + eta_t*(pi_t^((throat_mix.gamma - 1)/throat_mix.gamma) - 1);
display(tau_t);

% From Farokhi (2014), flow input to the turbine will be nearly choked,
% (M4 = 1) and a reasonable exit Mach number is 0.5 (M5 = 0.5)

M5 = 0.5;
M4 = 1;
T4 = throat_mix.T;
P4 = throat_mix.p;
gamma4 = throat_mix.gamma;
cp4 = throat_mix.cp / (throat_mix.N * throat_mix.MW);
throat_total = setStagnation(throat_mix);
Tt4 = throat_total.T;
Pt4 = throat_total.p;

print(throat_mix, throat_total);

PSFC_rkt = cp4 * Tt4 * (1 - tau_t);

fprintf('PSFC = %0.2f MW/kg/s\n', PSFC_rkt * 1e-6);
fprintf('     = %0.2f HP/lbm/hr\n', PSFC_rkt * 1e-6 * 0.169);

%% Calculate flow properties after the turbine

Tt5 = Tt4*tau_t;
Pt5 = Pt4*pi_t;
% Calculate gamma at Tt5
turbine_total = copy(throat_total);
turbine_total.setPressure(Pt5);
turbine_total.setTemperature(Tt5);
gamma5 = turbine_total.gamma;
fprintf('gamma_hot = %0.3f; gamma_cool = %0.3f\n', gamma4, gamma5);
turbine_exit_stagparam = 1 + ((gamma5 - 1)/2)*M5*M5;
T5 = Tt5 / turbine_exit_stagparam;
P5 = Pt5 / turbine_exit_stagparam;
turbine_exit = copy(turbine_total);
turbine_exit.setPressure(P5);
turbine_exit.setTemperature(T5);
turbine_exit.setProperties('mach', M5);
R5 = 8.314 / turbine_exit.MW;
print(turbine_exit, turbine_total);

[~, ~, ~, ~, Aconv] = flowisentropic(gamma5, M5);
[M6, T6, ~, ~, Adiv] = flowisentropic(gamma5, p_amb/Pt5, 'pres');
T6 = T6*Tt5;
Ve_rkt = M6 * sqrt(gamma5*R5*T6);

fprintf('Exit Nozzle %0.3f conv / %0.3f div\n', Aconv, Adiv);
fprintf('Exit Conditions %0.2f M / %d m/s\n', M6, Ve_rkt);

% Based on these exit conditions, we get approximately 15% of out total
% propulsive power from the rocket thrust

%% Fan Sizing

pi_f = 1.6; % ranges from 1.2 to 1.8 generally
gamma_air = 1.4;
tau_f = 1 + (pi_f^((gamma_air - 1)/gamma_air) - 1) / eta_f;

thrust_req_cruise = (dry_mass + prop_mass) * g / L_D;

cruise_stag_ratio = 1 + ((gamma_air-1)/2)*M_cruise*M_cruise;
Tt1 = T_amb*cruise_stag_ratio;
Pt1 = p_amb*(cruise_stag_ratio^(gamma_air/(gamma_air-1)));
Tt2 = Tt1*tau_f;
Pt2 = Pt1*pi_f;
rhot2 = (Pt2*1e5)/(R_air*Tt2); % convert bar to Pa for pressure
% We will assume constant axial velocity in this computation:
M2 = sqrt(1/( (Tt2/(M_cruise*M_cruise*T_amb)) - ((gamma_air - 1)/ 2)));
[~, T2, P2, ~, FanNozRatio] = flowisentropic(gamma_air, M2);
T2 = T2*Tt2;
P2 = P2*Pt2;
rho2 = P2*1e5/(R_air*T2);
airflow_per_fan_area = rho2*M2*sqrt(gamma_air*R_air*T2);
thrust_per_fan_area = 1e5*(P2 - p_amb);
cp_air = R_air*gamma_air/(gamma_air - 1);
power_per_fan_flow = cp_air*Tt1*(tau_f - 1);
power_per_fan_thrust = power_per_fan_flow * airflow_per_fan_area ...
    / thrust_per_fan_area;

fan_tsfc = PSFC_rkt * eta_shaft / power_per_fan_thrust;

overall_tsfc = fan_tsfc + Ve_rkt;

design_max_gross_thrust = thrust_req_cruise * thrust_margin;
design_max_fuelflow = design_max_gross_thrust / overall_tsfc;
design_max_fan_thrust = design_max_fuelflow * fan_tsfc;
fan_size = design_max_fan_thrust / thrust_per_fan_area;

fprintf('Fan sizing results in total fan area of %.1f m^2\n', fan_size);
% Though we used TSFC for ease of calculation, since fuel is primarily
% used to drive a turbine, we'll assume that PSFC is the constant value
PSFC_overall = design_max_gross_thrust * M_cruise * a_amb ...
    / design_max_fuelflow; % W/kg/s

%% Range calculation

GTOW = dry_mass + prop_mass;
max_fuel = prop_mass*(1-tank_fraction);
range = L_D*PSFC_overall*log(GTOW/(GTOW-max_fuel))/g; % meters
range_nmi = range / 1852;
fprintf('Calculated maximum range: %d nmi\n', round(range_nmi));

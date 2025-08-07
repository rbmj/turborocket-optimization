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

% Design Conditions
p_amb = 0.2391; % 35,000 ft, standard day
a_amb = 297; % m/s
M_cruise = 0.85;

%% Define mixture
mix = Mixture(sys);

%The problem we're going to solve here is part A of HW8 just as a test
fuel_mass_frac = 8; % percent
peroxide_strength = 70; % percent
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
mix.setTemperature(300); %K
mix.setPressure(75); %bar
mix.areaRatio = 2; % Need to put something to keep the toolbox happy

fprintf('Input Conditions: %.1f K / %.1f bar\n', ...
    mix.T, mix.p, mix.areaRatio);

[inj_mix, chamber_mix, throat_mix, ~] = solver.solve(mix);

fprintf('Chamber Conditions are %.1f K / %.1f bar\n', ...
    chamber_mix.T, chamber_mix.p);

%% Turbine Analysis

pi_t = 1/30;
eta_t = 0.9;
tau_t = 1 + eta_t*(pi_t^((throat_mix.gamma - 1)/throat_mix.gamma) - 1);
display(tau_t);

% From Farokhi (2014), flow input to the turbine will be nearly choked,
% and a reasonable exit Mach number is 0.5

T4 = throat_mix.T;
P4 = throat_mix.p;
gamma4 = throat_mix.gamma;
cp4 = throat_mix.cp / (throat_mix.N * throat_mix.MW);
throat_total = setStagnation(throat_mix);
Tt4 = throat_total.T;
Pt4 = throat_total.p;

print(throat_mix, throat_total);

PSFC = cp4 * Tt4 * (1 - tau_t);

fprintf('PSFC = %0.2f MW/kg/s\n', PSFC * 1e-6);
fprintf('     = %0.2f HP/lbm/hr\n', PSFC * 1e-6 * 0.169);

%% Calculate flow properties after the turbine

Tt5 = Tt4*tau_t;
Pt5 = Pt4*pi_t;
M5 = 0.5; % see above
% Use the approximation that Tt5 is close enough to T5 that gamma is const
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

%% Fan Computations

pi_f = 1.6; % Somewhat arbitrary, ranges from 1.2 to 1.8 generally
gamma_air = 1.4;



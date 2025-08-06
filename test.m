% Don't regen the database, but everything else should get cleared
% Otherwise weird state can persist in the mixture causing bad things
clearvars -except DB

import combustiontoolbox.databases.NasaDatabase;
import combustiontoolbox.rocket.*;
import combustiontoolbox.core.*;

DB = NasaDatabase();

sys = ChemicalSystem(DB);
solver = RocketSolver('problemType', 'ROCKET_IAC');
% Define mixture
% For some reason the toolbox uses bLb for liquid species
mix = Mixture(sys);

%The problem we're going to solve here is part A of HW8 just as a test
fuel_mass_frac = 50.75;

fuel = {'H2'};
oxidizer = {'O2bLb'};
fuel_molar_weight = DB.getProperty(fuel, 'W');
oxidizer_molar_weight = DB.getProperty(oxidizer, 'W');

set(mix, {'H2'}, 'fuel', fuel_mass_frac / fuel_molar_weight);
set(mix, {'O2bLb'}, 'oxidizer', (100 - fuel_mass_frac) / oxidizer_molar_weight);
mix = setTemperatureSpecies(mix, [200.1, 90]); % Same order as they were added
mix = setPressure(mix, 4967*0.0689476); % Bar
%Assume that the phase change enthalpy dominates and everything is at 90K
%setTemperature(mix, 90);
fprintf('Temperature set to %f\n', mix.T);
% Solver crashes if we don't specify an area ratio > 1:
mix.areaRatio = 2;

print(mix);

[inj_mix, chamber_mix, throat_mix, exit_mix] = solver.solve(mix);

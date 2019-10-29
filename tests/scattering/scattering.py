#!/usr/bin/env python3
# Copyright (C) 2019 Claire Hansel
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
# Copyright (C) 2019 Claire Hansel
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import simulation
import subprocess
import sys

def parameterDefaults(run_name):
    return {
        'rho_ion_initial_si': 1e26,
        'plasma_length_si': 100,
        'beam_energy_initial_gev': 10,
        'acceleration_gradient_gev_per_m': 0,
        'bennett_radius_initial_si': 2.5e-8,
        'cross_section_radius_si': 1e-8,
        'unperturbed_plasma_density_si': 1e24,
        'integration_tolerance': 1e-10,
        'vartheta_cutoff': 10,
        'ion_atomic_number': 1,
        'seed': 1563489705,
        'minimum_steps_per_betatron_period': 100,
        'particles_target': 100,
        'analysis_points_target': 100,
        'spline_points': 1000,
        'max_order': 3,
        'max_integration_depth': 15,
        'output_filename': 'data/{}_output'.format(run_name),
        'statistics_filename': 'data/{}_statistics'.format(run_name),
        'phase_space_filename': 'data/{}_phase_space'.format(run_name),
        'output_phase_space': True,
        'modified_bennett': True
    }

def run(run_name='run1', parameters={}, hoffman2=False, compute_processes=1):
    parameter_defaults = parameterDefaults(run_name)
    assert all(key in parameter_defaults for key in parameters)
    with open('data/{}_input'.format(run_name), 'w') as f:
        for key in parameter_defaults.keys():
            if key in parameters:
                value = parameters[key]
            else:
                value = parameter_defaults[key]
            f.write('{} = {}\n'.format(key, value))
    if hoffman2:
        subprocess.run(
            ['qsub', '-pe', 'dc_*', str(compute_processes + 1), 'run_epamss.sh', 'data/{}_input'.format(run_name)],
            check=True
        )
    else:
        subprocess.run(
            ['mpirun', '--oversubscribe', '-np', str(compute_processes + 1), 'epamss', 'data/{}_input'.format(run_name)],
            check=True
        )

def submitjobs():
    exponents = np.arange(-10, -5)
    for run_number, exponent in enumerate(exponents):
        runname = 'rsigma_{}'.format(exponent)
        run(runname, {'cross_section_radius_si': 10. ** exponent, 'minimum_steps_per_betatron_period': 100000}, hoffman2=True, compute_processes=63)
    exponents2 = np.arange(3, 7)
    for run_number, exponent in enumerate(exponents2):
        steps = 10 ** exponent
        runname = 'steps_{}'.format(run_number)
        run(runname, {'cross_section_radius_si': 1e-8, 'minimum_steps_per_betatron_period': 10 ** exponent}, hoffman2=True, compute_processes=63)

def analyze():
    sims1 = []
    labels1 = []
    exponents = np.arange(-10, -5)
    for run_number, exponent in enumerate(exponents):
        runname = 'rsigma_{}'.format(exponent)
        sims1.append(simulation.Simulation('data/{}_output'.format(runname)))
        labels1.append(r'$r_{\sigma} = 10^{' + str(exponent) + '}$')
    sims2 = []
    labels2 = []
    exponents2 = np.arange(3, 7)
    for run_number, exponent in enumerate(exponents2):
        steps = 10 ** exponent
        runname = 'steps_{}'.format(run_number)
        sims2.append(simulation.Simulation('data/{}_output'.format(runname)))
        labels2.append('$10^{' + str(exponent) + '}$ steps')

    simulation.Simulation.plotScatteringTest(sims1, labels1, filename='results/rsigma.png')
    simulation.Simulation.plotScatteringTest(sims2, labels2, filename='results/steps.png')

def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'run':
        submitjobs()
    elif len(sys.argv) > 1 and sys.argv[1] == 'analyze':
        analyze()

if __name__ == '__main__':
    main()

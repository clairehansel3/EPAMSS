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

rerun = True

def main1():
    compute_processes = 1
    sims = []
    labels = []
    exponents = np.arange(-10, -5)
    for run_number, exponent in enumerate(exponents):
        run = 'rsigma_{}'.format(run_number)
        cross_section_radius_si = 10. ** exponent
        if rerun:
            parameters = {
                'rho_ion_initial_si': 9.43718e23, # irrelevant
                'plasma_length_si': 100,
                'beam_energy_initial_gev': 10,
                'acceleration_gradient_gev_per_m': 0,
                'bennett_radius_initial_si': 2.37296e-7, # irrelevant
                'cross_section_radius_si': cross_section_radius_si,
                'unperturbed_plasma_density_si': 1e26,
                'integration_tolerance': 1e-10,
                'vartheta_cutoff': 10,
                'ion_atomic_number': 1,
                'minimum_steps_per_betatron_period': 100000, # steps
                'particles_target': 1000,
                'analysis_points_target': 1000,
                'spline_points': 1000,
                'max_order': 3,
                'max_integration_depth': 15,
                'output_filename': 'data/output_{}'.format(run),
                'statistics_filename': 'data/statistics_{}'.format(run),
                'phase_space_filename': 'data/phase_space_{}'.format(run),
                'output_phase_space': True,
                'modified_bennett': True # irrelevant
            }
            with open('data/input_{}'.format(run), 'w') as f:
                for key, value in parameters.items():
                    f.write('{} = {}\n'.format(key, value))
            subprocess.run(
                ['mpirun', '-np', str(compute_processes+1), 'epamss',
                'data/input_{}'.format(run)], check=True
            )
        sims.append(simulation.Simulation('data/output_{}'.format(run)))
        labels.append(r'$r_{\sigma} = 10^{' + str(exponent) + '}$')
    simulation.Simulation.plotScatteringTest(sims, labels, filename='results/rsigma.png')

def main2():
    compute_processes = 1
    sims = []
    labels = []
    exponents = np.arange(3, 7)
    for run_number, exponent in enumerate(exponents):
        steps = 10 ** exponent
        run = 'steps_{}'.format(run_number)
        if rerun:
            parameters = {
                'rho_ion_initial_si': 9.43718e23, # irrelevant
                'plasma_length_si': 100,
                'beam_energy_initial_gev': 10,
                'acceleration_gradient_gev_per_m': 0,
                'bennett_radius_initial_si': 2.37296e-7, # irrelevant
                'cross_section_radius_si': 1e-8,
                'unperturbed_plasma_density_si': 1e26,
                'integration_tolerance': 1e-10,
                'vartheta_cutoff': 10,
                'ion_atomic_number': 1,
                'minimum_steps_per_betatron_period': steps, # steps
                'particles_target': 1000,
                'analysis_points_target': 1000,
                'spline_points': 1000,
                'max_order': 3,
                'max_integration_depth': 15,
                'output_filename': 'data/output_{}'.format(run),
                'statistics_filename': 'data/statistics_{}'.format(run),
                'phase_space_filename': 'data/phase_space_{}'.format(run),
                'output_phase_space': True,
                'modified_bennett': True # irrelevant
            }
            with open('data/input_{}'.format(run), 'w') as f:
                for key, value in parameters.items():
                    f.write('{} = {}\n'.format(key, value))
            subprocess.run(
                ['mpirun', '-np', str(compute_processes+1), 'epamss',
                'data/input_{}'.format(run)], check=True
            )
        sims.append(simulation.Simulation('data/output_{}'.format(run)))
        labels.append('$10^{' + str(exponent) + '}$ steps')
    simulation.Simulation.plotScatteringTest(sims, labels, filename='results/steps.png')

if __name__ == '__main__':
    main1()
    main2()

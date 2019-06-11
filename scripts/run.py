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

import subprocess

Hoffman2 = False

def run():
    compute_processes = 1
    parameters = {
        'rho_ion_initial_si': 9.43718e23,
        'plasma_length_si': 1,
        'beam_energy_initial_gev': 1000,
        'acceleration_gradient_gev_per_m': 10,
        'bennett_radius_initial_si': 2.37296e-7,
        'cross_section_radius_si': 1e-8,
        'unperturbed_plasma_density_si': 9e21,
        'integration_tolerance': 1e-10,
        'vartheta_cutoff': 10,
        'ion_atomic_number': 1,
        'minimum_steps_per_betatron_period': 100,
        'particles_target': 10000,
        'analysis_points_target': 10,
        'spline_points': 1000,
        'max_order': 0,
        'max_integration_depth': 15,
        'output_filename': 'data/output',
        'statistics_filename': 'data/statistics',
        'phase_space_filename': 'data/phase_space',
        'output_phase_space': True,
        'modified_bennett': True
    }
    with open('data/input', 'w') as f:
        for key, value in parameters.items():
            f.write('{} = {}\n'.format(key, value))
    if Hoffman2:
        subprocess.run(
            ['qsub', '-pe', 'dc_*', str(compute_processes+1), 'run_epamss.sh'],
            check=True
        )
    else:
        subprocess.run(
            ['mpirun', '-np', str(compute_processes+1), 'epamss', 'data/input'],
            check=True
        )

if __name__ == '__main__':
    run()

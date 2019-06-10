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

Hoffman2 = True

def run():
    compute_processes = 63
    parameters = {
        'rho_ion_initial_si': 1e26,
        'plasma_length_si': 100,
        'beam_energy_initial_gev': 10,
        'acceleration_gradient_gev_per_m': 0,
        'bennett_radius_initial_si': 2.5e-8,
        'cross_section_radius_si': 1e-8,
        'unperturbed_plasma_density_si': 1e24,
        'integration_tolerance': 1e-15,
        'vartheta_cutoff': 10,
        'ion_atomic_number': 1,
        'minimum_steps_per_betatron_period': 100,
        'particles_target': 1000,
        'analysis_points_target': 1000,
        'spline_points': 1000,
        'max_order': 2,
        'max_integration_depth': 15,
        'output_filename': 'data/output',
        'statistics_filename': 'data/statistics',
        'phase_space_filename': 'data/phase_space',
        'output_phase_space': False,
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

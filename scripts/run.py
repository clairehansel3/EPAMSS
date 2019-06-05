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
    import os
    os.system('rm data/*')
    os.system('rm frames/*')
    os.system('rm results/*')
    compute_processes = 1
    parameters = {
        'rho_ion_initial_si': 1.3e26,
        'plasma_length_si': 1,
        'beam_energy_initial_gev': 8,
        'acceleration_gradient_gev_per_m': 8,
        'bennett_radius_initial_si': 2.5e-8,
        'cross_section_radius_si': 1e-8,
        'unperturbed_plasma_density_si': 9.3e23,
        'integration_tolerance': 1e-15,
        'vartheta_cutoff': 11,
        'ion_atomic_number': 1,
        'minimum_steps_per_betatron_period': 100,
        'particles_target': 53040,
        'analysis_points_target': 100,
        'spline_points': 952,
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
    import analyze
    analyze.analyze()

if __name__ == '__main__':
    run()

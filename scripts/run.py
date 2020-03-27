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

import matplotlib.pyplot as plt
import numpy as np
import simulation
import subprocess
import os.path

def parameterDefaults(run_name):
    return {
        'rho_ion_si': 1e26,
        'plasma_length_si': 1,
        'beam_energy_initial_gev': 10,
        'acceleration_gradient_gev_per_m': 10,
        'bennett_radius_initial_si': 170e-9,
        'cross_section_radius_si': 1e-8,
        'unperturbed_plasma_density_si': 1e24,
        'integration_tolerance': 1e-10,
        'vartheta_cutoff': 10,
        'drive_amplitude': 0,
        'drive_angular_frequency':0,
        'factor': 0,
        'ion_atomic_number': 1,
        'minimum_steps_per_betatron_period': 200,
        'particles_target': 50000,
        'analysis_points_target': 1000,
        'spline_points': 1000,
        'max_order': 3,
        'max_integration_depth': 15,
        'output_filename': 'data/{}_output'.format(run_name),
        'statistics_filename': 'data/{}_statistics'.format(run_name),
        'phase_space_filename': 'data/{}_phase_space'.format(run_name),
        'output_phase_space': False,
        'modified_bennett': True,
        'scattering2': True
    }

def run(run_name='default', parameters={}, hoffman2=False, compute_processes=1):
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
            ['mpirun', '-np', str(compute_processes + 1), 'epamss', 'data/{}_input'.format(run_name)],
            check=True
        )
    if os.path.isfile('data/{}_output'.format(run_name)):
        return simulation.Simulation(parameters['output_filename'] if 'output_filename' in parameters else parameter_defaults['output_filename'])


def main():
    s = run(run_name='facet2', hoffman2=True, compute_processes=63)
    if s is not None:
        s.plot2DEmittance()

if __name__ == '__main__':
    main()

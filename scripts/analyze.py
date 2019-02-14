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
import numpy.linalg as la
import matplotlib.pyplot as plt

parameter_types = {
    'maximum_ion_density': float,
    'plasma_length': float,
    'beam_energy': float,
    'electron_linear_density': float,
    'bennett_radius': float,
    'interaction_radius': float,
    'integration_tolerance': float,
    'vartheta_cutoff': float,
    'ion_atomic_number': int,
    'minimum_steps_per_betatron_period': int,
    'particles': int,
    'analysis_points_target': int,
    'spline_points': int,
    'seed': int,
    'max_order': int,
    'max_integration_depth': int,
    'output_filename': str,
    'statistics_filename': str,
    'phase_space_filename': str,
    'output_phase_space': bool,
    'ion_linear_density': float,
    'gamma': float,
    'alpha': float,
    'step_size': float,
    'cross_section': float,
    'minimum_angle': float,
    'betatron_frequency': float,
    'betatron_period': float,
    'max_scattering_r_div_a': float,
    'percent_with_scattering': float,
    'omega_on_axis': float,
    'steps': int,
    'stride': int,
    'actual_particles': int,
    'particles_per_process': int,
    'compute_processes': int,
    'analysis_points': int,
    'seconds_elapsed': float,
    'minutes_elapsed': float,
    'hours_elapsed': float,
    'approx_core_hours_elapsed': float
}

def getOutput(path):
    output = {}
    with open(path, 'r') as f:
        for line in f:
            key, _, value = line.split()
            output[key] = parameter_types[key](value)
    return output

def getZ(output):
    return np.linspace(0, output['plasma_length'], output['analysis_points'])

def getStatistics(output):
    raw_data = np.memmap(
        output['statistics_filename'],
        dtype=np.float64,
        mode='r',
        shape=(2, output['analysis_points'], 14)
    )
    means = np.empty((2, 4, output['analysis_points']))
    covariance_matrix = np.empty((2, 4, 4, output['analysis_points']))
    for i in (0, 1):
        means[i, 0, :] = raw_data[i, :, 0]
        means[i, 1, :] = raw_data[i, :, 1]
        means[i, 2, :] = raw_data[i, :, 2]
        means[i, 3, :] = raw_data[i, :, 3]
        covariance_matrix[i, 0, 0, :] = raw_data[i, :, 4]
        covariance_matrix[i, 0, 1, :] = covariance_matrix[i, 1, 0, :] = raw_data[i, :, 5]
        covariance_matrix[i, 0, 2, :] = covariance_matrix[i, 2, 0, :] = raw_data[i, :, 6]
        covariance_matrix[i, 0, 3, :] = covariance_matrix[i, 3, 0, :] = raw_data[i, :, 7]
        covariance_matrix[i, 1, 1, :] = raw_data[i, :, 8]
        covariance_matrix[i, 1, 2, :] = covariance_matrix[i, 2, 1, :] = raw_data[i, :, 9]
        covariance_matrix[i, 1, 3, :] = covariance_matrix[i, 3, 1, :] = raw_data[i, :, 10]
        covariance_matrix[i, 2, 2, :] = raw_data[i, :, 11]
        covariance_matrix[i, 2, 3, :] = covariance_matrix[i, 3, 2, :] = raw_data[i, :, 12]
        covariance_matrix[i, 3, 3, :] = raw_data[i, :, 13]
    return means, covariance_matrix

def getPhaseSpace(output):
    assert output['output_phase_space']
    phase_space = np.empty((2, 4, output['actual_particles'], output['analysis_points']))
    for i in range(output['compute_processes']):
        raw_data = np.memmap(
            '{}_{}'.format(output['phase_space_filename'], i + 1),
            dtype=np.float64,
            mode='r',
            shape=(2, output['analysis_points'], output['particles_per_process'], 4)
        )
        for j in (0, 1):
            for k in range(4):
                for step in range(output['analysis_points']):
                    phase_space[j, k, i*output['particles_per_process']:
                        (i+1)*output['particles_per_process'], step] = \
                        raw_data[j, step, :, k]
    return phase_space


def get2DEmittance(covariance_matrix):
    emittance = np.empty((2, 2, covariance_matrix.shape[3]))
    for i in (0, 1):
        for step in range(covariance_matrix.shape[3]):
            emittance[i, 0, step] = np.sqrt(la.det(covariance_matrix[i, 0:2, 0:2, step]))
            emittance[i, 1, step] = np.sqrt(la.det(covariance_matrix[i, 2:4, 2:4, step]))
    return emittance

def get4DEmittance(covariance_matrix):
    emittance = np.empty((2, covariance_matrix.shape[3]))
    for i in (0, 1):
        for step in range(covariance_matrix.shape[3]):
            emittance[i, step] = np.sqrt(la.det(covariance_matrix[i, :, :, step]))
    return emittance

def computeEmittance(x, vx):
    assert len(x) == len(vx)
    N = len(x)
    x_mean = 0.0
    vx_mean = 0.0
    for i in range(N):
        x_mean += x[i]
        vx_mean += vx[i]
    x_mean /= N
    vx_mean /= N
    m_x = 0.0
    m_vx = 0.0
    c_x = 0.0
    for i in range(N):
        delta_x = x[i] - x_mean
        delta_vx = vx[i] - vx_mean
        m_x += delta_x * delta_x
        m_vx += delta_vx * delta_vx
        c_x += delta_x * delta_vx
    val_x2 = m_x / N
    val_vx2 = m_vx / N
    val_xvx = c_x / N
    x_emit = np.sqrt(val_x2 * val_vx2 - val_xvx * val_xvx)
    return x_emit

def getEmittanceFromPhaseSpace(phase_space):
    emit = [[[], []], [[], []]]
    for i in (0, 1):
        for step in range(phase_space.shape[3]):
            emit[i][0].append(computeEmittance(phase_space[i, 0, :, step], phase_space[i, 1, :, step]))
            emit[i][1].append(computeEmittance(phase_space[i, 2, :, step], phase_space[i, 3, :, step]))
    return np.array(emit)

def getEnergyFromPhaseSpace(output, phase_space):
    energy = np.empty((2, output['particles'], output['analysis_points']))
    for i in (0, 1):
        for particle in range(output['actual_particles']):
            energy[i, particle, :] = phase_space[i, 1, particle, :] ** 2 + \
                phase_space[i, 3, particle, :] ** 2 + np.log(1 + (( \
                phase_space[i, 0, particle, :] ** 2 + \
                phase_space[i, 2, particle, :] ** 2) / ( \
                output['bennett_radius'] ** 2)))
    return energy

def plot4DEmittance(output, z, emit_4d):
    k_beta = np.sqrt(output['alpha']) / output['bennett_radius']
    plt.title('Emittance Growth')
    plt.plot(k_beta * z, 100 * ((emit_4d[0, :] / emit_4d[0, 0]) - 1), label='no scattering', color='blue')
    plt.plot(k_beta * z, 100 * ((emit_4d[1, :] / emit_4d[0, 0]) - 1), label='scattering', color='red')
    plt.xlabel(r'$k_{\beta} z$ [unitless]')
    plt.ylabel(r'Percent RMS 4D Emittance Growth [unitless]')
    plt.legend()
    plt.savefig('results/emit_4d.png')
    plt.cla()

def analyze():
    output = getOutput('data/output')
    z = getZ(output)
    means, covariance_matrix = getStatistics(output)
    emittance_2d = get2DEmittance(covariance_matrix)
    emittance_4d = get4DEmittance(covariance_matrix)
    if output['output_phase_space']:
        phase_space = getPhaseSpace(output)
        phase_space_emit = getEmittanceFromPhaseSpace(phase_space)
    plot4DEmittance(output, z, emittance_4d)

if __name__ == '__main__':
    analyze()

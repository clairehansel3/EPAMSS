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
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.signal

def getOutputDict(output_filename):
    """
    Takes the path to the output file as an argument and returns a dict
    containing the values of the parameters in the output file converted to the
    appropriate types (e.g. 'maximum_ion_density' is a python float and
    'max_integration_depth' is a python int).
    """
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
        'particles_target': int,
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
        'actual_analysis_points': int,
        'seconds_elapsed': float,
        'minutes_elapsed': float,
        'hours_elapsed': float,
        'approx_core_hours_elapsed': float
    }
    output_dict = {}
    with open(output_filename, 'r') as f:
        for line in f:
            key, _, value = line.split()
            if parameter_types[key] is bool:
                if value in ('true', 'True'):
                    output_dict[key] = True
                elif value in ('false', 'False'):
                    output_dict[key] = False
                else:
                    assert False
            else:
                output_dict[key] = parameter_types[key](value)
    return output_dict

def getZ(output_dict):
    """ Returns an array of the z positions of each analysis point """
    return np.linspace(0, output_dict['plasma_length'], output_dict['actual_analysis_points'])

def getStatistics(output_dict):
    """
    Reads the statistics file and returns the two arrays mean and
    covariance_matrix. For both arrays, the first dimension has size 2 and
    indicates whether scattering is enabled with 0 being without scattering and
    1 being with scattering. For the mean array, the second dimension has length
    4 and indicates which value the mean is of, with 0 being x, 1 being x', 2
    being y, and 3 being y'. For the covariance_matrix array, the second and
    third dimensions both have length 4 and together they indicate which element
    of the covariance matrix is wanted, e.g. (1, 2) corresponds to Cov(x', y).
    Finally for both arrays the last dimension indicates at which analysis point
    the data were taken at and has length output_dict['actual_analysis_points'].
    """
    raw_data = np.memmap(
        output_dict['statistics_filename'],
        dtype=np.float64,
        mode='r',
        shape=(2, output_dict['actual_analysis_points'], 14)
    )
    means = np.delete(np.swapaxes(raw_data, 1, 2), np.s_[4:], axis=1)
    covariance_matrix = np.empty((2, 4, 4, output_dict['actual_analysis_points']))
    covariance_matrix[:, 0, 0, :] = raw_data[:, :, 4]
    covariance_matrix[:, 0, 1, :] = covariance_matrix[:, 1, 0, :] = raw_data[:, :, 5]
    covariance_matrix[:, 0, 2, :] = covariance_matrix[:, 2, 0, :] = raw_data[:, :, 6]
    covariance_matrix[:, 0, 3, :] = covariance_matrix[:, 3, 0, :] = raw_data[:, :, 7]
    covariance_matrix[:, 1, 1, :] = raw_data[:, :, 8]
    covariance_matrix[:, 1, 2, :] = covariance_matrix[:, 2, 1, :] = raw_data[:, :, 9]
    covariance_matrix[:, 1, 3, :] = covariance_matrix[:, 3, 1, :] = raw_data[:, :, 10 ]
    covariance_matrix[:, 2, 2, :] = raw_data[:, :, 11]
    covariance_matrix[:, 2, 3, :] = covariance_matrix[:, 3, 2, :] = raw_data[:, :, 12]
    covariance_matrix[:, 3, 3, :] = raw_data[:, :, 13]
    return means, covariance_matrix

def getPhaseSpace(output_dict):
    """
    Reads the phase space files and returns an array containing the phase space
    data. The first dimension of the result has size 2 and indicates whether
    scattering is enabled with 0 being without scattering and 1 being with
    scattering. The second dimension has size 4 and indicates which phase space
    coordinate is wanted with 0 being x, 1 being x', 2 begin y, and 3 being y'.
    The third dimension has size output_dict['actual_particles'] and indicates
    the number of the particle. The fourth and last dimension has size
    output_dict['actual_analysis_points'] and indicates at which analysis point
    the data were taken.
    """
    assert output_dict['output_phase_space']
    phase_space = np.empty((2, 4, output_dict['actual_particles'], output_dict['actual_analysis_points']))
    for i in range(output_dict['compute_processes']):
        raw_data = np.memmap(
            '{}_{}'.format(output_dict['phase_space_filename'], i + 1),
            dtype=np.float64,
            mode='r',
            shape=(2, output_dict['actual_analysis_points'], output_dict['particles_per_process'], 4)
        )
        particle_start = i * output_dict['particles_per_process']
        particle_end = (i + 1) * output_dict['particles_per_process']
        phase_space[:, :, particle_start:particle_end, :] = np.transpose(raw_data, (0, 3, 2, 1))
    return phase_space

def get2DEmittances(covariance_matrix):
    """
    Computes the 2D emittances from the scattering matrix. The first dimension
    of the result has size 2 and indicates whether scattering is enabled with 0
    being without scattering and 1 being with scattering. The second dimension
    has size 2 and indicates which dimension is wanted with 0 being the x
    emittance and 1 being the y emittance. The third and last dimension has size
    output_dict['actual_analysis_points'] and indicates at which analysis point
    the data were taken.
    """
    emittances_2d = np.empty((2, 2, covariance_matrix.shape[3]))
    emittances_2d[:, 0, :] = np.sqrt(np.linalg.det(np.transpose(covariance_matrix[:, 0:2, 0:2, :], (0, 3, 1, 2))))
    emittances_2d[:, 1, :] = np.sqrt(np.linalg.det(np.transpose(covariance_matrix[:, 2:4, 2:4, :], (0, 3, 1, 2))))
    return emittances_2d

def get4DEmittance(covariance_matrix):
    """
    Computes the 4D transverse emittance from the scattering matrix. The first
    dimension of the result has size 2 and indicates whether scattering is
    enabled with 0 being without scattering and 1 being with scattering. The
    second dimension has size output_dict['actual_analysis_points'] and
    indicates at which analysis point the data were taken.
    """
    return np.sqrt(np.linalg.det(np.transpose(covariance_matrix, (0, 3, 1, 2))))

def getAverageEnergy(output_dict, phase_space):
    x = phase_space[:, 0, :, :]
    vx = phase_space[:, 1, :, :]
    y = phase_space[:, 2, :, :]
    vy = phase_space[:, 3, :, :]
    energy = vx ** 2 + vy ** 2 + output_dict['alpha'] * np.log(
        1 + ((x ** 2 + y ** 2) / (output_dict['bennett_radius'] ** 2)))
    avg_energy = np.mean(energy, axis=1)
    return avg_energy

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

def plotAverageEnergy(output_dict, z, avg_energy):
    k_beta = np.sqrt(output_dict['alpha']) / output_dict['bennett_radius']
    plt.title('Average Energy Growth')
    plt.plot(k_beta * z, 100 * ((avg_energy[0] / avg_energy[0, 0]) - 1), label='no scattering', color='blue')
    plt.plot(k_beta * z, 100 * ((avg_energy[1] / avg_energy[1, 0]) - 1), label='scattering', color='red')
    plt.xlabel(r'$k_{\beta} z$ [unitless]')
    plt.xlabel(r'Percent Average Energy Growth [unitless]')
    plt.legend()
    plt.savefig('results/avg_energy.png')
    plt.cla()

def plot4DEmittanceGrowth(output_dict, z, emittance_4d, smoothing_on, window_size=101, order=3):
    k_beta = np.sqrt(output_dict['alpha']) / output_dict['bennett_radius']
    ns_emit = np.copy(emittance_4d[0, :])
    s_emit = np.copy(emittance_4d[1, :])
    z2 = np.copy(z)
    avg_emit = np.mean(ns_emit)
    ns_emit /= avg_emit
    s_emit /= avg_emit
    if smoothing_on:
        ns_emit = scipy.signal.savgol_filter(ns_emit, window_size, order)
        s_emit = scipy.signal.savgol_filter(s_emit, window_size, order)
        ns_emit = ns_emit[window_size:-window_size]
        s_emit = s_emit[window_size:-window_size]
        z2 = z2[window_size:-window_size]
        plt.title('4D RMS Emittance Growth (smoothed, window_size={}, order={})'.format(window_size, order))
    else:
        plt.title('4D RMS Emittance Growth')
    plt.plot(k_beta * z2, ns_emit, label='no scattering', color='blue')
    plt.plot(k_beta * z2, s_emit, label='scattering', color='red')
    plt.xlim(0, k_beta * output_dict['plasma_length'])
    plt.xlabel(r'$k_{\beta} z$ [unitless]')
    plt.ylabel(r'Normalized $\epsilon_{4D}(z)$ [unitless]')
    plt.legend()
    if smoothing_on:
        plt.savefig('results/emit_4d_window_{}_order_{}.png'.format(window_size, order))
    else:
        plt.savefig('results/emit_4d.png')
    plt.cla()


def plot2DEmittanceGrowth(output_dict, z, emittances_2d, smoothing_on, window_size=101, order=3):
    k_beta = np.sqrt(output_dict['alpha']) / output_dict['bennett_radius']
    ns_x_emit = np.copy(emittances_2d[0, 0, :])
    ns_y_emit = np.copy(emittances_2d[0, 1, :])
    s_x_emit = np.copy(emittances_2d[1, 0, :])
    s_y_emit = np.copy(emittances_2d[1, 1, :])
    z2 = np.copy(z)
    avg_emit = np.mean(np.concatenate((ns_x_emit, ns_y_emit)))
    ns_x_emit /= avg_emit
    ns_y_emit /= avg_emit
    s_x_emit /= avg_emit
    s_y_emit /= avg_emit
    if smoothing_on:
        ns_x_emit = scipy.signal.savgol_filter(ns_x_emit, window_size, order)
        ns_y_emit = scipy.signal.savgol_filter(ns_y_emit, window_size, order)
        s_x_emit = scipy.signal.savgol_filter(s_x_emit, window_size, order)
        s_y_emit = scipy.signal.savgol_filter(s_y_emit, window_size, order)
        ns_x_emit = ns_x_emit[window_size:-window_size]
        ns_y_emit = ns_y_emit[window_size:-window_size]
        s_x_emit = s_x_emit[window_size:-window_size]
        s_y_emit = s_y_emit[window_size:-window_size]
        z = z[window_size:-window_size]
        plt.title('2D RMS Emittance Growth (smoothed, window_size={}, order={})'.format(window_size, order))
    else:
        plt.title('2D RMS Emittance Growth')
    plt.plot(k_beta * z, ns_x_emit, label='x (no scattering)', color='blue')
    plt.plot(k_beta * z, ns_y_emit, label='y (no scattering)', color='green')
    plt.plot(k_beta * z, s_x_emit, label='x (scattering)', color='red')
    plt.plot(k_beta * z, s_y_emit, label='y (scattering)', color='orange')
    plt.xlim(0, k_beta * output_dict['plasma_length'])
    plt.xlabel(r'$k_{\beta} z$ [unitless]')
    plt.ylabel(r'Normalized $\epsilon_{2D}(z)$ [unitless]')
    plt.legend()
    if smoothing_on:
        plt.savefig('results/emit_2d_window_{}_order_{}.png'.format(window_size, order))
    else:
        plt.savefig('results/emit_2d.png')
    plt.cla()

def generate2DEmittanceGrowthRates(z, emittances_2d):
    s_x_emit = np.copy(emittances_2d[1, 0, :])
    s_y_emit = np.copy(emittances_2d[1, 1, :])
    s_x_emit /= s_x_emit[0]
    s_y_emit /= s_y_emit[0]
    slope_x, intercept_x, r_value_x, p_value_x, std_err_x = scipy.stats.linregress(z, s_x_emit - 1);
    slope_y, intercept_y, r_value_y, p_value_y, std_err_y = scipy.stats.linregress(z, s_y_emit - 1);
    with open('results/growth_rate.txt', 'w') as f:
        f.write('g_x = {} +/- {} [%/m]\n'.format(100 * slope_x, 100 * std_err_x))
        f.write('g_y = {} +/- {} [%/m]\n'.format(100 * slope_y, 100 * std_err_y))

def analyze():
    output_dict = getOutputDict('data/output')
    z = getZ(output_dict)
    means, covariance_matrix = getStatistics(output_dict)
    emittance_4d = get4DEmittance(covariance_matrix)
    emittances_2d = get2DEmittances(covariance_matrix)
    plot4DEmittanceGrowth(output_dict, z, emittance_4d, False)
    plot4DEmittanceGrowth(output_dict, z, emittance_4d, True)
    plot2DEmittanceGrowth(output_dict, z, emittances_2d, False)
    plot2DEmittanceGrowth(output_dict, z, emittances_2d, True)
    generate2DEmittanceGrowthRates(z, emittances_2d)

if __name__ == '__main__':
    analyze()

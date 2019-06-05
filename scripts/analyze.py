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
import scipy.linalg
import plot

def getOutputDict(output_filename):
    """
    Takes the path to the output file as an argument and returns a dict
    containing the values of the parameters in the output file converted to the
    appropriate types (e.g. 'maximum_ion_density' is a python float and
    'max_integration_depth' is a python int).
    """
    parameter_types = {
        'rho_ion_initial_si': float,
        'plasma_length_si': float,
        'beam_energy_initial_gev': float,
        'acceleration_gradient_gev_per_m': float,
        'bennett_radius_initial_si': float,
        'cross_section_radius_si': float,
        'unperturbed_plasma_density_si': float,
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
        'modified_bennett': bool,
        'plasma_frequency_si': float,
        'plasma_angular_wavenumber_si': float,
        'plasma_skin_depth_si': float,
        'rho_ion_initial': float,
        'plasma_length': float,
        'gamma_initial': float,
        'gamma_prime': float,
        'bennett_radius_initial': float,
        'cross_section_radius': float,
        'delta': float,
        'gamma_final': float,
        'bennett_radius_final': float,
        'rho_ion_final': float,
        'betatron_frequency_final': float,
        'betatron_period_final': float,
        'betatron_frequency_final_si': float,
        'betatron_period_final_si': float,
        'step_size': float,
        'step_size_si': float,
        'gamma_minimum_angle': float,
        'omega_off_axis': float,
        'omega_on_axis_initial': float,
        'max_scattering_r_div_a_initial': float,
        'sigma_r_initial': float,
        'sigma_r_prime_initial': float,
        'steps': int,
        'stride': int,
        'compute_processes': int,
        'actual_particles': int,
        'particles_per_process': int,
        'actual_analysis_points': int,
        'seconds_elapsed': float,
        'minutes_elapsed': float,
        'hours_elapsed': float,
        'approx_core_hours_elapsed': float
    }
    output_dict = {}
    with open(output_filename, 'r') as f:
        for line in f:
            key, _, value, *_ = line.split()
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
        import os
        print(os.path.getsize('{}_{}'.format(output_dict['phase_space_filename'], i + 1)))
        print(2 * output_dict['actual_analysis_points'] * output_dict['particles_per_process'] * 4)
        print(8 * 2 * output_dict['actual_analysis_points'] * output_dict['particles_per_process'] * 4)
        raw_data = np.memmap(
            '{}_{}'.format(output_dict['phase_space_filename'], i + 1),
            dtype=np.float64,
            mode='r',
            shape=(2, output_dict['actual_analysis_points'], output_dict['particles_per_process'], 4)
        )
        particle_start = i * output_dict['particles_per_process']
        particle_end = (i + 1) * output_dict['particles_per_process']
        print('i = ', i, ', transposing')
        phase_space[:, :, particle_start:particle_end, :] = np.transpose(raw_data, (0, 3, 2, 1))
        print('done')
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

def getGamma(output_dict, z):
    return output_dict['gamma_initial'] + output_dict['gamma_prime'] * z

def getBennettRadius(output_dict, gamma):
    return output_dict['bennett_radius_initial'] * ((gamma / output_dict['gamma_initial']) ** (-1/4))

def getRhoIon(output_dict, gamma):
    return output_dict['rho_ion_initial'] * np.sqrt(gamma / output_dict['gamma_initial'])

def getEnergy(output_dict, phase_space, gamma, rho_ion, bennett_radius):
    x = phase_space[:, 0, :, :]
    px = phase_space[:, 1, :, :]
    y = phase_space[:, 2, :, :]
    py = phase_space[:, 3, :, :]
    a2 = bennett_radius * bennett_radius
    x2 = x * x
    y2 = y * y
    energy = (px * px + py * py) / (2 * gamma) + \
        (output_dict['ion_atomic_number'] / 4) * (x2 + y2 + rho_ion * a2 * np.log(a2 + x2 + y2))
    return energy

def computeEmittance(x, px):
    assert len(x) == len(px)
    N = len(x)
    x_mean = 0.0
    px_mean = 0.0
    for i in range(N):
        x_mean += x[i]
        px_mean += px[i]
    x_mean /= N
    px_mean /= N
    m_x = 0.0
    m_px = 0.0
    c_x = 0.0
    for i in range(N):
        delta_x = x[i] - x_mean
        delta_px = px[i] - px_mean
        m_x += delta_x * delta_x
        m_px += delta_px * delta_px
        c_x += delta_x * delta_px
    val_x2 = m_x / N
    val_px2 = m_px / N
    val_xpx = c_x / N
    x_emit = np.sqrt(val_x2 * val_px2 - val_xpx * val_xpx)
    return x_emit

def getEmittanceFromPhaseSpace(phase_space):
    emit = [[[], []], [[], []]]
    for i in (0, 1):
        for step in range(phase_space.shape[3]):
            emit[i][0].append(computeEmittance(phase_space[i, 0, :, step], phase_space[i, 1, :, step]))
            emit[i][1].append(computeEmittance(phase_space[i, 2, :, step], phase_space[i, 3, :, step]))
    return np.array(emit)

def analyze():
    output_dict = getOutputDict('data/output')
    phase_space = getPhaseSpace(output_dict)
    z = getZ(output_dict)
    gamma = getGamma(output_dict, z)
    bennett_radius = getBennettRadius(output_dict, gamma)
    rho_ion = getRhoIon(output_dict, gamma)
    means, covariance_matrix = getStatistics(output_dict)
    emittance_4d = get4DEmittance(covariance_matrix)
    emittances_2d = get2DEmittances(covariance_matrix)
    energy = getEnergy(output_dict, phase_space, gamma, rho_ion, bennett_radius)
    #plot.plotBennettParameters(output_dict, z, gamma, rho_ion, bennett_radius)
    plot.plotTraj(output_dict, z, phase_space)
    #plot.createHistogramMovie(output_dict, z, phase_space, bennett_radius, rho_ion, gamma)
    #plot.plotEnergy(output_dict, z, energy)
    #plot.plot4DEmittanceGrowth(output_dict, z, emittance_4d, False)
    #plot.plot2DEmittanceGrowth(output_dict, z, emittances_2d, False)
    #plot.generate2DEmittanceGrowthRates(output_dict, z, emittances_2d)

if __name__ == '__main__':
    analyze()

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
import os
import scipy.integrate

def getOutputDict(output_filename):
    '''
    Takes the path to the output file as an argument and returns a dict
    containing the values of the parameters in the output file converted to the
    appropriate types (e.g. 'plasma_length_si' is a python float and
    'ion_atomic_number' is a python int).
    '''
    parameter_types = {
        'rho_ion_si': float,
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
        'rho_ion_div_n0': float,
        'plasma_length': float,
        'gamma_initial': float,
        'gamma_prime': float,
        'bennett_radius_initial': float,
        'cross_section_radius': float,
        'delta': float,
        'gamma_final': float,
        'bennett_radius_final': float,
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
        'sigma_r': float,
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


def getStatistics(output_dict):
    '''
    Reads the statistics file and returns the two arrays mean and
    covariance_matrix. For both arrays, the first dimension has size 2 and
    indicates whether scattering is enabled with 0 being without scattering and
    1 being with scattering. For the mean array, the second dimension has length
    4 and indicates which value the mean is of, with 0 being x, 1 being px, 2
    being y, and 3 being py. For the covariance_matrix array, the second and
    third dimensions both have length 4 and together they indicate which element
    of the covariance matrix is wanted, e.g. (1, 2) corresponds to Cov(px, y).
    Finally for both arrays the last dimension indicates at which analysis point
    the data were taken at and has length output_dict['actual_analysis_points'].

    means: [scattering, coordinate, z]
    covariance_matrix: [scattering, coordinate, coordinate, z]
    '''
    raw_data = np.memmap(
        output_dict['statistics_filename'],
        dtype=np.float64,
        mode='r',
        shape=(2, output_dict['actual_analysis_points'], 14)
    )
    means = np.delete(np.swapaxes(raw_data, 1, 2), np.s_[4:], axis=1)
    covariance_matrix = np.empty((2, 4, 4,
        output_dict['actual_analysis_points']))
    covariance_matrix[:, 0, 0, :] = raw_data[:, :, 4]
    covariance_matrix[:, 0, 1, :] = raw_data[:, :, 5]
    covariance_matrix[:, 1, 0, :] = raw_data[:, :, 5]
    covariance_matrix[:, 0, 2, :] = raw_data[:, :, 6]
    covariance_matrix[:, 2, 0, :] = raw_data[:, :, 6]
    covariance_matrix[:, 0, 3, :] = raw_data[:, :, 7]
    covariance_matrix[:, 3, 0, :] = raw_data[:, :, 7]
    covariance_matrix[:, 1, 1, :] = raw_data[:, :, 8]
    covariance_matrix[:, 1, 2, :] = raw_data[:, :, 9]
    covariance_matrix[:, 2, 1, :] = raw_data[:, :, 9]
    covariance_matrix[:, 1, 3, :] = raw_data[:, :, 10]
    covariance_matrix[:, 3, 1, :] = raw_data[:, :, 10]
    covariance_matrix[:, 2, 2, :] = raw_data[:, :, 11]
    covariance_matrix[:, 2, 3, :] = raw_data[:, :, 12]
    covariance_matrix[:, 3, 2, :] = raw_data[:, :, 12]
    covariance_matrix[:, 3, 3, :] = raw_data[:, :, 13]
    return means, covariance_matrix

def getPhaseSpace(output_dict, gamma):
    '''
    Reads the phase space files and returns an array containing the phase space
    data. The first dimension of the result has size 2 and indicates whether
    scattering is enabled with 0 being without scattering and 1 being with
    scattering. The second dimension has size 4 and indicates which phase space
    coordinate is wanted with 0 being x, 1 being px, 2 begin y, and 3 being py.
    The third dimension has size output_dict['actual_particles'] and indicates
    the number of the particle. The fourth and last dimension has size
    output_dict['actual_analysis_points'] and indicates at which analysis point
    the data were taken.

    [scattering, coordinate, particle, z]
    '''
    assert output_dict['output_phase_space']
    phase_space = np.empty((2, 4, output_dict['actual_particles'],
        output_dict['actual_analysis_points']))
    for i in range(output_dict['compute_processes']):
        raw_data = np.memmap(
            '{}_{}'.format(output_dict['phase_space_filename'], i + 1),
            dtype=np.float64,
            mode='r',
            shape=(2, output_dict['actual_analysis_points'],
                   output_dict['particles_per_process'], 4)
        )
        particle_start = i * output_dict['particles_per_process']
        particle_end = (i + 1) * output_dict['particles_per_process']
        print('i = ', i, ', transposing')
        phase_space[:, :, particle_start:particle_end, :] = np.transpose(
            raw_data, (0, 3, 2, 1))
        print('done transposing')
    # undo coordinate transform
    sqrt_gamma = np.sqrt(gamma)
    for i in range(2):
        for j in range(output_dict['actual_particles']):
            phase_space[i, 0, j] /= sqrt_gamma
            phase_space[i, 1, j] = sqrt_gamma * phase_space[i, 1, j] - 0.5 * \
                phase_space[i, 0, j] * output_dict['gamma_prime']
            phase_space[i, 2, j] /= sqrt_gamma
            phase_space[i, 3, j] = sqrt_gamma * phase_space[i, 3, j] - 0.5 * \
                phase_space[i, 2, j] * output_dict['gamma_prime']
    return phase_space

def get2DEmittances(covariance_matrix):
    '''
    Computes the 2D emittances from the scattering matrix. The first dimension
    of the result has size 2 and indicates whether scattering is enabled with 0
    being without scattering and 1 being with scattering. The second dimension
    has size 2 and indicates which dimension is wanted with 0 being the x
    emittance and 1 being the y emittance. The third and last dimension has size
    output_dict['actual_analysis_points'] and indicates at which analysis point
    the data were taken.

    [scattering, coordinate, z]
    '''
    emittances_2d = np.empty((2, 2, covariance_matrix.shape[3]))
    emittances_2d[:, 0, :] = np.sqrt(np.linalg.det(np.transpose(
        covariance_matrix[:, 0:2, 0:2, :], (0, 3, 1, 2))))
    emittances_2d[:, 1, :] = np.sqrt(np.linalg.det(np.transpose(
        covariance_matrix[:, 2:4, 2:4, :], (0, 3, 1, 2))))
    return emittances_2d

def get4DEmittance(covariance_matrix):
    '''
    Computes the 4D transverse emittance from the scattering matrix. The first
    dimension of the result has size 2 and indicates whether scattering is
    enabled with 0 being without scattering and 1 being with scattering. The
    second dimension has size output_dict['actual_analysis_points'] and
    indicates at which analysis point the data were taken.

    [scattering, z]
    '''
    return np.sqrt(np.linalg.det(np.transpose(covariance_matrix, (0, 3, 1, 2))))

def getEnergy(output_dict, phase_space, gamma, bennett_radius, rho_ion_div_n0):
    '''
    Computes the energy (divided by m_e c^2) from the phase space. The first
    dimension of the result has size 2 and indicates whether scattering is
    enabled with 0 being without scattering and 1 being with scattering. The
    second dimension has size output_dict['actual_particles'] and indicates the
    number of the particle. The third and final dimension has size
    output_dict['actual_analysis_points'] indicates at which analysis point the
    data were taken.

    [scattering, particle, z]
    '''
    kinetic = 0.5 * (phase_space[:, 1] ** 2 + phase_space[:, 3] ** 2) / gamma
    r2 = phase_space[:, 0] ** 2 + phase_space[:, 2] ** 2
    potential = 0.25 * output_dict['ion_atomic_number'] * (output_dict['delta'] * r2 + (bennett_radius ** 2) * rho_ion_div_n0 * np.log(1 + r2 / (
        bennett_radius ** 2)))
    return kinetic + potential

class Simulation(object):

    def __init__(self, output_filename, read_phase_space=True):
        self._output_dict = getOutputDict(output_filename)
        self.z = np.linspace(0, self['plasma_length'],
            self['actual_analysis_points'])
        self.z_si = np.linspace(0, self['plasma_length_si'],
            self['actual_analysis_points'])
        self.gamma = self['gamma_initial'] + self['gamma_prime'] * self.z
        self.bennett_radius = self['bennett_radius_initial'] * ((self.gamma /
            self['gamma_initial']) ** (-1/4))
        self.bennett_radius_si = self['bennett_radius_initial_si'] * ((
            self.gamma / self['gamma_initial']) ** (-1/4))
        self.means, self.covariance_matrix = getStatistics(self._output_dict)
        if read_phase_space and self['output_phase_space']:
            self.phase_space = getPhaseSpace(self._output_dict, self.gamma)
            self.energy = getEnergy(self._output_dict, self.phase_space,
                self.gamma, self.bennett_radius, self['rho_ion_div_n0'])
        self.emittance_4d_si = get4DEmittance(self.covariance_matrix) * \
            (self['plasma_skin_depth_si'] ** 2)
        self.emittances_2d_si = get2DEmittances(self.covariance_matrix) * \
            self['plasma_skin_depth_si']

    def __getitem__(self, key):
        return self._output_dict[key]

    def plot2DEmittance(self, filename='results/2Demittance.png'):
        plt.plot(self.z_si, self.emittances_2d_si[0, 0, :] * 1e9,
            label='x (scattering off)', color='blue')
        plt.plot(self.z_si, self.emittances_2d_si[0, 1, :] * 1e9,
            label='y (scattering off)', color='green')
        plt.plot(self.z_si, self.emittances_2d_si[1, 0, :] * 1e9,
            label='x (scattering on)', color='red')
        plt.plot(self.z_si, self.emittances_2d_si[1, 1, :] * 1e9,
            label='y (scattering on)', color='orange')
        plt.xlim(0, self['plasma_length_si'])
        plt.xlabel(r'$z$ (m)')
        plt.ylabel(r'$\epsilon_N$ (nm)')
        plt.legend()
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
        plt.clf()

    def plotEnergy(self, filename='results/energy.png'):
        plt.plot(self.z_si, self.energy[0, 0])# * (self.gamma ** (1/4)))
        plt.xlim(0, self['plasma_length_si'])
        plt.xlabel(r'$z$ (m)')
        plt.ylabel(r'$\mathcal{H} / m_e c^2$')
        plt.savefig(filename)
        plt.clf()

    @staticmethod
    def plotScatteringTest(sims, labels, filename='results/scatteringtest.png'):
        # check assertions
        parameters = ('unperturbed_plasma_density_si', 'plasma_length_si',
            'beam_energy_initial_gev')
        for parameter in parameters:
            assert all(s[parameter] == sims[0][parameter] for s in sims)
        for s in sims:
            assert s['ion_atomic_number'] == 1
            assert s['acceleration_gradient_gev_per_m'] == 0

        for i in range(len(sims)):
            s = sims[i]
            label = labels[i]
            thx = np.std(np.arctan(s.phase_space[1, 1] / s.gamma), axis=0)
            plt.plot(s.z_si, thx, label=label)
            #plt.plot(s.z_si, thy, label='simulation {} (y)'.format(label))
        # get theoretical values
        z_values = np.linspace(0, sims[0]['plasma_length_si'], 10000)
        x = (1.6726219e-27 / 630.4) * sims[0]['unperturbed_plasma_density_si'] \
            * z_values
        theory = 26.61458558 * np.sqrt(x) * (1 + 0.038 * np.log(x)) / \
            sims[0]['gamma_initial']
        plt.plot(z_values, theory, label='Theory')
        plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 1))
        plt.xlabel(r'$z$ (m)')
        plt.ylabel(r'$\sigma_{\theta_x}$')
        plt.legend(loc='upper left')
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
        plt.clf()

    def createDistributionMovie(self):
        os.system('rm results/*.mp4')
        points = 1000
        bins = 100
        total_frames = 2 * self['actual_analysis_points']
        sigma = self.bennett_radius_si * np.sqrt(np.pi * 2.8179403227e-15 * \
            self['ion_atomic_number'] * self['rho_ion_si'] / (2 * self.gamma))
        for i, scattering in enumerate(('ns', 's')):
            for j in range(self['actual_analysis_points']):
                frame_number = j + i * self['actual_analysis_points']
                print('rendering frame ', frame_number, ' of ', total_frames)
                # get phase space coordinates
                x = self.phase_space[i, 0, :, j]
                vx = self.phase_space[i, 1, :, j] / self.gamma[j]
                y = self.phase_space[i, 2, :, j]
                vy = self.phase_space[i, 3, :, j] / self.gamma[j]
                r = np.sqrt(x ** 2 + y ** 2)
                th = np.arctan2(y, x)
                vr = (x * vx + y * vy) / r
                rvth = (x * vy - y * vx) / r
                # plot r distribution
                r_div_a_max = 5
                r_values = np.linspace(0, r_div_a_max * self.bennett_radius[0],
                    points)
                r = r[r < r_div_a_max * self.bennett_radius[0]]
                def unnormalizedPDF(r):
                    return r * np.exp(-r ** 2 / (2 * self['sigma_r'] ** 2)) / \
                        ((1 + ((r / self.bennett_radius[j]) ** 2)) ** 2)
                normalization_factor = scipy.integrate.quadrature(
                    unnormalizedPDF, 0, r_div_a_max * self.bennett_radius[j],
                    miniter=5)[0]
                r_pdf = unnormalizedPDF(r_values) / normalization_factor
                r_pdf_unmodified = 2 * self.bennett_radius[j] ** 2 * r_values \
                    / ((self.bennett_radius[j] ** 2 + r_values ** 2) ** 2)
                plt.title(r'z = {} m'.format(self.z_si[j]))
                plt.hist(r * self['plasma_skin_depth_si'], density=True,
                    bins=bins, range=(0, r_div_a_max *
                    self.bennett_radius_si[j]))
                plt.plot(r_values * self['plasma_skin_depth_si'], r_pdf *
                    self['plasma_angular_wavenumber_si'], label='modified')
                plt.plot(r_values * self['plasma_skin_depth_si'],
                    r_pdf_unmodified * self['plasma_angular_wavenumber_si'],
                    label='unmodified')
                plt.xlabel(r'$r$ [m]')
                plt.xlim(0, r_div_a_max * self.bennett_radius_si[j])
                #plt.ylim(0, 35000000 * 2)
                plt.legend()
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.savefig('frames/frame_r_{}_{}'.format(scattering, j))
                plt.clf()
                # plot theta distribution
                th_values = np.linspace(-np.pi, np.pi, points)
                plt.title(r'z = {} m'.format(self.z_si[j]))
                plt.hist(th, density=True, bins=bins, range=(-np.pi, np.pi))
                plt.plot(th_values, np.ones_like(th_values)/(2 * np.pi))
                plt.xlabel(r'$\theta$ [rad]')
                plt.xlim(-np.pi, np.pi)
                plt.ylim(0, 0.25)
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.savefig('frames/frame_th_{}_{}'.format(scattering, j))
                plt.clf()
                # plot r' distribution
                vr_values = np.linspace(-4 * sigma[0], 4 * sigma[0], points)
                plt.title(r'z = {} m'.format(self.z_si[j]))
                plt.hist(vr, density=True, bins=bins, range=(-4 * sigma[j],
                    4 * sigma[j]))
                plt.plot(vr_values, np.exp(-vr_values ** 2 / (2 * sigma[j] ** \
                    2)) / (np.sqrt(2 * np.pi) * sigma[j]))
                plt.xlabel(r'$\frac{dr}{dz}$ [unitless]')
                plt.xlim(-4 * sigma[j], 4 * sigma[j])
                #plt.ylim(0, 7000 * 2)
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.savefig('frames/frame_vr_{}_{}'.format(scattering, j))
                plt.clf()
                # plot r' distribution
                rvth_values = np.linspace(-4 * sigma[0], 4 * sigma[0], points)
                plt.title(r'z = {} m'.format(self.z_si[j]))
                plt.hist(rvth, density=True, bins=bins, range=(-4 * sigma[j], 4 *
                    sigma[j]))
                plt.plot(rvth_values, np.exp(-rvth_values ** 2 / (2 * sigma[j] ** \
                    2)) / (np.sqrt(2 * np.pi) * sigma[j]))
                plt.xlabel(r'$r \frac{d\theta}{dz}$ [unitless]')
                plt.xlim(-4 * sigma[j], 4 * sigma[j])
                #plt.ylim(0, 7000 * 2)
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.savefig('frames/frame_rvth_{}_{}'.format(scattering, j))
                plt.clf()
                # 4d
                plt.subplot(221)
                plt.hist(r * self['plasma_skin_depth_si'], density=True,
                    bins=bins, range=(0, r_div_a_max *
                    self.bennett_radius_si[j]), label='simulation')
                plt.plot(r_values * self['plasma_skin_depth_si'], r_pdf *
                    self['plasma_angular_wavenumber_si'], label='theory')
                plt.xlabel(r'$r$ [m]')
                plt.xlim(0, r_div_a_max * self.bennett_radius_si[j])
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.legend()
                plt.subplot(222)
                plt.hist(th, density=True, bins=bins, range=(-np.pi, np.pi), label='simulation')
                plt.plot(th_values, np.ones_like(th_values)/(2 * np.pi), label='theory')
                plt.xlabel(r'$\theta$ [rad]')
                plt.xlim(-np.pi, np.pi)
                plt.ylim(0, 0.3)
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.legend()
                plt.subplot(223)
                # plot r' distribution
                vr_values = np.linspace(-4 * sigma[0], 4 * sigma[0], points)
                plt.hist(vr, density=True, bins=bins, range=(-4 * sigma[j],
                    4 * sigma[j]), label='simulation')
                plt.plot(vr_values, np.exp(-vr_values ** 2 / (2 * sigma[j] ** \
                    2)) / (np.sqrt(2 * np.pi) * sigma[j]), label='theory')
                plt.xlabel(r'$\frac{dr}{dz}$ [unitless]')
                plt.xlim(-4 * sigma[j], 4 * sigma[j])
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.legend()
                plt.subplot(224)
                # plot r' distribution
                rvth_values = np.linspace(-4 * sigma[0], 4 * sigma[0], points)
                plt.hist(rvth, density=True, bins=bins, range=(-4 * sigma[j], 4 *
                    sigma[j]), label='simulation')
                plt.plot(rvth_values, np.exp(-rvth_values ** 2 / (2 * sigma[j] ** \
                    2)) / (np.sqrt(2 * np.pi) * sigma[j]), label='theory')
                plt.xlabel(r'$r \frac{d\theta}{dz}$ [unitless]')
                plt.xlim(-4 * sigma[j], 4 * sigma[j])
                plt.gca().get_xaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits(
                    (0, 1))
                plt.legend()
                plt.tight_layout()
                plt.savefig('frames/combined_{}_{}'.format(scattering, j), dpi=500)
                plt.clf()
            os.system('ffmpeg -i \'frames/frame_r_{}_%d.png\' -vcodec libx264 -'
                'vf scale=640:-2,format=yuv420p results/movie_r_{}.mp4'.format(
                scattering, scattering))
            os.system('ffmpeg -i \'frames/frame_th_{}_%d.png\' -vcodec libx264 '
                '-vf scale=640:-2,format=yuv420p results/movie_th_{}.mp4'
                .format(scattering, scattering))
            os.system('ffmpeg -i \'frames/frame_vr_{}_%d.png\' -vcodec libx264 '
                '-vf scale=640:-2,format=yuv420p results/movie_vr_{}.mp4'
                .format(scattering, scattering))
            os.system('ffmpeg -i \'frames/frame_rvth_{}_%d.png\' -vcodec libx26'
                '4 -vf scale=640:-2,format=yuv420p results/movie_rvth_{}.mp4'
                .format(scattering, scattering))

import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.integrate
import scipy.signal

def plotBennettParameters(output_dict, z, gamma, rho_ion, bennett_radius):
    plt.plot(z, gamma)
    plt.title('Lorentz Factor')
    plt.ylim((output_dict['gamma_initial'], output_dict['gamma_final']))
    plt.xlabel(r'$z$ [unitless]')
    plt.ylabel(r'$\gamma$ [unitless]')
    plt.savefig('results/lorentz_factor.png')
    plt.cla()
    plt.plot(z, 1 + rho_ion)
    plt.title('Peak Density')
    plt.ylim((output_dict['rho_ion_initial'], output_dict['rho_ion_final']))
    plt.xlabel(r'$z$ [unitless]')
    plt.ylabel(r'peak ion density [units of the plasma density]')
    plt.savefig('results/peak_density.png')
    plt.cla()
    plt.plot(z, bennett_radius)
    plt.title('Bennett Radius')
    plt.ylim((output_dict['bennett_radius_final'], output_dict['bennett_radius_initial']))
    plt.xlabel(r'$z$ [unitless]')
    plt.ylabel(r'$a$ [unitless]')
    plt.savefig('results/bennett_radius.png')
    plt.cla()

def plotTraj(output_dict, z, phase_space):
    plt.title('x')
    plt.plot(z, phase_space[0, 0, 1, :], label='x')
    plt.plot(z, phase_space[0, 2, 1, :], label='y')
    plt.legend()
    plt.savefig('results/x.png')
    plt.cla()
    plt.title('p')
    plt.plot(z, phase_space[0, 1, 1, :], label='px')
    plt.plot(z, phase_space[0, 3, 1, :], label='py')
    plt.legend()
    plt.savefig('results/p.png')
    plt.cla()

def normal_pdf(x, mu, sigma):
    return 1.0 / (sigma * (2.0 * np.pi)**(1/2)) * np.exp(-(x - mu)**2 / (2.0 * (sigma**2)))

def createHistogramMovie(output_dict, z, phase_space, bennett_radius, rho_ion, gamma):
    r_div_a_max = 10
    sigma = output_dict['plasma_skin_depth_si'] * bennett_radius * np.sqrt(np.pi * 2.8179403227e-15 * output_dict['ion_atomic_number'] * rho_ion * output_dict['unperturbed_plasma_density_si'] / (2 * gamma))
    for i, scattering in enumerate(('ns',)):
        for j in range(output_dict['actual_analysis_points']):
            print('frame {} of {} ({})'.format(j, output_dict['actual_analysis_points'], scattering))
            x = phase_space[i, 0, :, j]
            vx = phase_space[i, 1, :, j] / gamma[j]
            y = phase_space[i, 2, :, j]
            vy = phase_space[i, 3, :, j] / gamma[j]
            r = np.sqrt(x ** 2 + y ** 2)
            th = np.arctan2(y, x)
            vr = (x * vx + y * vy) / r
            rvth = (x * vy - y * vx) / r
            rs = np.linspace(0, r_div_a_max * bennett_radius[i], 1000)
            ths = np.linspace(-np.pi, np.pi, 1000)
            vrs = np.linspace(-4 * sigma[j], 4 * sigma[j], 1000)
            rvths = np.linspace(-4 * sigma[j], 4 * sigma[j], 1000)
            def unnormalizedpdf(r2):
                return r2 * np.exp(-r2 ** 2 / (2 * output_dict['sigma_r_initial'] ** 2)) / \
                    ((1 + ((r2 / bennett_radius[j]) ** 2)) ** 2)
            normalization_factor = scipy.integrate.quadrature(unnormalizedpdf, 0, r_div_a_max * bennett_radius[j], miniter=5)[0]
            theoretical_particle_density = unnormalizedpdf(rs) / normalization_factor
            theoretical_particle_density2 = 2 * bennett_radius[j] * bennett_radius[j] * rs / ((bennett_radius[j] * bennett_radius[j] + rs * rs) ** 2)
            plt.title(r'z = {}'.format(z[j]))
            plt.hist(r, density=True, bins=100, range=(0, r_div_a_max * bennett_radius[j]))
            plt.plot(rs, theoretical_particle_density, label='modified')
            plt.plot(rs, theoretical_particle_density2, label='unmodified')
            #plt.xlim(0.01, 0.04)
            #plt.ylim(0, 30)
            plt.legend()
            plt.savefig('frames/movie_r_{}'.format(j))
            plt.cla()

            plt.title(r'z = {}'.format(z[j]))
            plt.hist(th, density=True, bins=100, range=(-np.pi, np.pi))
            plt.plot(ths, np.ones_like(ths)/(2 * np.pi))
            plt.savefig('frames/movie_th_{}'.format(j))
            plt.cla()

            plt.title(r'z = {}'.format(z[j]))
            plt.hist(vr, density=True, bins=100, range=(-4 * sigma[j], 4 * sigma[j]))
            plt.plot(vrs, normal_pdf(vrs, 0, sigma[j]))
            plt.savefig('frames/movie_vr_{}'.format(j))
            plt.cla()

            plt.title(r'z = {}'.format(z[j]))
            plt.hist(rvth, density=True, bins=100, range=(-4 * sigma[j], 4 * sigma[j]))
            plt.plot(rvths, normal_pdf(rvths, 0, sigma[j]))
            plt.savefig('frames/movie_rvth_{}'.format(j))
            plt.cla()

    os.system('ffmpeg -i \'frames/movie_r_%d.png\' -vcodec libx264 -vf scale=640:-2,format=yuv420p results/movie_r.mp4')
    os.system('ffmpeg -i \'frames/movie_th_%d.png\' -vcodec libx264 -vf scale=640:-2,format=yuv420p results/movie_th.mp4')
    os.system('ffmpeg -i \'frames/movie_vr_%d.png\' -vcodec libx264 -vf scale=640:-2,format=yuv420p results/movie_vr.mp4')
    os.system('ffmpeg -i \'frames/movie_rvth_%d.png\' -vcodec libx264 -vf scale=640:-2,format=yuv420p results/movie_rvth.mp4')

def plotEnergy(output_dict, z, energy):
    plt.title('energy')
    plt.plot(z, energy[0, 0, :])
    plt.savefig('results/energy.png')
    plt.cla()

def plot4DEmittanceGrowth(output_dict, z, emittance_4d, smoothing_on, window_size=101, order=3):
    ns_emit = np.copy(emittance_4d[0, :])
    s_emit = np.copy(emittance_4d[1, :])
    z2 = np.copy(z)
    #avg_emit = np.mean(ns_emit)
    #ns_emit /= avg_emit
    #s_emit /= avg_emit
    if smoothing_on:
        ns_emit = scipy.signal.savgol_filter(ns_emit, window_size, order)
        s_emit = scipy.signal.savgol_filter(s_emit, window_size, order)
        ns_emit = ns_emit[window_size:-window_size]
        s_emit = s_emit[window_size:-window_size]
        z2 = z2[window_size:-window_size]
        plt.title('4D RMS Emittance Growth (smoothed, window_size={}, order={})'.format(window_size, order))
    else:
        plt.title('4D RMS Emittance Growth')
    plt.plot(z2, output_dict['plasma_skin_depth_si'] * output_dict['plasma_skin_depth_si'] * ns_emit, label='no scattering', color='blue')
    plt.plot(z2, output_dict['plasma_skin_depth_si'] * output_dict['plasma_skin_depth_si'] * s_emit, label='scattering', color='red')
    plt.xlim(0, output_dict['plasma_length'])
    plt.xlabel(r'$k_p z$ [unitless]')
    plt.ylabel(r'$\epsilon_{4D}(z)$ [m^2]')
    plt.legend()
    if smoothing_on:
        plt.savefig('results/emit_4d_window_{}_order_{}.png'.format(window_size, order))
    else:
        plt.savefig('results/emit_4d.png')
    plt.cla()

def plot2DEmittanceGrowth(output_dict, z, emittances_2d, smoothing_on, window_size=101, order=3):
    ns_x_emit = np.copy(emittances_2d[0, 0, :])
    ns_y_emit = np.copy(emittances_2d[0, 1, :])
    s_x_emit = np.copy(emittances_2d[1, 0, :])
    s_y_emit = np.copy(emittances_2d[1, 1, :])
    z2 = np.copy(z)
    #avg_emit = np.mean(np.concatenate((ns_x_emit, ns_y_emit)))
    #ns_x_emit /= avg_emit
    #ns_y_emit /= avg_emit
    #s_x_emit /= avg_emit
    #s_y_emit /= avg_emit
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
    plt.plot(z, output_dict['plasma_skin_depth_si'] * ns_x_emit, label='x (no scattering)', color='blue')
    plt.plot(z, output_dict['plasma_skin_depth_si'] * ns_y_emit, label='y (no scattering)', color='green')
    plt.plot(z, output_dict['plasma_skin_depth_si'] * s_x_emit, label='x (scattering)', color='red')
    plt.plot(z, output_dict['plasma_skin_depth_si'] * s_y_emit, label='y (scattering)', color='orange')
    plt.xlim(0, output_dict['plasma_length'])
    plt.xlabel(r'$k_p z$ [unitless]')
    plt.ylabel(r'$\epsilon_{2D}(z)$ [m]')
    plt.legend()
    if smoothing_on:
        plt.savefig('results/emit_2d_window_{}_order_{}.png'.format(window_size, order))
    else:
        plt.savefig('results/emit_2d.png')
    plt.cla()

def generate2DEmittanceGrowthRates(output_dict, z, emittances_2d):
    s_x_emit = np.copy(emittances_2d[1, 0, :])
    s_y_emit = np.copy(emittances_2d[1, 1, :])
    #s_x_emit /= s_x_emit[0]
    #s_y_emit /= s_y_emit[0]
    slope_x, intercept_x, r_value_x, p_value_x, std_err_x = scipy.stats.linregress(z, output_dict['plasma_skin_depth_si'] * s_x_emit);
    slope_y, intercept_y, r_value_y, p_value_y, std_err_y = scipy.stats.linregress(z, output_dict['plasma_skin_depth_si'] * s_y_emit);
    with open('results/growth_rate.txt', 'w') as f:
        f.write('g_x = {} +/- {}\n'.format(slope_x, std_err_x))
        f.write('g_y = {} +/- {}\n'.format(slope_y, std_err_y))

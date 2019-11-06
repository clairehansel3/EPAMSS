import run
import sys
import simulation
import numpy as np

exponents1 = np.arange(-10, -6)
exponents2 = np.arange(1, 7)

def submit_jobs():
    for run_number, exponent in enumerate(exponents1):
        runname = 'rsigma_{}'.format(exponent)
        run.run(
            run_name=runname,
            parameters={
                'rho_ion_si': 1e26,
                'plasma_length_si': 100,
                'acceleration_gradient_gev_per_m': 0,
                'unperturbed_plasma_density_si': 1e26,
                'particles_target': 1000,
                'cross_section_radius_si': 10. ** exponent,
                'minimum_steps_per_betatron_period': 10 ** 4
            },
            hoffman2=False,
            compute_processes=7
        )
    for run_number, exponent in enumerate(exponents2):
        runname = 'steps_{}'.format(exponent)
        run.run(
            run_name=runname,
            parameters={
                'rho_ion_si': 1e26,
                'plasma_length_si': 100,
                'acceleration_gradient_gev_per_m': 0,
                'unperturbed_plasma_density_si': 1e26,
                'particles_target': 1000,
                'cross_section_radius_si': 10. ** -8,
                'minimum_steps_per_betatron_period': 10 ** exponent
            },
            hoffman2=False,
            compute_processes=7
        )

def analyze():
    sims1 = []
    labels1 = []
    for run_number, exponent in enumerate(exponents1):
        runname = 'rsigma_{}'.format(exponent)
        sims1.append(simulation.Simulation('data/{}_output'.format(runname)))
        labels1.append(r'$r_{\sigma} = 10^{' + str(exponent) + '}$')
    simulation.Simulation.plotScatteringTest(sims1, labels1, filename='results/rsigmalog.png')
    sims2 = []
    labels2 = []
    for run_number, exponent in enumerate(exponents2):
        runname = 'steps_{}'.format(exponent)
        sims2.append(simulation.Simulation('data/{}_output'.format(runname)))
        labels2.append('$10^{' + str(exponent) + '}$ steps')
    simulation.Simulation.plotScatteringTest(sims2, labels2, filename='results/stepslog.png')

if sys.argv[1] == 'run':
    submit_jobs()
elif sys.argv[1] == 'analyze':
    analyze()

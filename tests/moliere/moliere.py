import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 14})

varthetas = []
pdfs = []

with open('pdf') as f:
    for line in f:
        vartheta, *pdf = list(float(i) for i in line.split())
        varthetas.append(vartheta)
        pdfs.append(pdf)

varthetas = np.array(varthetas)
pdfs = np.array(pdfs).T
plt.plot(varthetas, pdfs[0], label='order 0')
plt.plot(varthetas, pdfs[1], label='order 1')
plt.plot(varthetas, pdfs[2], label='order 2')
plt.plot(varthetas, pdfs[3], label='order 3')
plt.xlim(0, 5)
plt.ylim(0, 1)
plt.legend()
plt.xlabel(r"$\vartheta$")
plt.ylabel(r"$\mathrm{PDF}(\vartheta)$")
plt.savefig('pdf.png', dpi=500)
plt.cla()

samples = []

with open('samples') as f:
    for line in f:
        samples.append(list(float(i) for i in line.split()))

samples = np.array(samples).T

for order in (0, 1, 2, 3):
    plt.plot(varthetas, pdfs[order], label='pdf', color='C3')
    plt.hist(samples[order], color='C0', density=True, bins=50, range=(0, 5))
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    plt.xlabel(r"$\vartheta$")
    plt.ylabel('probability density')
    plt.savefig('samples_order_{}.png'.format(order), dpi=500)
    plt.cla()

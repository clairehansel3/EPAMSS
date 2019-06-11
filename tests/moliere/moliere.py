import matplotlib.pyplot as plt
import numpy as np

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
plt.ylabel(r"$PDF(\vartheta)$")
plt.savefig('pdf.png')
plt.cla()

samples = []

with open('samples') as f:
    for line in f:
        samples.append(list(float(i) for i in line.split()))

samples = np.array(samples).T

for order in (0, 1, 2, 3):
    plt.title('order {}'.format(order))
    plt.plot(varthetas, pdfs[order], label='pdf')
    plt.hist(samples[order], density=True, bins=100, range=(0, 5))
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    plt.legend()
    plt.xlabel(r"$\vartheta$")
    plt.savefig('samples_order_{}.png'.format(order))
    plt.cla()

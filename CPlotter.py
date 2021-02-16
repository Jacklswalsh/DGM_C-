import matplotlib.pyplot as plt
import csv
import numpy as np

x = []
y = []
y2 = []
y3 = []




with open('/Users/Jacklswalsh/Desktop/Desktop/C++/PhD/DGM_c-/output.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
        y2.append(float(row[2]))
        y3.append(float(row[3]))

    xReal = np.arange(-1,1,2/1000)
    sigma = 0.2
    fReal0 = np.exp(-np.log(2) * np.power((xReal + 0.5), 2) / sigma ** 2)
    fReal1 = np.exp(-np.log(2) * np.power((xReal + 0), 2) / sigma ** 2)
    fReal2 = np.exp(-np.log(2) * np.power((xReal - 0.5), 2) / sigma ** 2)

plt.subplots()
plt.scatter(x,y, label='T = 0')
plt.scatter(x,y2, label='T = 0.5')
plt.scatter(x,y3, label='T = 1')


plt.plot(xReal, fReal0, color='b', label='T (Real) = {}'.format(0))
plt.plot(xReal, fReal1, color='r', label='T (Real) = {}'.format(0.5))
plt.plot(xReal, fReal2, color='k', label='T (Real) = {}'.format(1))

plt.xlabel('x')
plt.ylabel('y')
plt.title('Discontinuous Galerkin Method \nWave Equation Benchmark Soln \n (C++ Implementation)')
plt.legend()
plt.show()
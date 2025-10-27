from cell import Cell
import matplotlib.pyplot as plt
from utils import make_data

dt, time, current, voltage, soc, ocv, T = make_data(short=True)

cell1 = Cell()
cell1.soc = soc
cell1.T = T
mean, sigma = cell1.solve(current=current, time=time)

# Plot the solution
fig, ax = plt.subplots()
ax.plot(time, voltage, label='Measured', c='b')
ax.fill_between(time, mean - 2 * sigma, mean + 2 * sigma, alpha=0.3, color='r')
ax.plot(time, mean, c='r', label='Simulated')
fig.legend()
ax.set_xlabel('t/[s]')
plt.show()
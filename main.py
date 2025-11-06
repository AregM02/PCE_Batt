from cell import Cell
import matplotlib.pyplot as plt
from utils import make_data, load_validation, add_trace

time, current, voltage, soc, T, C_nom = load_validation()

cell = Cell(initial_soc=soc[0], capacity=C_nom, capacity_unc=C_nom*0.01)
cell.solve(current=current, time=time, temperature=T)

# Plot the solution
fig, ax = plt.subplots()
ax.plot(time, voltage, label='Measurement', c='b')
add_trace(ax, time, cell.voltage.mean, cell.voltage.std, name='voltage', c='r')
fig.legend()
ax.set_xlabel('t/[s]')
plt.show()
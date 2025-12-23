from cell import Cell
import matplotlib.pyplot as plt
from utils import load_validation, add_trace, RMSE

time, current, measurement, soc, T, C_nom = load_validation()

cell = Cell(initial_soc=soc[0], capacity=C_nom, capacity_unc=C_nom*0.002)
cell.solve(current=current, time=time, temperature=T)

# # Load cached PCE simulation at 25deg (saves time)
# cell.voltage.mean, cell.voltage.std = unify_measurements(cell.voltage.mean, cell.voltage.std, *load_cached('vpce_25'))

# -----PLOTS-----
fig, ax = plt.subplots(figsize = (24, 16))
add_trace(ax, time, measurement[0], measurement[1], c='b', legend=['Real Mean', '95% CI'])
add_trace(ax, time, cell.voltage.mean, cell.voltage.std, c='r', legend=['Model Mean', '95% CI'])
ax.legend(loc = 'upper right')
ax.set_xlabel('Time / [s]')
ax.set_ylabel('Voltage / [V]')
ax.grid(True, alpha = 0.3)
# plt.savefig("profile.pdf", dpi=600)
plt.show()

# Detailed comparison
print(f"RMSE: {1000*RMSE(cell.voltage.mean, measurement[0]):.2f} mV")
print(f"STD RMSE: {1000*RMSE(cell.voltage.std, measurement[1]):.2f} mV")
fig, ax = plt.subplots()
ax.plot(measurement[1], c='b', label='Measured Std.')
ax.plot(cell.voltage.std, c='r', label='Simulated Std.')
ax.legend()
ax.grid(True, alpha=0.3)
# plt.savefig("stds.pdf", dpi=600)
plt.show()
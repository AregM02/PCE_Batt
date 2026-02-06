import matplotlib.pyplot as plt
from src.utils import load_validation, add_trace, RMSE, plot_pack_voltage, plot_current_imbalance, compare_pack_voltage,compare_current_imbalance
from src.core import PCPack, Cell, MCPack

time, current, measurement, soc, temperature, C_nom, C_nom_unc = load_validation(short=True)
cell_params = dict(initial_soc=soc[0], capacity=C_nom, capacity_unc=C_nom*C_nom_unc)

pcpack = PCPack("1-p(3)-p(10)", cell_params=cell_params)
pcpack.solve(current=current, time=time, temperature=temperature)

mcpack = MCPack("2-p(3)-p(10)", cell_params=cell_params)
mcpack.solve(current=current, time=time, temperature=temperature)


compare_pack_voltage(pack1=pcpack, pack2=mcpack, labels=["pcpack", "mcpack"])
compare_current_imbalance(pack1=pcpack, pack2=mcpack, labels=["pcpack", "mcpack"])

# # -----PLOTS-----
# fig, ax = plt.subplots(figsize = (24, 16))
# add_trace(ax, time, measurement[0], measurement[1], c='b', legend=['Real Mean', '95% CI'])
# add_trace(ax, time, mcpack.voltage.mean, mcpack.voltage.std, c='r', legend=['Model Mean', '95% CI'])
# ax.legend(loc = 'upper right')
# ax.set_xlabel('Time / [s]')
# ax.set_ylabel('Voltage / [V]')
# ax.grid(True, alpha = 0.3)
# # plt.savefig("profile.pdf", dpi=600)
# plt.show()

# # Detailed comparison
# print(f"RMSE: {1000*RMSE(pack.voltage.mean, measurement[0]):.2f} mV")
# print(f"STD RMSE: {1000*RMSE(pack.voltage.std, measurement[1]):.2f} mV")
# fig, ax = plt.subplots()
# ax.plot(measurement[1], c='b', label='Measured Std.')
# ax.plot(pack.voltage.std, c='r', label='Simulated Std.')
# ax.plot()
# # ax.plot(cell.voltage.std/measurement[1], c='black', label='Simulated Std./Measured Std.', linewidth = 0.7)
# ax.grid(True, alpha=0.3)
# # plt.savefig("stds.pdf", dpi=600)
# plt.show()
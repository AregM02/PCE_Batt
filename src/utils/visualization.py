import matplotlib.pyplot as plt
import numpy as np
from src.utils import create_logger

# Initialize logger for the module
logger = create_logger("Visualization")


def add_trace(ax, x, mean, sigma, c, legend=[None, None], lw=0.7):
    if not legend[0]:
        ax.plot(x, mean, c=c, linewidth=lw)
    else:
        ax.plot(x, mean, c=c, label=legend[0], linewidth=lw)
    if not legend[1]: 
        ax.fill_between(x, mean - 2 * sigma, mean + 2 * sigma,
                        alpha=0.5, color=c, linewidth=lw)
    else:
        ax.fill_between(x, mean - 2 * sigma, mean + 2 * sigma,
                        alpha=0.5, color=c, label=legend[1], linewidth=lw)


def plot_pack_voltage(pack, num_samples: int = 5):
    """
    Plots the pack voltage results directly from a Pack object using pack.time.
    """
    if pack.voltage.mean is None:
        logger.error("Pack has not been solved yet. Call pack.solve() before plotting.")
        return

    time = pack.time 
    voltage_ns = pack.voltage
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    add_trace(ax, time, voltage_ns.mean, voltage_ns.std, c='tab:blue', 
              legend=["Mean Pack Voltage", r"95% Confidence Interval (2$\sigma$)"])
    
    if voltage_ns.samples is not None:
        indices = np.random.choice(pack.n_virtual_packs, min(num_samples, pack.n_virtual_packs), replace=False)
        for i in indices:
            ax.plot(time, voltage_ns.samples[i, :], color='tab:blue', alpha=0.15, lw=0.5)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Terminal Voltage [V]")
    ax.set_title(f"Pack Voltage Response: {pack.topology_str}")
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_current_imbalance(pack):
    """
    Visualizes current distribution using pack.time.
    """
    if not hasattr(pack, 'time') or pack.time is None:
        logger.error("Pack object is missing the 'time' attribute.")
        return

    if not pack.currents.per_block_std:
        logger.warning("No parallel blocks found or pack not solved.")
        return

    time = pack.time
    currents_ns = pack.currents
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    num_blocks = len(currents_ns.per_block_std)
    colors = plt.cm.plasma(np.linspace(0, 0.7, num_blocks))
    
    for i, std_series in enumerate(currents_ns.per_block_std):
        ax.plot(time, std_series, color=colors[i], lw=1.5, 
                label=fr"Block {i} Distribution ($\sigma_I$)")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Current Std Dev [A]")
    ax.set_title(f"Current Distribution: {pack.topology_str}")
    ax.legend()
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.show()
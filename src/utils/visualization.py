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


# TODO: this can be moved into the Pack class for more intutitive calls
def plot_pack_voltage(pack, num_samples: int = 5):
    """
    Plots the pack voltage results directly from a Pack object using pack.time.
    """
    if pack.voltage.mean is None:
        logger.error("Pack has not been simulated yet. Call pack.solve() before plotting.")
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


def compare_pack_voltage(pack1, pack2, labels=["Pack 1", "Pack 2"], num_samples: int = 5):
    """
    Plots two pack voltage results side-by-side for direct comparison.
    """
    # Check if both packs are simulated
    for p, label in zip([pack1, pack2], labels):
        if p.voltage.mean is None:
            print(f"Error: {label} has not been simulated yet.")
            return

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    
    for i, (pack, ax) in enumerate(zip([pack1, pack2], axes)):
        time = pack.time 
        voltage_ns = pack.voltage
        
        # Plot mean and confidence interval using the helper
        add_trace(ax, time, voltage_ns.mean, voltage_ns.std, c='tab:blue', 
                  legend=[f"Mean {labels[i]}", r"95% CI (2$\sigma$)"])
        
        # Plot random samples for the "identity" spread
        if voltage_ns.samples is not None:
            indices = np.random.choice(pack.n_virtual_packs, 
                                       min(num_samples, pack.n_virtual_packs), 
                                       replace=False)
            for idx in indices:
                ax.plot(time, voltage_ns.samples[idx, :], color='tab:blue', alpha=0.15, lw=0.5)

        ax.set_xlabel("Time [s]")
        if i == 0:
            ax.set_ylabel("Terminal Voltage [V]")
        
        ax.set_title(f"{labels[i]} ({pack.topology_str})")
        ax.legend(loc='lower left')
        ax.grid(True, alpha=0.3)

    plt.suptitle("Pack Voltage Response Comparison", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def compare_current_imbalance(pack1, pack2, labels=["Pack 1", "Pack 2"]):
    """
    Visualizes current distribution (standard deviation) side-by-side.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    
    for i, (pack, ax) in enumerate(zip([pack1, pack2], axes)):
        if not hasattr(pack, 'time') or pack.time is None or not pack.currents.per_block_std:
            ax.text(0.5, 0.5, f"No Current Data for {labels[i]}", 
                    ha='center', va='center', transform=ax.transAxes)
            continue

        time = pack.time
        currents_ns = pack.currents
        num_blocks = len(currents_ns.per_block_std)
        colors = plt.cm.plasma(np.linspace(0, 0.7, num_blocks))
        
        for b_idx, std_series in enumerate(currents_ns.per_block_std):
            ax.plot(time, std_series, color=colors[b_idx], lw=1.5, 
                    label=fr"Block {b_idx} ($\sigma_I$)")

        ax.set_xlabel("Time [s]")
        if i == 0:
            ax.set_ylabel("Current Std Dev [A]")
        
        ax.set_title(f"{labels[i]} Current Spread")
        ax.legend(loc='best')
        ax.grid(True, alpha=0.2)

    plt.suptitle("Current Imbalance Comparison (Parallel Blocks)", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
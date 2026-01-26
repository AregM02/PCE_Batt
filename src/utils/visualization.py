import matplotlib.pyplot as plt


def add_trace(ax, x, mean, sigma, c, legend = [None, None], lw = 0.7):
    if not legend[0]:
        ax.plot(x, mean, c=c, linewidth = lw)
    else:
        ax.plot(x, mean, c=c, label = legend[0], linewidth = lw)
    if not legend[1]: 
        ax.fill_between(x, mean - 2 * sigma, mean + 2 * sigma,
                        alpha=0.5, color=c, linewidth = lw)
    else:
        ax.fill_between(x, mean - 2 * sigma, mean + 2 * sigma,
                        alpha=0.5, color=c, label = legend[1], linewidth = lw)


def plot_dist(pce, joint, variable_num, samples_n=10000000):
        samples_pce = joint.sample(samples_n)
        x = pce[-1](**{f'q{i}':samples_pce[i] for i in range(variable_num)})
        plt.hist(x, bins=50, edgecolor='black')
        plt.show()
        quit()
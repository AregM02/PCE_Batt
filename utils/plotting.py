import matplotlib.pyplot as plt

def add_trace(ax, x, mean, sigma, name, c):
    ax.plot(x, mean, c=c, label=f'{name.capitalize()}: mean')
    ax.fill_between(x, mean - 2 * sigma, mean + 2 * sigma, alpha=0.3, color=c, label=f'{name.capitalize()}: 95% conf. interval')


def plot_dist(pce, joint, variable_num, samples_n=10000000):
        samples_pce = joint.sample(samples_n)
        x = pce[-1](**{f'q{i}':samples_pce[i] for i in range(variable_num)})
        plt.hist(x, bins=50, edgecolor='black')
        plt.show()
        quit()
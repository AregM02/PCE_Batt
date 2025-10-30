import matplotlib.pyplot as plt

def add_trace(ax, x, mean, sigma, name, c):
    ax.plot(x, mean, c=c, label=f'{name.capitalize()}: mean')
    ax.fill_between(x, mean - 2 * sigma, mean + 2 * sigma, alpha=0.3, color=c, label=f'{name.capitalize()}: 95% conf. interval')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

def get_fig_axes():
    fig = plt.figure(figsize = (10, 9))
    spec = gs.GridSpec(2, 4, figure = fig)
    ax1 = fig.add_subplot(spec[0, :2])
    ax2 = fig.add_subplot(spec[0, 2:])
    ax3 = fig.add_subplot(spec[1, 1:3])
    return (fig, spec, ax1, ax2, ax3)
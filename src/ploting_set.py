import matplotlib as mpl
from matplotlib import pyplot as plt

def acs_plot_style(use_style: bool = True) -> None:
    """
    Set the plot style to ACS.
    Parameters
    ----------
    use_style : bool
        If True, set the plot style to ACS. If False, reset to default.
    """
    if use_style:
        plt.style.use('acs')
    else:
        # Font settings
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = ['Arial']
        mpl.rcParams['font.size'] = 8

        # Figure settings
        mpl.rcParams['figure.figsize'] = (3.5, 2.5)
        mpl.rcParams['figure.dpi'] = 300
        mpl.rcParams['figure.facecolor'] = 'white'

        # Axes settings
        mpl.rcParams['axes.labelsize'] = 8
        mpl.rcParams['axes.titlesize'] = 9
        mpl.rcParams['axes.linewidth'] = 0.8
        mpl.rcParams['axes.edgecolor'] = 'black'
        mpl.rcParams['axes.spines.top'] = True
        mpl.rcParams['axes.spines.right'] = True
        mpl.rcParams['axes.grid'] = False

        # Ticks settings
        mpl.rcParams['xtick.direction'] = 'out'
        mpl.rcParams['ytick.direction'] = 'out'
        mpl.rcParams['xtick.major.size'] = 3
        mpl.rcParams['ytick.major.size'] = 3
        mpl.rcParams['xtick.major.width'] = 0.8
        mpl.rcParams['ytick.major.width'] = 0.8
        mpl.rcParams['xtick.labelsize'] = 9
        mpl.rcParams['ytick.labelsize'] = 9

        # Lines and markers settings
        mpl.rcParams['lines.linewidth'] = 1.2
        mpl.rcParams['lines.markersize'] = 4

        # Legend settings
        mpl.rcParams['legend.fontsize'] = 7
        mpl.rcParams['legend.frameon'] = False

        # Savefig settings
        mpl.rcParams['savefig.dpi'] = 600
        mpl.rcParams['savefig.transparent'] = False
        mpl.rcParams['savefig.bbox'] = 'tight'
        mpl.rcParams['savefig.pad_inches'] = 0.05

##################################################################################################

#http://bkanuka.com/articles/native-latex-plots/


def figsize(scale):
    from numpy import sqrt
    
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean*1.2              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size
    
def first_run():
    import matplotlib as mpl
    mpl.use('pgf')

    pgf_with_latex = {                      # setup matplotlib to use latex for output
        "pgf.texsystem": "xelatex",        # change this if using xetex or lautex
        "text.usetex": True,                # use LaTeX to write all text
        "font.family": "serif",
        "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
        "font.sans-serif": [],
        "font.monospace": [],
        "axes.labelsize": 10,               # LaTeX default is 10pt font.
        "font.size": 10,
        "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
        "pgf.preamble": [
            r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
            r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
            ]
        }

    mpl.rcParams.update(pgf_with_latex)
    
    pgf_with_rc_fonts = {"font.serif": []}# use latex default serif font

    mpl.rcParams.update(pgf_with_rc_fonts)

# I make my own newfig and savefig functions
def newfig(width):
    from matplotlib.pyplot import clf, figure
    clf() #limpa a figura atual
    fig = figure(figsize=figsize(width))
    return fig

def savefig(filename):
    from matplotlib.pyplot import savefig
    savefig('{}.pgf'.format(filename))
    savefig('{}.pdf'.format(filename))

if __name__ == '__savefig__':
    savefig()
elif __name__ == '__newfig__':
    newfig()
elif __name__ == '__first_run__':
    first_run()
elif __name__ == '__figsize__':
    figsize()
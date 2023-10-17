from matplotlib import pyplot as plt
from matplotlib import colors


# ######################
# ##### Color Maps #####
# ######################

def initialize():
    fig, axes = plt.subplots(figsize=(8,8))
    plt.close()



# Initialize a new plot
def newplot(width = 8, height = 8, fontsize = 20, style = 'serif', fontset = "cm", auto_layout = True):


    fig, axes = plt.subplots(figsize=(width,height))

    plt.rcParams['figure.figsize'] = (width,height)
    plt.rcParams['font.family'] = style
    plt.rcParams['mathtext.fontset']= fontset
    plt.rcParams['figure.autolayout'] = auto_layout
    plt.rcParams['font.size'] = str(fontsize)

    return fig, axes


# Red
red_cdict = {'red':   ((0.0,  0.0, 0.0),
                       (1.0,  1.0, 1.0)),

             'green': ((0.0,  0.0, 0.0),
                       (1.0,  0.0, 0.0)),

             'blue':  ((0.0,  0.0, 0.0),
                       (1.0,  0.0, 0.0)), }

red_cmap = colors.LinearSegmentedColormap('custom', red_cdict)


# Yellow
yellow_cdict = {'red':   ((0.0,  0.0, 0.0),
                          (1.0,  1.0, 1.0)),

                'green': ((0.0,  0.0, 0.0),
                          (1.0,  1.0, 1.0)),

                'blue':  ((0.0,  0.0, 0.0),
                          (1.0,  0.0, 0.0)), }

yellow_cmap = colors.LinearSegmentedColormap('custom', yellow_cdict)


# Green
green_cdict = {'red':   ((0.0,  0.0, 0.0),
                         (1.0,  0.0, 0.0)),

               'green': ((0.0,  0.0, 0.0),
                         (1.0,  1.0, 1.0)),

               'blue':  ((0.0,  0.0, 0.0),
                         (1.0,  0.0, 0.0)), }

green_cmap = colors.LinearSegmentedColormap('custom', green_cdict)


# Cyan
cyan_cdict = {'red':   ((0.0,  0.0, 0.0),
                        (1.0,  0.0, 0.0)),

              'green': ((0.0,  0.0, 0.0),
                        (1.0,  1.0, 1.0)),

              'blue':  ((0.0,  0.0, 0.0),
                        (1.0,  1.0, 1.0)), }

cyan_cmap = colors.LinearSegmentedColormap('custom', cyan_cdict)


# Blue
blue_cdict = {'red':   ((0.0,  0.0, 0.0),
                        (1.0,  0.0, 0.0)),

              'green': ((0.0,  0.0, 0.0),
                        (1.0,  0.0, 0.0)),

              'blue':  ((0.0,  0.0, 0.0),
                        (1.0,  1.0, 1.0)), }

blue_cmap = colors.LinearSegmentedColormap('custom', blue_cdict)


# magenta
magenta_cdict = {'red':   ((0.0,  0.0, 0.0),
                           (1.0,  1.0, 1.0)),

                 'green': ((0.0,  0.0, 0.0),
                           (1.0,  0.0, 0.0)),

                 'blue':  ((0.0,  0.0, 0.0),
                           (1.0,  1.0, 1.0)), }

magenta_cmap = colors.LinearSegmentedColormap('custom', magenta_cdict)


cmaps = {"red": red_cmap,
         "yellow": yellow_cmap,
         "green": green_cmap,
         "cyan": cyan_cmap,
         "blue": blue_cmap,
         "magenta": magenta_cmap}


energy_strings = {0.01 : "10 MeV",
                  0.1 : "100 MeV",
                  10 : "10 GeV",
                  63 : "$m_H/2$",
                  1500 : "1.5 TeV",
                  5000 : "5 TeV",}
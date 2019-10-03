#!/usr/bin/env python3
import multipanel
import pandas as pd
import numpy as np
from matplotlib import cm
import matplotlib.patches as mpatches

import matplotlib as mpl
mpl.use("pgf")
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
         r"\usepackage{units}",         # load additional packages
         r"\usepackage{mathspec}",         # load additional packages
         r"\setmainfont{FreeSans}",         # load additional packages
         r"\setmathsfont(Digits,Latin,Greek)[Uppercase=Italic,Lowercase=Italic]{FreeSans}",
#         r"\setmathsfont{[STIXMath-Regular.otf]}",
#         r"\setmainfont{DejaVu Serif}", # serif font via preamble
         ]
}
mpl.rcParams.update(pgf_with_custom_preamble)
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"

#### relatedness coefficients diploid






##### get the data #####
##### produced by generate_battleground_data.r #####
data = pd.read_csv("battleground_data_vary_df.csv"
        ,sep=";")

gensys = "dip"
data = data.loc[(data["genetic_system"] == gensys)]


######## make the plot ########

the_fig = multipanel.MultiPanel(
        panel_widths=[1,1,1]
        ,panel_heights=[1,1]
        ,filename="effect_philopatry_and_l_" + gensys + ".pdf"
        ,width=15
        ,height=10
        )

offspring_color = "black"
offspring_line = "solid"

maternal_color = "red"
maternal_line = (0,(1,1))

paternal_color = "blue"
paternal_line = (0,(2,2))

padumnal_color = "#009dff"
padumnal_line = (0,(6,4))
lwd_padumnal = 0.75

madumnal_color = "#ff007e"
madumnal_line = (0,(9,1))
lwd_madumnal = 0.75

xlim=[-0.05,1.05]
ylim=[1,5]

lvals = [ 1.0, 0.5 ]
lwd = 0.75

titles = [r"Only local mating, $\ell = 1.0$", r"Half of matings with nonlocal males, $\ell = 0.5"]

for col_i, l in enumerate(lvals):

    the_axis = the_fig.start_block(
            row=0
            ,col=col_i)

    # offspring expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"].isnull())
            ]

    subset = subset.sort_values(by="df")
        
    the_axis.plot(
                    subset["df"]
                    ,subset["af"]
                    ,color=offspring_color
                    ,linewidth=lwd
                    ,linestyle=offspring_line
                    ,label=r"Offspring")

    # maternal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "M")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["af"]
                    ,color=maternal_color
                    ,linewidth=lwd
                    ,linestyle=maternal_line
                    ,label=r"Maternal effect")

    # paternal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "P")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["af"]
                    ,color=paternal_color
                    ,linewidth=lwd
                    ,linestyle=paternal_line
                    ,label=r"Paternal effect")

    # padumnal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "PI")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["af"]
                    ,color=padumnal_color
                    ,linewidth=lwd_padumnal
                    ,linestyle=padumnal_line
                    ,label=r"Paternal imprint")

    # madumnal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "MI")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["af"]
                    ,color=madumnal_color
                    ,linewidth=lwd_madumnal
                    ,linestyle=madumnal_line
                    ,label=r"Maternal imprint")

    title = titles[col_i]

    ylab = ""

    if col_i == 0:
        ylab = r"Female age at maturity, $a_{\text{f}}$"


    if col_i > 0:
        the_axis.legend()

    # end the figure
    the_fig.end_block(
            the_axis
            ,xlim=xlim
            ,ylim=ylim
            ,y_ticks_minor = 5
            ,x_ticks_minor = 5
            ,x_ticks_major_multiple = 0.2
            ,y_ticks_major_multiple = 0.5
            ,xticks=False
            ,yticks=col_i == 0
            ,title=title
            ,ylabel=ylab
            ,xlabel=""
            ,loc_title=True
            ,loc_title_pos=[0.05,0.05]
            )

    the_axis = the_fig.start_block(
            row=1
            ,col=col_i)

    # offspring expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"].isnull())
            ]

    subset = subset.sort_values(by="df")

    assert(subset.shape[0] > 0)
        
    the_axis.plot(
                    subset["df"]
                    ,subset["am"]
                    ,color=offspring_color
                    ,linewidth=lwd
                    ,linestyle=offspring_line
                    ,label=r"Offspring")

    # maternal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "M")
            ]

    print(subset[["df","am"]])

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["am"]
                    ,color=maternal_color
                    ,linewidth=lwd
                    ,linestyle=maternal_line
                    ,label=r"Maternal effect")
    
    # paternal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "P")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["am"]
                    ,color=paternal_color
                    ,linewidth=lwd
                    ,linestyle=paternal_line
                    ,label=r"Paternal")

    # padumnal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "PI")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["am"]
                    ,color=padumnal_color
                    ,linewidth=lwd_padumnal
                    ,linestyle=padumnal_line
                    ,label=r"Padumnal")

    # madumnal expression
    subset = data.loc[
            (data["l"] == l)
            &
            (data["expression"] == "MI")
            ]

    subset = subset.sort_values(by="df")

    the_axis.plot(
                    subset["df"]
                    ,subset["am"]
                    ,color=madumnal_color
                    ,linewidth=lwd_madumnal
                    ,linestyle=madumnal_line
                    ,label=r"Madumnal")

    ylab = ""
    xlab = ""

    if col_i == 0:
        ylab = r"Male age at maturity, $a_{\text{m}}$"

#    if col_i == 1:
#        xlab = r"Female vs male dispersal, $d_{\text{\thinspace f}} = 1 - d_{\text{m}}$"


    # end the figure
    the_fig.end_block(
            the_axis
            ,xlim=xlim
            ,ylim=ylim
            ,y_ticks_minor = 5
            ,x_ticks_minor = 5
            ,x_ticks_major_multiple = 0.2
            ,y_ticks_major_multiple = 0.5
            ,xticks=True
            ,yticks=col_i == 0
            ,title=""
            ,ylabel=ylab
            ,xlabel=xlab
            ,loc_title=True
            ,loc_title_pos=[0.05,0.9]
            )

the_fig.fig.text(x=0.35,
        y=0.06,
        s=r"Female vs male dispersal probability, $d_{\text{f}} = 1 - d_{\text{m}}$",
        fontsize=16,
        verticalalignment="center",
        horizontalalignment="center",
        )


        


the_fig.close(tight=True)

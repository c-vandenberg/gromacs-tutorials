#!/usr/bin/env python3
import sys
import os
import numpy
from numpy import ndarray
from typing import Generator

sys.path.append(
    os.getenv('GROMACS_PYTHON_BASE_PATH')
)

from line_graph import LineGraph


def main():
    base_dir: str = os.path.dirname(__file__)
    md_processed_data_path: str = os.path.join(base_dir, '../../../data/processed/')

    # Ignore comment lines in `mindist.xvg`. Comment lines start with '@', '#' or '&'
    min_dist_gen: Generator = (line for line in open(os.path.join(md_processed_data_path, 'mindist.xvg'))
                               if not line[0] in ('@', '#', '&'))

    min_dist_time: ndarray
    distance: ndarray
    min_dist_time, distance = numpy.genfromtxt(min_dist_gen).T

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack((min_dist_time, distance))],
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Distance$ (nm)',
        x_lim=(0, 1000),
        y_lim=(1, 3),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ Distance\ To\ Periodic\ Boundary\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of coagulation factor Xa distance to periodic boundary during MD '
                    r'simulation, as a function of time.',
        figure_text_font_size=17.5,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=2
    )

    rmsd_time: ndarray
    rmsd: ndarray
    rmsd_time, rmsd = numpy.loadtxt(
        os.path.join(md_processed_data_path, 'rmsd_xray.xvg')
    ).T

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack((rmsd_time, rmsd))],
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ns)',
        y_label=r'$RMSD$ (nm)',
        y_lim=(0.05, 0.225),
        x_lim=(0.0, 1.0),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ Backbone\ RMSD\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 2}$ Evolution of coagulation factor Xa backbone RMSD during MD simulation, as a '
                    r'function of time. Reference structure is energy minimized backbone.',
        figure_text_font_size=15,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=2
    )

    radius_gyration: ndarray = numpy.loadtxt(
        os.path.join(md_processed_data_path, 'gyrate.xvg')
    ).T

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack((radius_gyration[0], radius_gyration[1]))],
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$R_g$ (nm)',
        y_lim=(1.85, 1.95),
        x_lim=(0, 1000),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ Radius\ of\ Gyration\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 3}$ Evolution of coagulation factor Xa radius of gyration ($R_g$) during MD '
                    r'simulation, as a function of time.',
        figure_text_font_size=15,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=2
    )

    h_bonds: ndarray = numpy.loadtxt(
        os.path.join(md_processed_data_path, 'hbnum.xvg')
    ).T

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack((h_bonds[0], h_bonds[1]))],
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$N$ (Hydrogen Bonds)',
        y_lim=(0.0, 5.0),
        x_lim=(-10, 1000),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ Chain\ 1-Chain\ 2\ Hydrogen\ Bonds\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 4}$ Evolution of coagulation factor chain 1-chain 2 hydrogen bonds during MD '
                    r'simulation, as a function of time.',
        figure_text_font_size=15,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=2
    )


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
import sys
import os
import numpy

from typing import List

sys.path.append(
    os.getenv('GROMACS_PYTHON_BASE_PATH')
)

from modules.line_graph import LineGraph


def main():
    base_dir: str = os.path.dirname(__file__)
    em_processed_data_path: str = os.path.join(base_dir, '../../../data/processed/')

    time, potential_energy = numpy.loadtxt(
        os.path.join(em_processed_data_path, 'potential.xvg')
    ).T

    # Multiply time by emstep to get time in picoseconds
    time *= 0.01

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack((time, potential_energy))],
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Potential\ Energy$ (KJ/mol)',
        x_lim=(0, 10),
        y_lim=(-630000, -155000),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ Potential\ Energy\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of coagulation factor Xa system temperature during NVT ensemble '
                    r'equilibration, as a function of time',
        figure_text_font_size=15,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=2
    )


if __name__ == '__main__':
    main()

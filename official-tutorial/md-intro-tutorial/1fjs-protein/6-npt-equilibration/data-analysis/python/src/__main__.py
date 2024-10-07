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
    npt_equilibration_processed_data_path: str = os.path.join(base_dir, '../../../data/processed/')

    time, pressure = numpy.loadtxt(
        os.path.join(npt_equilibration_processed_data_path, 'pressure.xvg')
    ).T

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack([time, pressure])],
        figure_size=(18, 6),
        line_colours=['cyan'],
        x_label='$t$ (ps)',
        y_label=r'$Pressure$ (Bar)',
        x_lim=(0, 100),
        y_lim=(-5250, 500),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ System\ Pressure\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of coagulation factor Xa system pressure during NVT ensemble '
                    r'equilibration, as a function of time.',
        figure_text_font_size=20,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=3
    )

    _, density = numpy.loadtxt(
        os.path.join(npt_equilibration_processed_data_path, 'density.xvg')
    ).T

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack([time, density])],
        figure_size=(18, 6),
        line_colours=['cyan'],
        x_label='$t$ (ps)',
        y_label=r'$Density$ (kg/m$^{3}$)',
        x_lim=(0, 100),
        y_lim=(1000, 1050),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ System\ Density\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 2}$ Evolution of coagulation factor Xa system density during NVT ensemble '
                    r'equilibration, as a function of time.',
        figure_text_font_size=20,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=3,
    )


if __name__ == '__main__':
    main()

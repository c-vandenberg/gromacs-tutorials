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
    time, pressure = numpy.loadtxt(
        '../../pressure.xvg'
    ).T

    # Combine time and potential energy
    pressure_vs_time: List[numpy.ndarray] = [
        numpy.vstack((time, pressure))
    ]

    pressure_vs_time_line_graph = LineGraph()
    pressure_vs_time_line_graph.single_line_graph(
        data_arrays=pressure_vs_time,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Pressure$ (Bar)',
        y_lim=(-5250, 500),
        x_lim=(0, 100),
        graph_title=r'$\bf{Protein\ 1FJS\ Pressure\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of protein 1FJS pressure as a function of time during NPT ensemble '
                    r'equilibration',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )

    _, density = numpy.loadtxt(
        '../../density.xvg'
    ).T

    # Combine time and potential energy
    density_vs_time_data: List[numpy.ndarray] = [
        numpy.vstack((time, density))
    ]

    density_vs_time_line_graph = LineGraph()
    density_vs_time_line_graph.single_line_graph(
        data_arrays=density_vs_time_data,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Density$ (kg/m$^{3}$)',
        y_lim=(1000, 1050),
        x_lim=(0, 100),
        graph_title=r'$\bf{Protein\ 1FJS\ Density\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of protein 1FJS density as a function of time during NPT ensemble '
                    r'equilibration',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()

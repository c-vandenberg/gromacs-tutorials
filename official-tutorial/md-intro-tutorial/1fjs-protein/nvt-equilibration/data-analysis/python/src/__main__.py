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
    time, temperature = numpy.loadtxt(
        '../../temperature.xvg'
    ).T

    # Subtracts the first element from every element in the time array, ensuring the first element is 0/the time starts
    # at zero
    time -= time[0]

    # Combine time and potential energy
    temperature_vs_time_data: List[numpy.ndarray] = [
        numpy.vstack((time, temperature))
    ]

    temperature_vs_timestep_line_graph = LineGraph()
    temperature_vs_timestep_line_graph.single_line_graph(
        data_arrays=temperature_vs_time_data,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Temperature$ (K)',
        y_lim=(295, 305),
        x_lim=(0, 80),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ System\ Temperature\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of coagulation factor Xa system temperature as a function of time during '
                    r'NVT ensemble equilibration',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()

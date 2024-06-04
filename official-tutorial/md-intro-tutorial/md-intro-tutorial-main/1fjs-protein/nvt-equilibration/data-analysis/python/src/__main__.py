#!/usr/bin/env python3
import sys
import os
import pandas
import numpy
import matplotlib.pyplot as pyplot

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
    ave_pot_energy_vs_time: List[numpy.ndarray] = [
        numpy.vstack((time, temperature))
    ]

    ave_pot_energy_vs_timestep_line_graph = LineGraph()
    ave_pot_energy_vs_timestep_line_graph.single_line_graph(
        data_arrays=ave_pot_energy_vs_time,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Temperature$ (K)',
        y_lim=(295, 305),
        x_lim=(0, 80),
        graph_title=r'$\bf{Protein\ 1FJS\ Temperature\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of protein 1FJS temperature as a function of time during NVT ensemble '
                    r'equilibration',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()

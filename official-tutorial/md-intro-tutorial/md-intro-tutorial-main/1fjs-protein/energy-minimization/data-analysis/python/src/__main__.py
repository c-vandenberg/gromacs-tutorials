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
    time, potential_energy = numpy.loadtxt(
        '../../potential.xvg'
    ).T

    # Multiply time by emstep to get time in picoseconds
    time *= 0.01

    # Combine time and potential energy
    ave_pot_energy_vs_time: List[numpy.ndarray] = [
        numpy.vstack((time, potential_energy))
    ]

    ave_pot_energy_vs_timestep_line_graph = LineGraph()
    ave_pot_energy_vs_timestep_line_graph.single_line_graph(
        data_arrays=ave_pot_energy_vs_time,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$Potential\ Energy$ (KJ/mol)',
        y_lim=(-630000, -155000),
        x_lim=(0, 10),
        graph_title=r'$\bf{Protein\ 1FJS\ Potential\ Energy\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of protein 1FJS average potential energy as a function of time '
                    r'during steepest descent minimization',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=2
    )


if __name__ == '__main__':
    main()

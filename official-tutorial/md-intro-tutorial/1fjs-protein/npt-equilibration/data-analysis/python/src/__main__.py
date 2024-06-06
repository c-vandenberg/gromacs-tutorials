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
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ System\ Pressure\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of coagulation factor Xa system pressure as a function of time during '
                    r'NPT ensemble equilibration',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )

    _, density = numpy.loadtxt(
        os.path.join(npt_equilibration_processed_data_path, 'density.xvg')
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
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ System\ Density\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 2}$ Evolution of coagulation factor Xa system density as a function of time during NPT '
                    r'ensemble equilibration',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()

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
    nvt_equilibration_processed_data_path: str = os.path.join(base_dir, '../../../data/processed/')

    time, temperature = numpy.loadtxt(
        os.path.join(nvt_equilibration_processed_data_path, 'temperature.xvg')
    ).T

    # Subtracts the first element from every element in the time array, ensuring the first element is 0/the time starts
    # at zero
    time -= time[0]

    LineGraph.single_line_graph(
        data_arrays=[numpy.vstack((time, temperature))],
        figure_size=(18, 6),
        line_colours=['cyan'],
        x_label='$t$ (ps)',
        y_label=r'$Temperature$ (K)',
        x_lim=(0, 80),
        y_lim=(295, 305),
        graph_title=r'$\bf{Coagulation\ Factor\ Xa\ System\ Temperature\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of coagulation factor Xa system temperature during NVT ensemble '
                    r'equilibration, as a function of time.',
        figure_text_font_size=20,
        figure_text_x_coord=0.5,
        figure_text_y_coord=-0.0005,
        font_size=20,
        tick_label_size=20,
        line_width=3
    )


if __name__ == '__main__':
    main()

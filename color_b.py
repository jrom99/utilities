"""
DESCRIPTION:
color_b selection='ca',mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1

USAGE:
color_b selection[, gradient]

PARAMS:
selection (PyMOL selection)
    name of the PyMOL object to color according to b_factor

gradient (string; defaults to red_white_blue)
    name of the gradient
"""

from pymol import cmd
import re

_GRADIENTS = {
    "bgr": "blue green red",
    "rgb": "red green blue",
    "bwr": "blue_white_red",
    "rwb": "red_white_blue",
    "bmr": "blue magenta red",
    "rmb": "red magenta blue",
    "rw": "red white",
    "wr": "white red",
    "ry": "red_yellow",
    "yr": "yellow_red",
    "gw": "green white",
    "wg": "white green",
    "gy": "green_yellow",
    "yg": "yellow_green",
    "gray": "black white",
    "gray_rev": "white black",
    "blue_green": "blue_green",
    "blue_magenta": "blue_magenta",
    "blue_red": "blue_red",
    "blue_white_green": "blue_white_green",
    "blue_white_magenta": "blue_white_magenta",
    "blue_white_red": "blue_white_red",
    "blue_white_yellow": "blue_white_yellow",
    "blue_yellow": "blue_yellow",
    "cbmr": "cbmr",
    "cyan_magenta": "cyan_magenta",
    "cyan_red": "cyan_red",
    "cyan_white_magenta": "cyan_white_magenta",
    "cyan_white_red": "cyan_white_red",
    "cyan_white_yellow": "cyan_white_yellow",
    "cyan_yellow": "cyan_yellow",
    "gcbmry": "gcbmry",
    "green_blue": "green_blue",
    "green_magenta": "green_magenta",
    "green_red": "green_red",
    "green_white_blue": "green_white_blue",
    "green_white_magenta": "green_white_magenta",
    "green_white_red": "green_white_red",
    "green_white_yellow": "green_white_yellow",
    "green_yellow": "green_yellow",
    "green_yellow_red": "green_yellow_red",
    "magenta_blue": "magenta_blue",
    "magenta_cyan": "magenta_cyan",
    "magenta_green": "magenta_green",
    "magenta_white_blue": "magenta_white_blue",
    "magenta_white_cyan": "magenta_white_cyan",
    "magenta_white_green": "magenta_white_green",
    "magenta_white_yellow": "magenta_white_yellow",
    "magenta_yellow": "magenta_yellow",
    "rainbow": "rainbow",
    "rainbow2": "rainbow2",
    "rainbow2_rev": "rainbow2_rev",
    "rainbow_cycle": "rainbow_cycle",
    "rainbow_cycle_rev": "rainbow_cycle_rev",
    "rainbow_rev": "rainbow_rev",
    "red_blue": "red_blue",
    "red_cyan": "red_cyan",
    "red_green": "red_green",
    "red_white_blue": "red_white_blue",
    "red_white_cyan": "red_white_cyan",
    "red_white_green": "red_white_green",
    "red_white_yellow": "red_white_yellow",
    "red_yellow": "red_yellow",
    "red_yellow_green": "red_yellow_green",
    "rmbc": "rmbc",
    "yellow_blue": "yellow_blue",
    "yellow_cyan": "yellow_cyan",
    "yellow_cyan_white": "yellow_cyan_white",
    "yellow_green": "yellow_green",
    "yellow_magenta": "yellow_magenta",
    "yellow_red": "yellow_red",
    "yellow_white_blue": "yellow_white_blue",
    "yellow_white_green": "yellow_white_green",
    "yellow_white_magenta": "yellow_white_magenta",
    "yellow_white_red": "yellow_white_red",
    "yrmbcg": "yrmbcg",
}


def color_b(selection, gradient="red_white_blue"):
    """
    DESCRIPTION:
        color_b selection='ca',mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1

    USAGE:
        color_b selection[, gradient]

    PARAMS:
        selection (PyMOL selection)
            name of the PyMOL object to color according to b_factor

        gradient (string; defaults to red_white_blue)
            name of the gradient
    """

    if gradient not in _GRADIENTS:
        raise KeyError(f"{gradient!s} not found")
    cmd.spectrum("b", _GRADIENTS[gradient], selection)


cmd.extend("color_b", color_b)
cmd.color_b = color_b
cmd.auto_arg[0]['color_b'] = [ cmd.object_sc, 'object', '']
cmd.auto_arg[1]['color_b'] = [lambda: cmd.Shortcut([*_GRADIENTS]), 'params', '']

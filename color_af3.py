"""
USAGE:
color_af3 selection

PARAMS:
selection (PyMOL selection)
    name of the PyMOL object to color according to b_factor
"""

from pymol import cmd
import re


def color_af3(selection, use_ca="no"):
    """USAGE:
        color_af3 selection[, use_ca=1]
    """
    if str(use_ca) in ("yes", "on", "1"):
        prefix = "name CA and "
    elif str(use_ca) in ("no", "off", "0"):
        prefix = ""
    else:
        raise ValueError("Invalid boolean value for argument 'use_ca'")

    cmd.color("0xef821e", f"{prefix}({selection})")
    cmd.color("0xf6ed12", f"{prefix}({selection}) and b > 50")
    cmd.color("0x10cff1", f"{prefix}({selection}) and b > 70")
    cmd.color("0x106dff", f"{prefix}({selection}) and b > 90")


cmd.extend("color_af3", color_af3)
cmd.color_af3 = color_af3
cmd.auto_arg[0]['color_af3'] = [ cmd.object_sc, 'object', '']
cmd.auto_arg[1]['color_af3'] = [ lambda: cmd.Shortcut(["off", "on"]), 'params', '']

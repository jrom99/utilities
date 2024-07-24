"""
USAGE:
color_af3 selection

PARAMS:
selection (PyMOL selection)
    name of the PyMOL object to color according to b_factor
"""

from pymol import cmd
import re

def color_af3(selection):
    """USAGE:
        color_af3 selection
    """
    cmd.color("0xef821e", f"name CA and ({selection})")
    cmd.color("0xf6ed12", f"name CA and ({selection}) and b > 50")
    cmd.color("0x10cff1", f"name CA and ({selection}) and b > 70")
    cmd.color("0x106dff", f"name CA and ({selection}) and b > 90")

cmd.extend("color_af3", color_af3)
cmd.color_af3 = color_af3
cmd.auto_arg[0]['color_af3'] = [ cmd.object_sc, 'object', '']



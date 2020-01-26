#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>

from pymol import cmd

cmd.set_key("CTRL-W", cmd.reinitialize)
cmd.set_key("CTRL-D", cmd.quit)

# coding: utf8
"""
                            PyMol Toolbox
"""

import pymol
from pymtools.commands import PymToolCommand
from pymtools.base import CustomLogging

pymol.pymol_argv = ['pymol', '-qxc']


def main():
    """
    Pymol toolbox
    """
    mylog = CustomLogging(desc=__doc__)

    command = PymToolCommand(custom_logging=mylog)

    command.run()

if __name__ == '__main__':
    main()

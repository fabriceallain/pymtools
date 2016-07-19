# coding: utf8
"""
                            PyMol Toolbox
"""

import os
import pymol
import logging
import argparse
from pymtools.commands import PymToolCommand
from pymtools.colorbyrmsd import colorbyrmsd
from pymtools.average3d import avg_states
from pymtools.base import CustomLogging
from pymol import cmd, CmdException

pymol.pymol_argv = ['pymol', '-qxc']


def settings(desc):
    """
    Parse args and options
    :return: 2-tuple (opts, args)
    """
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-o", "--output", default="", help="Output path")
    parser.add_argument("--pretty", dest="pretty", action="store_true",
                        default=False,
                        help="Nice representation and colors")
    parser.add_argument("--prefix", help="output prefix", default="prot")

    argaver = argparse.ArgumentParser(add_help=False)

    argaver.add_argument("pdb", help="pdb structure with several states")
    argaver.add_argument("--sel", default="name CA",
                         help="Subset of atoms on which to operate (e.g., "
                              "'name CA and resi 20-50')")
    argaver.add_argument("--states", default=[1, 0],
                         help="number of the first state to include in "
                              "averaging and number of the last state to "
                              "include in averaging. A value of '0' (zero) "
                              "means use the last state in the object.")
    argaver.add_argument("--fit", choices=['yes', 'no'], default='no',
                         help="Use the pymol function intra_fit() to fit all "
                              "the states of 'object & object_sel'")
    argaver.add_argument("--pair", choices=[0, 1], default=1,
                         help="Calculate average pairwise RMSD for each residue"
                              " position")

    argvisu = argparse.ArgumentParser(add_help=False)

    argvisu.add_argument("struct", help="structure file [.pdb]")

    argalign = argparse.ArgumentParser(add_help=False)

    argalign.add_argument("target", help="target structure file [.pdb]")
    argalign.add_argument("mob",
                          help="mobile structure file [.pdb]")
    argalign.add_argument("--align", dest="align", type=int, choices=(0, 1),
                          default=1,
                          help="superpose structures before RMSD calculation")
    argalign.add_argument("--guide", dest="guide", default=False,
                          action="store_true",
                          help="Only use C-alpha atoms")
    argalign.add_argument("--method", dest="method", default="tmalign",
                          choices=('super', "align", 'tmalign', 'cealign'),
                          help="Method to match atoms")

    sp = parser.add_subparsers(dest='command')

    sp.add_parser('align', parents=[argalign],
                  help="Align 2 pdb structures")

    sp.add_parser('visu', parents=[argvisu],
                  help="Pymol visualisation")

    sp.add_parser('aver', parents=[argaver],
                  help="Average 3d states of a pdb structure")

    args = parser.parse_args()

    return args, parser.prog


def pympretty(rainbow=False):
    """
    High quality color scheme for pymol
    :return:
    """
    pymol.cmd.show_as('cartoon')
    if rainbow:
        pymol.cmd.spectrum()
    pymol.cmd.set("ray_trace_mode", 3)
    pymol.cmd.bg_color("white")
    pymol.cmd.set("antialias", 2)

    # Remove ligand and heteroatoms
    pymol.cmd.remove("resn HOH")
    pymol.cmd.remove("resn HET")
    pymol.cmd.orient()
    pymol.cmd.center()
    # pymol.cmd.ray(800, 800)


def main():
    """
    Pymol toolbox
    """
    mylog = CustomLogging(desc=__doc__)

    command = PymToolCommand(custom_logging=mylog)

    command.run()

    # if args.command == "align":
    #
    #     logger.info("Load %s target file", args.target)
    #     try:
    #         cmd.load(args.target)
    #     except CmdException:
    #         logger.error("Can't load %s", args.target)
    #
    #     logger.info("Load %s mobile file", args.mob)
    #     try:
    #         cmd.load(args.mob)
    #     except CmdException:
    #         logger.error("Can't load %s", args.mob)
    #
    #     logger.info("Align and color by RMSD mobile against target structure")
    #     mob = os.path.basename(os.path.splitext(args.mob)[0])
    #     target = os.path.basename(os.path.splitext(args.target)[0])
    #     logger.debug("mobile prefix: %s", mob)
    #     logger.debug("target prefix: %s", target)
    #     colorbyrmsd(mob, target, quiet=0, guide=1 if args.guide else 0,
    #                 method=args.method)
    # elif args.command == "visu":
    #     logger.info("Load %s mobile file", args.struct)
    #     try:
    #         cmd.load(args.struct)
    #     except CmdException:
    #         logger.error("Can't load %s", args.struct)
    # elif args.command == "aver":
    #     cmd.set('orthoscopic', 1)
    #     # Show all states ?
    #     # cmd.set('all_states', 1)
    #     logger.info("Load structure file %s" % args.pdb)
    #     try:
    #         cmd.load(args.pdb)
    #     except CmdException:
    #         logger.error("Can't load %s", args.pdb)
    #     mypdb = os.path.splitext(os.path.basename(args.pdb))[0]
    #     cmd.hide('everything', mypdb)
    #     cmd.show('ribbon', mypdb)
    #     cmd.show('spheres', '%s and name CA' % mypdb)
    #     cmd.alter(mypdb, 'vdw=0.5')
    #     cmd.rebuild(mypdb, 'spheres')
    #     # Fit to the first state by default (which is the best one according to
    #     #  aria order)
    #     objfit = avg_states(mypdb, object_sel=args.sel, first=1, last=0,
    #                         newobj='fit_ALL', fitverdict=args.fit, writefiles=0)
    #     cmd.cartoon('putty', objfit)
    #     cmd.show('cartoon', objfit)
    #     cmd.spectrum('b', 'blue_white_red', objfit)
    # elif args.command == "cnstr":
    #     # Check if violations file (csv) or tbl file or contact file (txt)
    #     # If violations file, show restraints with colors according to
    #     # contact_5, contact_8, viol, valid, group labels (argsparse option)
    #     # Otherwise, fetch contribution list and plot contacts or restraints
    #     # In this case, we can give an extra pdb file correspondint to the
    #     # native state to label
    #     pass
    #
    # if args.pretty and args.command != 'aver':
    #     logger.info("Ray tracing pymol png")
    #     spec = True if args.command == "visu" else False
    #     pympretty(rainbow=spec)
    #
    # logger.info("Saving pymol png and session")
    # pymol.cmd.png("%s.png" % os.path.join(args.output, args.prefix))
    # pymol.cmd.save("%s.pse" % os.path.join(args.output, args.prefix))
    #
    # logger.info("Closing pymol")
    # pymol.cmd.quit()

if __name__ == '__main__':
    main()

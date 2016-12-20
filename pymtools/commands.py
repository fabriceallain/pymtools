# coding=utf-8
"""
                        Command line class
"""
import os
import pymol
import logging
import argparse as argp
from .base import Settings
from pymol import cmd, CmdException
from .colorbyrmsd import colorbyrmsd
from .average3d import avg_states
from .cnstr import Cnstr
from .align_all import align_all


LOG = logging.getLogger(__name__)


def check_file(prospective_file):

    """

    :param prospective_file:
    """
    if not os.path.exists(prospective_file):
        raise argp.ArgumentTypeError("readable_file:'{0}' is not a valid "
                                     "path".format(prospective_file))
    if not os.access(prospective_file, os.R_OK):
        raise argp.ArgumentTypeError("readable_file:'{0}' is not a readable "
                                     "file".format(prospective_file))


class FileChecker(argp.Action):

    """
    File checker for argparse action
    """

    def __init__(self, *args, **kwargs):
        super(FileChecker, self).__init__(*args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if type(values) == list:
            for prospective_file in values:
                check_file(prospective_file)
        elif type(values) == str:
            check_file(values)
        setattr(namespace, self.dest, values)


class PymToolSettings(Settings):

    """
    Main settings class
    """
    DESC = u"PyMol Toolbox"
    COMMANDS = ("align", "visu", "aver", "cnstr", "multialign")
    COMMANDESCS = (u"Align 2 pdb structures",
                   u"Pymol visualisation",
                   u"Average 3d states of a pdb structure",
                   u"Visualize constrains on a pdb structure",
                   u"Align multiple structures to one reference state")
    CONFPATH = "conf/pymtools.ini"

    def __init__(self):
        super(PymToolSettings, self).__init__(self.COMMANDS)
        self.load_config(self.CONFPATH, pkg=True)

    def load_config(self, configpath, **kwargs):
        """
        Load configuration file
        :param configpath:
        :param kwargs:
        """
        LOG.debug("Loading default config file")
        super(PymToolSettings, self).load_config(configpath, **kwargs)


class Command(object):
    """
    Argparse interface enable multiple commands
    """

    def __init__(self, log=None):
        self._settings = None
        parser = self._create_argparser()
        self._args = parser.parse_args()
        if log and hasattr(log, 'welcome'):
            log.welcome()
        self.command = self._args.command
        self._update_logger(log)
        self._update_settings()

    @property
    def settings(self):
        """
        settings
        """
        raise NotImplementedError

    def _create_argparser(self):
        parser = argp.ArgumentParser(
            description=self.settings.DESC,
            formatter_class=argp.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-o", "--output",
                            type=str, help="Output directory", required=True)
        parser.add_argument("-c", "--conf", action=FileChecker,
                            dest="conf_file", default=None,
                            help="configuration file")
        parser.add_argument("--log", action="store_true",
                            default=False, help="Generate log files")
        # Create subcommands
        self._create_subparsers(parser.add_subparsers(dest="command"))
        return parser

    def _create_subparsers(self, parser):
        """
        Generate subcommands
        :param parser: argparser object
        :return:
        """
        for index, command in enumerate(self.settings.COMMANDS):
            # Create subparser defined in command list
            subcommand = getattr(self, "_" + command + "_parser")(
                desc=self.settings.COMMANDESCS[index])
            parser.add_parser(command, parents=[subcommand])

    def _update_logger(self, log):
        if log:
            if hasattr(self._args, "output"):
                if not self._args.log:
                    # Don't generate log files
                    LOG.removeHandler("info_file_handler")
                    LOG.removeHandler("error_file_handler")
                    LOG.removeHandler("debug_file_handler")
                log.set_outdir(self._args.output)

    def _update_settings(self):
        if self._args.conf_file:
            LOG.info("Updating settings with conf file")
            self.settings.load_config(self._args.conf_file)
        # Update settings associated to command section
        LOG.info("Updating %s settings with args" % self.command)
        getattr(self.settings, self.command).args.update(self._args.__dict__)

    def run(self):
        """
        run command
        """
        LOG.info("Run %s command" % self.command)
        getattr(self, self.command)(getattr(self.settings, self.command))


class PymToolCommand(Command):
    """
    Argparse interface for pymtools
    """

    @property
    def settings(self):
        """
        settings
        :return:
        """
        if not self._settings:
            self._settings = PymToolSettings()
        return self._settings

    def __init__(self, custom_logging=None):
        super(PymToolCommand, self).__init__(custom_logging)

    def _create_argparser(self):
        parser = super(PymToolCommand, self)._create_argparser()
        parser.add_argument("--prefix", help="output prefix", default="prot")
        return parser

    @staticmethod
    def _align_parser(desc=None):
        parser = argp.ArgumentParser(description=desc,
                                     add_help=False)
        group = parser.add_argument_group('required arguments')
        group.add_argument("target", help="target structure file [.pdb]")
        group.add_argument("mob",
                           help="mobile structure file [.pdb]")
        group.add_argument("--pretty", dest="pretty", action="store_true",
                           default=False,
                           help="Nice representation and colors")
        return parser

    @staticmethod
    def _multialign_parser(desc=None):
        parser = argp.ArgumentParser(description=desc,
                                     add_help=False)
        group = parser.add_argument_group('required arguments')
        group.add_argument("target", help="target structure file [.pdb]")
        group.add_argument("mob",
                           help="Structure file(s) [.pdb]",
                           nargs="+")
        group.add_argument("--pretty", dest="pretty", action="store_true",
                           default=False,
                           help="Nice representation and colors")
        return parser

    @staticmethod
    def _visu_parser(desc=None):
        parser = argp.ArgumentParser(description=desc,
                                     add_help=False)
        group = parser.add_argument_group('required arguments')
        group.add_argument("struct", help="structure file [.pdb]",
                           action=FileChecker)
        group.add_argument("--pretty", dest="pretty", action="store_true",
                           default=False,
                           help="Nice representation and colors")
        return parser

    @staticmethod
    def _aver_parser(desc=None):
        parser = argp.ArgumentParser(description=desc,
                                     add_help=False)
        group = parser.add_argument_group('required arguments')
        group.add_argument("pdb", help="pdb structure with several states")
        group.add_argument("--sel", default="name CA",
                           help="Subset of atoms on which to operate (e.g., "
                                "'name CA and resi 20-50')")
        group.add_argument("--states", default=[1, 0],
                           help="number of the first state to include in "
                                "averaging and number of the last state to "
                                "include in averaging. A value of '0' (zero) "
                                "means use the last state in the object.")
        return parser

    @staticmethod
    def _cnstr_parser(desc=None):
        parser = argp.ArgumentParser(description=desc, add_help=False)

        group = parser.add_argument_group('required arguments')
        group.add_argument(
            "constrains", action=FileChecker,
            help="NMR constrains in .tbl, .txt or in .csv format. For .txt and "
                 ".csv format, the file should follow the same syntax as in "
                 "examples/data folder.")
        group.add_argument("pdb", action=FileChecker, help="PDB file")
        group.add_argument("-r", "--ref", action=FileChecker,
                           help="Reference pdb file")
        group.add_argument("-g", "--group", dest="group", default=None,
                           help="Name of the field if csv file given")
        group.add_argument("-w", "--writefiles", dest="writefiles",
                           action="store_true", default=False,
                           help="Write output txt files")
        return parser

    @staticmethod
    def align(settings):
        """
        align command
        :param settings:
        """
        LOG.info("Load %s target file", settings.args.get('target'))
        try:
            cmd.load(settings.args.get('target'))
        except CmdException:
            LOG.error("Can't load %s", settings.args.get('target'))

        LOG.info("Load %s mobile file", settings.args.get('mob'))
        try:
            cmd.load(settings.args.get('mob'))
        except CmdException:
            LOG.error("Can't load %s", settings.args.get('mob'))
            raise

        LOG.info("Align and color by RMSD mobile against target structure "
                 "with method %s" % settings.config.get('method'))
        mob = os.path.basename(os.path.splitext(settings.args.get('mob'))[0])
        target = os.path.basename(os.path.splitext(
            settings.args.get('target'))[0])
        LOG.debug("mobile prefix: %s", mob)
        LOG.debug("target prefix: %s", target)
        colorbyrmsd(mob, target, quiet=0,
                    guide=settings.config.get("guide"),
                    method=settings.config.get('method'))

    @staticmethod
    def multialign(settings):
        """
        align command
        :param settings:
        """
        LOG.info("Load %s target file", settings.args.get('target'))
        try:
            cmd.load(settings.args.get('target'))
        except CmdException:
            LOG.error("Can't load %s", settings.args.get('target'))

        LOG.info("Load %s structure file(s)", settings.args.get('mob'))
        for pdbfile in settings.args.get('mob'):
            try:
                cmd.load(pdbfile)
            except CmdException:
                LOG.error("Can't load %s", pdbfile)
                raise

        target = os.path.basename(os.path.splitext(
            settings.args.get('target'))[0])
        LOG.debug("target prefix: %s", target)
        align_all(target, method=settings.config.get('method'))

    @staticmethod
    def visu(settings):
        """
        visu command
        :param settings:
        """
        LOG.info("Load %s mobile file", settings.args.get('struct'))
        try:
            cmd.load(settings.args.get('struct'))
        except CmdException:
            LOG.error("Can't load %s", settings.args.get('struct'))

    @staticmethod
    def aver(settings):
        """
        average pdb command
        :param settings:
        """
        cmd.set('orthoscopic')
        if settings.config.get('all_states', 0):
            cmd.set('all_states')
        LOG.info("Load structure file %s" % settings.args.get('pdb'))
        try:
            cmd.load(settings.args.get('pdb'))
        except CmdException:
            LOG.error("Can't load %s", settings.args.get('pdb'))
        mypdb = os.path.splitext(os.path.basename(settings.args.get('pdb')))[0]
        cmd.hide('everything', mypdb)
        cmd.show('ribbon', mypdb)
        cmd.show('spheres', '%s and name CA' % mypdb)
        cmd.alter(mypdb, 'vdw=%.2f' % settings.config.get('vdwspheres_radius'))
        cmd.rebuild(mypdb, 'spheres')
        # Fit to the first state by default (which is the best one according to
        #  aria order)
        objfit = avg_states(mypdb, object_sel=settings.args.get('sel'),
                            first=settings.config.get('first'),
                            last=settings.config.get('last'),
                            newobj=settings.config.get('newobj'),
                            fitverdict=settings.config.get('fit'),
                            verb=settings.config.get('verb'),
                            pairs=settings.config.get('pairs'),
                            writefiles=settings.config.get('writefiles'))
        cmd.cartoon('putty', objfit)
        cmd.show('cartoon', objfit)
        cmd.spectrum('b', 'blue_white_red', objfit)

    @staticmethod
    def cnstr(settings):
        """
        constraints visu command
        :param settings:
        """
        Cnstr(settings=settings).run()

    @staticmethod
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

    def save(self, settings):
        """
        Save pymol session
        :param settings:
        """
        LOG.info("Saving pymol png and session")
        if settings.args.get('pretty'):
            LOG.info("Ray tracing pymol png")
            spec = True if self.command == "visu" else False
            self.pympretty(rainbow=spec)
        pymol.cmd.png("%s.png" % os.path.join(settings.args.get('output'),
                                              settings.args.get('prefix')))
        pymol.cmd.save("%s.pse" % os.path.join(settings.args.get('output'),
                                               settings.args.get('prefix')))

    def run(self):
        # Avoid GUI opening
        """
        Run & avoid GUI opening
        """
        pymol.finish_launching()
        # call method relative to args.command
        super(PymToolCommand, self).run()
        # Save image and session
        self.save(getattr(self.settings, self.command))
        LOG.info("Closing pymol")
        pymol.cmd.quit()

# coding: utf8
"""
                            Base functions
"""
import os
import re
import ast
import json
import logging
import collections
import logging.config
import pkg_resources as pkgr
from configparser import SafeConfigParser


logger = logging.getLogger(__name__)


class Setting:

    def __init__(self):
        self.config = collections.defaultdict()
        self.args = collections.defaultdict()

    def __repr__(self):
        return "Setting object\n    config: %s\n    args  : %s" % (self.config,
                                                                   self.args)


class Settings(object):

    def __init__(self, sections):
        self._sections = set(sections)
        for section in self._sections:
            setattr(self, section, Setting())

    def load_config(self, configpath, pkg=False):
        """
        Use ConfigParser module to load config sections
        :param pkg: file is inside the package (configpath is the relative path
        inside the package)
        :param configpath:
        :return:
        """
        if not pkg and not os.path.exists(configpath):
            logger.error("Configuration file not found (%s)" % configpath)
            from errno import ENOENT
            raise OSError(ENOENT)
        config = SafeConfigParser(allow_no_value=True)
        if pkg:
            with pkgr.resource_stream(__name__, configpath) as conf:
                config.readfp(conf)
        else:
            config.read(configpath)
        logger.debug(config)
        for section in config.sections():
            if hasattr(self, section):
                tmp = format_dict(dict(config.items(section)))
                getattr(self, section).config.update(tmp)
                logger.debug("%s config updated" % section)
                logger.debug("%s.%s : %s" % (self.__class__.__name__, section,
                                             getattr(self, section)))
            else:
                logger.warning("Unknow config section %s" % section)

    def write_config(self, filename):
        # Ecrit les config de toutes les sections dans un autre fichier
        logger.info("Writing .ini file (%s)" % filename)
        config = SafeConfigParser(allow_no_value=True)
        iniout = open(filename, mode="w")
        for section in self._sections:
            config.add_section(section)
            if hasattr(self, section):
                for opt in getattr(self, section).config:
                    config.set(section, str(opt),
                               str(getattr(self, section).config.get(opt)))
        config.write(iniout)

    def __repr__(self):
        return "<Settings object>\n    sections: %s" % self._sections


class CustomLogging(object):
    """
    Custom configuration for logging
    """
    # default_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),
    #                             "conf/logging.json")
    default_file = "conf/logging.json"

    def __init__(self, level=logging.INFO, desc=None):
        """

        :param level:
        :param desc:
        :return:
        """
        # TODO: detect path log filenames and makedirs if not exists
        logging.basicConfig(level=level)
        if desc:
            self.msg = desc.strip()
        else:
            self.msg = ""
        self.config = self.default_config()

    def update_msg(self, desc):
        """

        :param desc:
        :return:
        """
        if type(self.msg) == list:
            self.msg += desc
            self.msg = " - ".join(self.msg)
        elif type(self.msg) == str:
            self.msg = " - ".join((self.msg, desc.capitalize()))

    def default_config(self):
        """

        :return:
        """
        # with open(self.default_file, 'rt') as f:
        with pkgr.resource_stream(__name__, self.default_file) as f:
            config = json.load(f)
            logging.config.dictConfig(config)
        return config

    def set_outdir(self, outdir):
        """
        Create log directory and change log files location
        :param outdir: path output directory
        :return:
        """
        outdir = os.path.join(outdir, "log") if "log" not in outdir else outdir
        if not os.path.exists(os.path.abspath(outdir)):
            os.makedirs(outdir)
        if outdir and "handlers" in self.config:
            for hand in self.config["handlers"]:
                if "filename" in self.config["handlers"][hand]:
                    self.config["handlers"][hand]["filename"] = \
                        os.path.join(outdir, os.path.basename(
                            self.config["handlers"][hand]["filename"]))
            logging.config.dictConfig(self.config)

    def welcome(self):
        """

        :return:
        """
        desc = '''
================================================================================

{:^80}

================================================================================
'''.format(self.msg)
        print(desc)


def format_str(string):
    """
    Convert str in bool, float, int or str
    :param string:
    :return:
    """
    if re.search(r"^\s*(true)\s*$", string, re.I):
        return True
    elif re.search(r"^\s*(false)\s*$", string, re.I):
        return False
    elif re.search(r"^\s*\d+\s*$", string):
        return int(string)
    elif re.search(r"^[\s\d-]+\.\d+\s*$", string):
        return float(string)
    elif "," in string:
        return string.split(',')
    elif "+" in string:
        return string.split('+')
    elif re.search(r"[/\w]+", string):
        return string
    else:
        if string:
            try:
                ev_str = ast.literal_eval(string)
            except ValueError:
                logger.error("Don't understand given string %s. Please check "
                             "format." % string)
                return None
            except SyntaxError:
                logger.error("Given string %s is not a valid expression" % string)
                return None
            return ev_str
        else:
            return None


def format_dict(indict):
    """

    :param indict:
    :return:
    """
    for key in indict:
        indict[key] = format_str(indict[key])
    return indict
"""
                            Pymol module reading NMR restraints
"""

from __future__ import print_function

import sys
import os
import logging
from pymol import cmd, CmdException
from .reader import CnstrFile


logger = logging.getLogger(__name__)


class Cnstr(object):
    """
    Read restraints and show them on pdb struct
    """
    def __init__(self, settings):
        self.settings = settings
        self._pdbname = os.path.basename(os.path.splitext(
            self.settings.args.get('pdb'))[0])
        self._cnstrfile = CnstrFile(self.settings.args.get('constrains'))
        self._reflag = False
        self._refname = None

    def load_cnstr(self):
        self._cnstrfile.load()
        names = []
        eval_method = self.settings.config.get('eval_method', 'contact')
        contact_type = self.settings.config.get('contact_type')

        if self.settings.args.get('writefiles', False):
            outfile = open("%s.out" % os.path.join(
                self.settings.args.get('output'),
                self.settings.args.get('prefix')), 'w')
            extra = "\t".join(("lower", "target", "upper")) if \
                eval_method == "bound" else ""
            outfile.write("\t".join(
                ("resid1", "resid2", "atm1", "atm2",
                 "dist_ref" if self._reflag else "dist", eval_method,
                 extra, "spec")))
        else:
            outfile = None

        for idline in self._cnstrfile.lines:
            # For each restraint
            line = self._cnstrfile.lines[idline]
            spec = line.get('spec', 'noname')
            resids = (line.get('resid1', ''), line.get('resid2', ''))
            atms = (line.get('atm1', ''), line.get('atm2', ''))
            target = line.get('target', None)
            upper = line.get('upper', None)
            lower = line.get('lower', None)
            group = line.get('group', None)

            if target:
                upper = float(target) + float(upper) if upper else None
                lower = float(target) - float(lower) if lower else None

            extra = "\t".join((str(lower), str(target), str(upper))) if \
                eval_method == "bound" else ""
            # TODO: add a violation tolerance ?

            if '' in atms:
                logger.warning('No atms related to cnsrt no %d' % idline)
            if '' in resids:
                logger.error('No cnstr at line %d' % idline)
                continue

            # TODO: condition below can be simplified ...
            if contact_type == "all":
                sel1 = '%s and resi %s and name %s' % (
                    self._pdbname, resids[0], atms[0]) if atms[0] \
                    else '%s and resi %d' % (self._pdbname, resids[0])
                sel2 = '%s and resi %s and name %s' % (
                    self._pdbname, resids[1], atms[1]) if atms[1] \
                    else '%s and resi %d' % (self._pdbname, resids[1])
                sel3 = '{ref} and resi {resid1} and name {atm1}'.format(
                    resid1=resids[0], atm1=atms[0], ref=self._refname)
                sel4 = '{ref} and resi {resid2} and name {atm2}'.format(
                    resid2=resids[1], atm2=atms[1], ref=self._refname)
            else:
                sel1 = '{pdb} and resi {resid} and name CA'.format(
                    pdb=self._pdbname, resid=resids[0])
                sel2 = '{pdb} and resi {resid} and name CA'.format(
                    pdb=self._pdbname, resid=resids[1])
                sel3 = '{ref} and resi {resid1} and name CA'.format(
                    resid1=resids[0], ref=self._refname)
                sel4 = '{ref} and resi {resid2} and name CA'.format(
                    resid2=resids[1], ref=self._refname)

            if self._reflag:
                # If structure given, we look at reference distance
                dist = cmd.get_distance(sel3, sel4)
            else:
                dist = cmd.get_distance(sel1, sel2)

            # TODO: condition below can be simplified
            # Group ca be decided by a contact treshold in the structure or in
            # the reference (if given). It can also be if distances in the
            # structure (or ref if given) lie inside bounds defined in the
            # cnstr file
            if not group:
                if eval_method == 'contact' and self._reflag:
                    group = "TP" if dist < float(self.settings.config.get(
                        'contact', 5.0)) else "FP"
                elif eval_method == 'contact' and not self._reflag:
                    group = "unviol" if dist < float(self.settings.config.get(
                        'contact', 5.0)) else "viol"
                elif eval_method == "bound" and self._reflag:
                    if lower and upper:
                        group = "TP" if lower <= dist <= upper else "FP"
                    else:
                        logger.error("No bound available for restraint no %d" %
                                     idline)
                        continue
                elif eval_method == "bound" and not self._reflag:
                    if lower and upper:
                        group = "unviol" if lower <= dist <= upper else "viol"
                    else:
                        logger.error("No bound available for restraint no %d" %
                                     idline)
                        continue
                else:
                    logger.critical(
                        "Evaluation method given is not supported (%s)"
                        % eval_method)
                    sys.exit(1)

            name = "_".join(('mobile', spec, group))
            nameref = "_".join(('ref', spec, group)) if self._reflag \
                else None
            if name not in names:
                names.append(name)
                if self._reflag:
                    names.append(nameref)
            try:
                cmd.distance(name, sel1, sel2)
                if self._reflag:
                    cmd.distance(nameref, sel3, sel4)
            except CmdException:
                logger.error("Wrong selection")

            if outfile:
                if contact_type == "all" or (
                                contact_type == "min" and atms == ("CA", "CA")):
                    outfile.write("\n%s" % "\t".join((resids[0], resids[1],
                                                      atms[0], atms[1],
                                                      str(dist), group, extra,
                                                      spec)))

        for name in names:
            if name.endswith("_FP") or name.endswith("_viol"):
                color = self.settings.config.get("fp_color", 'red')
                cmd.color(color, name)
                if 'mobile' in name and self._reflag:
                    # If reference structure given, we look first at reference
                    # contacts
                    cmd.disable(name)
            else:
                color = self.settings.config.get("tp_color", 'green')
                cmd.color(color, name)
                if 'mobile' in name and self._reflag:
                    cmd.disable(name)

        cmd.show_as(self.settings.config.get('view', 'cartoon'))
        cmd.dss()
        cmd.show('dashes')
        if contact_type == "all":
            cmd.show('sticks')
            cmd.set('stick_radius', self.settings.config.get('stick_radius',
                                                             0.05))
            cmd.set('dash_width', self.settings.config.get('dash_width', 0.5))
            cmd.set('dash_gap', self.settings.config.get('dash_gap', 0.5))
            cmd.set('stick_transparency', self.settings.config.get('stick_alpha',
                                                                   0.8))
            # cmd.hide('labels')
        cmd.color(self.settings.config.get('color', 'gray'), self._pdbname)
        if self._reflag:
            cmd.spectrum(selection=self._refname)

        if outfile:
            outfile.close()

    def run(self):
        logger.info("Load structure file %s" % self.settings.args.get('pdb'))
        try:
            cmd.load(self.settings.args.get('pdb'), object=self._pdbname)
        except CmdException:
            logger.error("Can't load %s", self.settings.args.get('pdb'))
        if self.settings.args.get('ref'):
            self._refname = os.path.basename(os.path.splitext(
                self.settings.args.get('ref'))[0])
            logger.info("Load reference structure file")
            try:
                cmd.load(self.settings.args.get('ref'), object=self._refname)
                self._reflag = True
            except CmdException:
                logger.error("Can't load %s", self._refname)
            method = self.settings.config.get('ali_method', 'align')
            mobile = self._pdbname
            target = self._refname
            try:
                align = cmd.keyword[method][0]
                logger.info("Align %s on %s" % (mobile, target))
                align(mobile, target, object="_aln")
                cmd.delete("_aln")
            except Exception as msg:
                logger.error('Wrong alignment method (%s)' % method)
                raise CmdException(msg)
        logger.info("Load constrains file %s" %
                    self.settings.args.get('constrains'))
        # read cnstr file
        self.load_cnstr()

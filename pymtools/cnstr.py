# coding=utf-8
"""
                            Pymol module reading NMR restraints
"""

from __future__ import print_function

import sys
import os
import logging
from pymol import cmd, CmdException
from .reader import CnstrFile


LOG = logging.getLogger(__name__)


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

    def _contact_label(self, eval_method, dist, lower, upper, idline):
        if eval_method == 'contact' and self._reflag:
            return "TP" if dist < float(self.settings.config.get(
                'contact', 5.0)) else "FP"
        elif eval_method == 'contact' and not self._reflag:
            return "unviol" if dist < float(self.settings.config.get(
                'contact', 5.0)) else "viol"
        elif eval_method == "bound" and self._reflag:
            if lower and upper:
                return "TP" if lower <= dist <= upper else "FP"
            else:
                LOG.error("No bound available for restraint no %d" %
                          (idline + 1))
        elif eval_method == "bound" and not self._reflag:
            if lower and upper:
                return "unviol" if lower <= dist <= upper else "viol"
            else:
                LOG.error("No bound available for restraint no %d" %
                          (idline + 1))
                return
        else:
            LOG.critical(
                "Evaluation method given is not supported (%s)"
                % eval_method)
            sys.exit(1)

    def load_cnstr(self):
        """
        Load constraint file
        """
        self._cnstrfile.load()

        names = []
        eval_method = self.settings.config.get('eval_method', 'contact')
        contact_type = self.settings.config.get('contact_type')
        groupname = self.settings.args.get('group')

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

        target, upper, lower, spec = None, None, None, None
        dist = []
        for idline in self._cnstrfile.lines:
            # For each restraint
            # print(idline + 1 == len(self._cnstrfile.lines))
            # If outfile and assignflag OR last line, save oldresults
            line = self._cnstrfile.lines[idline]

            # If assignflag, initialize all cnstr vars
            # Else if ambiflag, add vars to cnstr
            # Else return error
            resids = (line.get('resid1', ''), line.get('resid2', ''))
            atms = (line.get('atm1', ''), line.get('atm2', ''))
            ambiflag = line.get('ambiflag', False)
            assignflag = line.get('assignflag', False)
            if assignflag:
                dist = []
                spec = line.get('spec', 'noname')
                target = line.get('target', None)
                upper = line.get('upper', None)
                lower = line.get('lower', None)
                if target:
                    upper = float(target) + float(upper) if upper else None
                    lower = float(target) - float(lower) if lower else None

            if groupname and line.get(groupname) in ("False", "True"):
                group = groupname if line.get(groupname) == "True" \
                    else "un" + groupname
            else:
                group = line.get('group', None)

            extra = "\t".join((str(lower), str(target), str(upper))) if \
                eval_method == "bound" else ""
            # TODO: add a violation tolerance ?

            if '' in atms:
                LOG.warning('No atms related to cnsrt at line %s' %
                            (idline + 1))
            if '' in resids:
                LOG.error('No cnstr at line %s' % (idline + 1))
                continue

            # TODO: condition below can be simplified ...
            # Setting pymol selections
            if (contact_type == "all") or \
                    (contact_type == "min" and atms == ("CA", "CA")):
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
                LOG.info("Selecting restraint at line %s" % (idline + 1))
            else:
                LOG.debug("Ignoring line %s", idline + 1)
                LOG.debug("%s", line)
                continue
            # If ambig flag, append new dist to dist list
            # Otherwise initialize dist list with new dist
            if self._reflag:
                # If structure given, we look at reference distance
                dist = cmd.get_distance(sel3, sel4)
            else:
                dist = cmd.get_distance(sel1, sel2)

            # TODO: condition below can be simplified
            # Group can be decided by a contact treshold in the structure or in
            # the reference (if given). It can also be if distances in the
            # structure (or ref if given) lie inside bounds defined in the
            # cnstr file
            if not group:
                group = self._contact_label(eval_method, dist, lower, upper,
                                            idline)

            # PYMOL cmd.distance calls
            name = "_".join(('mobile', spec, group))
            nameref = "_".join(('ref', spec, group)) if self._reflag \
                else None
            if name not in names:
                names.append(name)
                if self._reflag:
                    names.append(nameref)
            try:
                LOG.debug("Distance cmd on %s and %s with groupname %s", sel1,
                          sel2, name)
                cmd.distance(name, sel1, sel2)
                if self._reflag:
                    LOG.debug("Distance cmd on %s and %s with groupname %s",
                              sel3,
                              sel4, nameref)
                    cmd.distance(nameref, sel3, sel4)
            except CmdException:
                LOG.error("Wrong selection")

            # OUTPUT line
            if outfile:
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
            cmd.set('stick_transparency', self.settings.config.get(
                'stick_alpha', 0.8))
            # cmd.hide('labels')
        cmd.color(self.settings.config.get('color', 'gray'), self._pdbname)
        if self._reflag:
            cmd.spectrum(selection=self._refname)

        if outfile:
            outfile.close()

    def run(self):
        """
        Load structure and restraint files and save pymol session with projected
        restraints.
        """
        LOG.info("Load structure file %s" % self.settings.args.get('pdb'))
        try:
            cmd.load(self.settings.args.get('pdb'), object=self._pdbname)
        except CmdException:
            LOG.error("Can't load %s", self.settings.args.get('pdb'))

        if self.settings.args.get('ref'):
            LOG.info("Load reference structure file")
            self._refname = os.path.basename(os.path.splitext(
                self.settings.args.get('ref'))[0])
            try:
                cmd.load(self.settings.args.get('ref'), object=self._refname)
                self._reflag = True
            except CmdException:
                LOG.error("Can't load %s", self._refname)
            method = self.settings.config.get('ali_method', 'align')
            mobile = self._pdbname
            target = self._refname
            try:
                align = cmd.keyword[method][0]
                LOG.info("Align %s on %s" % (mobile, target))
                align(mobile, target, object="_aln")
                cmd.delete("_aln")
            except Exception as msg:
                LOG.error('Wrong alignment method (%s)' % method)
                raise CmdException(msg)
        LOG.info("Load constrains file %s" %
                 self.settings.args.get('constrains'))
        # read cnstr file
        self.load_cnstr()

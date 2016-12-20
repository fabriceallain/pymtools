# coding=utf-8
"""
original Written by Jules Jacobsen (jacobsen@ebi.ac.uk). Feel free to do
whatever you like with this code.
extensively modified by Robert L. Campbell (rlc1@queensu.ca)
"""
import logging
from pymol import cmd, CmdException

LOG = logging.getLogger(__name__)


def align_all(target=None, mobile_selection='name ca',
              target_selection='name ca', cutoff=2, cycles=5, cgo_object=0,
              method='align'):
    """
  Aligns all models in a list to one target

  usage:
    align_all [target][target_selection=name ca][mobile_selection=name ca][cutoff=2][cycles=5][cgo_object=0][method='align']

        where method can be align, super or cealign

        where target specifies the model id you want to align all others against,
        and target_selection, mobile_selection, cutoff and cycles are options
        passed to the align or super command.

    Options for method='align' or method='super':

      By default the selection is all C-alpha atoms and the cutoff is 2 and the
      number of cycles is 5.
      Setting cgo_object to 1, will cause the generation of an alignment object for
      each object.  They will be named like <object>_on_<target>, where <object> and
      <target> will be replaced by the real object and target names.

    Example:
      align_all target=name1, mobile_selection=c. b & n. n+ca+c+o,target_selection=c. a & n. n+ca+c+o

  """
    cutoff = int(cutoff)
    cycles = int(cycles)
    cgo_object = int(cgo_object)

    object_list = cmd.get_names()
    object_list.remove(target)

    rmsd = {}
    rmsd_list = []
    for i in range(len(object_list)):
        if cgo_object:
            objectname = 'align_%s_on_%s' % (object_list[i], target)
            if method == 'align':
                rms = cmd.align('%s & %s' % (object_list[i], mobile_selection),
                                '%s & %s' % (target, target_selection),
                                cutoff=cutoff, cycles=cycles, object=objectname)
            elif method == 'super':
                rms = cmd.super('%s & %s' % (object_list[i], mobile_selection),
                                '%s & %s' % (target, target_selection),
                                cutoff=cutoff, cycles=cycles, object=objectname)
            elif method == 'cealign':
                rmsdict = cmd.cealign('%s & %s' % (target, target_selection),
                                      '%s & %s' % (
                                          object_list[i], mobile_selection))
                rms = [rmsdict['RMSD'], rmsdict['alignment_length'], 1, 0, 0]
            else:
                msg = "only 'align', 'super' and 'cealign' are accepted as " \
                      "methods"
                raise CmdException(msg)
        else:
            if method == 'align':
                rms = cmd.align('%s & %s' % (object_list[i], mobile_selection),
                                '%s & %s' % (target, target_selection),
                                cutoff=cutoff, cycles=cycles)
            elif method == 'super':
                rms = cmd.super('%s & %s' % (object_list[i], mobile_selection),
                                '%s & %s' % (target, target_selection),
                                cutoff=cutoff, cycles=cycles)
            elif method == 'cealign':
                rmsdict = cmd.cealign('%s & %s' % (target, target_selection),
                                      '%s & %s' % (
                                          object_list[i], mobile_selection))
                rms = [rmsdict['RMSD'], rmsdict['alignment_length'], 1, 0, 0]
            else:
                msg = "only 'align', 'super' and 'cealign' are accepted as " \
                      "methods"
                raise CmdException(msg)

        rmsd[object_list[i]] = (rms[0], rms[1])
        rmsd_list.append((object_list[i], rms[0], rms[1]))

    def mycmp(x, y):
        return (x > y) - (x < y)
    rmsd_list.sort(lambda x, y: mycmp(x[1], y[1]))
    # loop over dictionary and print out matrix of final rms values
    LOG.info("Aligning against: %s", target)
    for object_name in object_list:
        LOG.info("%s: %6.3f using %d atoms",
                 object_name, rmsd[object_name][0], rmsd[object_name][1])

    LOG.info("\nSorted from best match to worst:")
    for r in rmsd_list:
        LOG.info("%s: %6.3f using %d atoms", r[0], r[1], r[2])


cmd.extend('align_all', align_all)

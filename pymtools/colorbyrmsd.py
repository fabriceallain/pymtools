# coding: utf8
"""
http://pymolwiki.org/index.php/ColorByRMSD

Original Authors: Shivender Shandilya; Jason Vertrees
Complete rewrite by Thomas Holder

License: BSD-2-Clause
"""
from __future__ import print_function
import logging
# from pymol import cmd, CmdException
from tmalign import *


LOG = logging.getLogger(__name__)


def colorbyrmsd(mobile, target, doalign=1, dopretty=1, guide=1, method='super',
                quiet=1):
    """
DESCRIPTION

    Align two structures and show the structural deviations in color to more
    easily see variable regions.

    Colors each mobile/target atom-pair by distance (the name is a bit
    misleading).

    Modifies the B-factor columns in your original structures.

ARGUMENTS

    mobile = string: atom selection for mobile atoms

    target = string: atom selection for target atoms

    doalign = 0 or 1: Superpose selections before calculating distances
    {default: 1}

    dopretty = 0 or 1: Show nice representation and colors {default: 1}

EXAMPLE

    fetch 1ake 4ake, async=0
    remove chain B
    colorbyrmsd 1ake, 4ake
    """
    from chempy import cpv

    doalign, dopretty = int(doalign), int(dopretty)
    guide, quiet = int(guide), int(quiet)
    aln, seleboth = '_aln', 'objSelBoth'

    try:
        align = cmd.keyword[method][0]
    except Exception as msg:
        LOG.error('Wrong alignment method (%s)' % method)
        raise CmdException(msg)

    if guide:
        mobile = '(%s) and guide' % mobile
        target = '(%s) and guide' % target

    try:
        if doalign:
            # superpose
            LOG.info("Align %s with %s" % (mobile, target))
            align(mobile, target, quiet=0)

        # get alignment without superposing
        if method in ('super', 'align'):
            align(mobile, target, cycles=0, transform=0, object=aln)
        else:
            align(mobile, target, transform=0, object=aln)
    except Exception as msg:
        LOG.error(' Error: Alignment with method %s failed' % method)
        raise CmdException(msg)

    cmd.select(seleboth, '(%s) or (%s)' % (mobile, target))

    idx2coords = dict()
    cmd.iterate_state(-1, seleboth, 'idx2coords[model,index] = (x,y,z)',
                      space=locals())

    if cmd.count_atoms('?' + aln, state=1) == 0:
        # this should ensure that "aln" will be available as selectable object
        cmd.refresh()

    LOG.info("Saving RMS in b  factor columns")
    b_dict = dict()
    for col in cmd.get_raw_alignment(aln):
        assert len(col) == 2
        b = cpv.distance(idx2coords[col[0]], idx2coords[col[1]])
        for idx in col:
            b_dict[idx] = b

    cmd.alter(seleboth, 'b = b_dict.get((model, index), -1)', space=locals())

    if dopretty:
        cmd.orient(seleboth)
        cmd.show_as('cartoon', 'byobj ' + seleboth)
        cmd.color('gray', seleboth)
        cmd.spectrum('b', 'blue_red', seleboth + ' and b > -0.5')
        # cmd.color('forest', target)

    if not quiet:
        LOG.info(" ColorByRMSD: Minimum Distance: %.2f" % (
            min(b_dict.values())))
        LOG.info(" ColorByRMSD: Maximum Distance: %.2f" % (
            max(b_dict.values())))
        LOG.info(" ColorByRMSD: Average Distance: %.2f" % (
            sum(b_dict.values()) / len(b_dict)))

    # cmd.delete(aln)
    # cmd.delete(seleboth)

cmd.extend('colorbyrmsd', colorbyrmsd)

# tab-completion of arguments
cmd.auto_arg[0]['colorbyrmsd'] = cmd.auto_arg[0]['align']
cmd.auto_arg[1]['colorbyrmsd'] = cmd.auto_arg[1]['align']

# vi: ts=4:sw=4:smarttab:expandtab

"""
Pymol module adapted from nmr_cnstr
See more here: http://www.pymolwiki.org/index.php/nmr_cnstr
################################################################################
# Pymol Script: For visualizing the NMR constrains (DYANA & CNS format),       #
#   on top of the calculated structure.                                        #
#               Author: Evangelos Papadopoulos.                                #
#  previous affiliation: Dept. of Biochemistry and Biophysics,                 #
#                       Arrhenius Laboratories,                                #
#                       Stockholm University                                   #
#                       SE-106 91 Stockholm, Sweden                            #
#                email:evangelos.papadopoulos@gmail.com                        #
#                NOTES: This is a preliminary version.                         #
#                                                                              #
#     Reference: please if you find this script usefull add the following      #
#     reference:                                                               #
#     NMR Solution Structure of the Peptide Fragment 1-30, Derived from        #
#     Unprocessed Mouse Doppel                                                 #
#     Protein, in DHPC Micelles. Biochemistry. 2006 Jan 10;45(1):159-166.      #
#     PMID: 16388591                                                           #
#                                                                              #
################################################################################
"""

from __future__ import print_function
from pymol import cmd


def upl(fname):

    f = open(fname, 'r')
    i = 1
    upline = f.readline()

    while upline != '':

        print(upline, i)
        cnsline = upline.split()
        cmd.dist('upl' + str(i), 'i. ' + cnsline[0] + ' & n. ' + cnsline[2],
                 'i. ' + cnsline[3] + ' & n. ' + cnsline[5])
        upline = f.readline()
        i += 1
#
    f.close()
    cmd.hide('labels')
    cmd.set('dash_gap', 0.05)
    cmd.do("orient")
    cmd.set('movie_delay', 1500)


def cns(fname):

    f = open(fname, 'r')
    i = 1
    upline = f.readline()
    print(upline, i)
    while upline != '':
        if upline == '\n':
            upline = f.readline()
            continue
        cnsline = upline.split()
        print(cnsline, i)
        if cnsline[0] == 'assign':
            print('CNS')
            if cnsline[5] == 'HB*':
                print('CNS***')
            cmd.dist('upl' + str(i), 'i. ' + cnsline[2] + ' & n. ' +
                     cnsline[5], 'i. ' + cnsline[7] + ' & n. ' + cnsline[10])
        i += 1
        upline = f.readline()
        print('*' + upline + '*', i)

    f.close()
    cmd.set('dash_gap', 0.05)
    cmd.hide('labels')
    cmd.do("orient")
    cmd.set('movie_delay', 1500)

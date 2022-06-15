#!/usr/bin/env python
#
# This is a self-test of the pandeia-repo BIT translation code.
#
# Assume:  the pyetc test repo contains true/correct translations
#          from apt queries to peng files under: pyetc/test/engine/bit
#          covering all nightly tests of these queries.  This may not
#          be all the queries but it (I believe) is the full unique set.
#
# Assume:  each peng file in the pyetc/test/engine/bit tree includes the
#          apt query used to generate it as a comment at the top
#
# The Plan:  find all of these files and for every one, read the apt
#            query, use pandeia repo code to generate the peng, verify
#            the the output peng is the same as what was in pyetc's repo
#
# How to use: Change the repo_dir variable in main function to the pyetc
#             repo in your local machine. 
#             Run with `python bit_test_query_to_peng.py`.

import os, sys
from pandeia.engine.helpers.bit import bit
from pandeia.engine.helpers.bit.pyetc_util import write_dict, fast_write_dict

from pandeia_test import utils


def create_peng_ref_obj(fname):
    """ read in a single peng and get all of its meta info to make a dict """
    retval = {}
    f = open(fname, 'r')
    retval['buf'] = f.read().strip() # no need for final newline
    f.close()

    pathexp = 'pyetc/test/engine/bit/'
    idx = fname.index(pathexp) # let throw on err
    retval['fname'] = fname # e.g. '/Users/bob/dev/pyetc/test/engine/bit/acs/imag/sbc/sbc/f115lp/f5v_mag15v.peng'
    retval['testname'] = fname[idx+len(pathexp):] # e.g. 'acs/imag/sbc/sbc/f115lp/f5v_mag15v.peng'

    lines = retval['buf'].split('\n')
    assert len(lines) > 5, 'Unexpected peng format for: '+fname
    # header is like:
    # 'url query: type=bot&magnitude=15&config=ACS/SBC&filter0=F115LP&aperture=SBC&requestCount=2&request0=BrightestPixelRate&request1=ImageRate&stellarType=F5V&magType, test name: acs/imag/sbc/sbc/f115lp/f5v_mag15v.peng'
    header = lines[0]
    assert header[0] == '#', 'Unexpected 1st line of peng: '+lines[0]
    parts = header.split(',')
    query = parts[0].split(' ')[-1]
    retval['query'] = query

    # final sanity check
    tname = parts[-1].split(' ')[-1]
    tname = tname.strip("'")
    assert tname == retval['testname'], f'test name from fname: {retval["testname"]} != test name from header: {tname}'

    return retval


def find_all_pyetc_bit_pengs(pyetcrepo_test_engine_bit_dir):
    """ Return a list of dict-objects, one for each apt query found at the top comment of
    a peng file in the pyetc/test/engine/bit repo. """

    # check input
    for inst in ('acs','cos','stis','wfc3ir','wfc3uvis'):
        assert os.path.exists(pyetcrepo_test_engine_bit_dir+'/'+inst), 'Cannot find: '+pyetcrepo_test_engine_bit_dir+'/'+inst

    # all engine/bit pengs from pyetc/test
    pengs = utils.rglob(pyetcrepo_test_engine_bit_dir,'*.peng')
    retval = [create_peng_ref_obj(p) for p in pengs]
    return retval


def main():

    # First find all of our "truth" files.
    #
    # list[dict{fname: str, query: str, buf: str}]
    repo_dir = '/Users/sontag/dev/pyetc/test/engine/bit'
    bit_peng_files = find_all_pyetc_bit_pengs(repo_dir)

    # loop over them
    for peng_obj in bit_peng_files:

        print('')
        # use pandeia code to convert to peng
        orig_file = peng_obj['fname']
        query = peng_obj['query']
        testname = peng_obj['testname']
        pweb_dict, peng_dict = bit.query_to_peng(query)

        # pweb currently gets dropped on floor
        new_file = orig_file+'_regen'

        # write new file
        if os.path.exists(new_file):
            os.remove(new_file)
        fast_write_dict(peng_dict, new_file, comment=f"'url query: {query}, test name: {testname}'")
        print('Orig: ', orig_file)
        print('Wrote:', new_file)

        # diff it for pass/fail but don't capture (for now)
        if sys.platform=='darwin':
            diffcmd = "/opt/sw/bin/diff"
        else:
            diffcmd = "/usr/bin/diff"
        assert os.path.exists(diffcmd), 'Cannot find: '+diffcmd
        # !!! expect MJD diffs for now !!!
        diffcmd += " --ignore-matching-lines='obsmode.*,MJD.59' "+orig_file+" "+new_file

        x = os.system(diffcmd)
        assert x == 0, 'New file is different: '+new_file
        os.remove(new_file)

    print(f'\nJust looped over {len(bit_peng_files)} in {repo_dir}')
    print('All tests passed')


#
# main routine
#
if __name__=='__main__': # in case something else imports this file

    main()

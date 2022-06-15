#!/usr/bin/env python
#
# This is a regen tool to go from APT query to all files that can
# come of it:  .pweb, .peng, .jeng
#
# Assume:  the pyetc test repo contains the latest correct APT query sets
#

import os, re, sys
from pandeia.engine.helpers.bit import bit
from pandeia.engine.helpers.bit.pyetc_util import write_dict, fast_write_dict
from pandeia.engine.helpers.peng import peng3jeng


def get_all_queries(repo_dir_to_apt_queries):
    retval = {}
    for inst in ('acs', 'cos', 'stis', 'wfc3'):
        queryfname = repo_dir_to_apt_queries+'/'+inst+'.txt'
        assert os.path.exists(queryfname), 'Could not find: '+queryfname
        with open(queryfname, 'r') as f:
            iqueries = [q.strip() for q in f.readlines() if len(q.strip()) > 0]
            retval[inst] = iqueries
    return retval


any_unsafe = re.compile('[^a-zA-Z0-9_=\.\-]')
def filenamify(qstr):
    """
    Convert some chars to things that work well in file/dirnames, but do so
    in a way that can be fully back-convertable in case we want to convert back.
    Note if the conversion back is done, it muct be done in the reverse order!
    """

    fixed = qstr

    assert '_d_' not in fixed, 'Cannot use "_d_" as a replacement for "." in: '+qstr
    fixed = fixed.replace('.','_d_')

    assert '.' not in fixed, 'Cannot use "." as a replacement for "&" in: '+qstr
    fixed = fixed.replace('&','.')

    assert '_s_' not in fixed, 'Cannot use "_s_" as a replacement for "/" in: '+qstr
    fixed = fixed.replace('/','_s_')

    # make sure the rest is just alpha and numeric and dash/equal/under/dot
    assert any_unsafe.search(fixed) == None, 'Need to handle this query: '+qstr
    return fixed


def write_out_files_for_queries(queries_dict, out_dir, overwrite=True):
    """ Fo each query, dump all the files """
    # logic to show output but reduce debugging
    last_shown = ''
    # go by instrument
    insts = sorted(queries_dict.keys())
    for inst in insts:
        # then by each query
        for query in queries_dict[inst]:
            # output
            line = f'{inst}: {query}'
            line = line[0:line.find('requestCount')-1]+' ...'
            if line != last_shown:
                print(line)
                last_shown = line
            # dirs and fnames
            basedir = out_dir+'/'+inst
            if not os.path.isdir(basedir):
                os.mkdir(basedir)
            basenm = basedir+'/'+filenamify(query)
            out_pweb = basenm+'.pweb'
            out_peng = basenm+'.peng'
            out_jeng = basenm+'.jeng'
            for fname in (out_pweb, out_peng, out_jeng):
                if os.path.exists(fname):
                    if overwrite:
                        os.remove(fname)
                    else:
                        raise IOError('File exists: '+fname)
            # pweb & peng: do the conversions, write em
            pweb_dict, peng_dict = bit.query_to_peng(query)
            fast_write_dict(pweb_dict, out_pweb, comment=f'apt query: {query}')
            fast_write_dict(peng_dict, out_peng, comment=f'apt query: {query}')
            # jeng: convert and write
            if inst == 'stis' and 'centralWavelength' not in query: # some dont work yet
                jeng_dict = peng3jeng.pyetc_in_dict_to_pandeia_in_dict(peng_dict, None, pweb_dict=pweb_dict)
                fast_write_dict(jeng_dict, out_jeng, comment=f'apt query: {query}')


def main():

    repo_dir_to_apt_queries = '/Users/sontag/dev/pyetc/pyetc/etc_web/bit/apt_queries'
    out_dir = '/Users/sontag/dev/bit_pandeia'

    qdict = get_all_queries(repo_dir_to_apt_queries)
    write_out_files_for_queries(qdict, out_dir)

    print('Done!')


#
# main routine
#
if __name__=='__main__': # in case something else imports this file

    main()

#! /usr/bin/env python
"""
This is a collection of pyetc utils and stubs to get the
pyetc code working (minimally) under pandeia as needed.
"""

import os, pprint, re, six, sys
import numpy as np
from numpy import nan, inf, NaN, Inf
from collections import OrderedDict
import ast

class mock_log:
    """ Mock up a "log" pkg that simply matches what the old pyetc code needs.
        The intent is to NOT bring in any more pyetc code than necessary. """

    def __init__(self, be_verbose=False):
        self.VERBOSE_FLAG = be_verbose

    def info(self, str1, str2='', str3=''):
        print(str1, str2, str3)

    def error(self, str1, str2='', str3='', str4=''):
        print(str1, str2, str3, str4, file=sys.stderr)

    def write(self, str1='', eol='\n'):
        print(str1, end=eol, flush=True)

    verbose = info # oddly this is the way it works in pyetc; caller checks VERBOSE_FLAG
    warning = error # these seem to be the same basically

log = mock_log(False)

# rather than take PYSYN_CDBS from a pyetc config system, pull from env, but don't make a critical error on import
PYSYN_CDBS = os.environ.get('PYSYN_CDBS','PYSYN_CDBS_is_not_set' )

# pyetc's dynamic_import
PACKAGE_RE = re.compile("[A-Za-z_0-9.]+")
def dynamic_import(package):
    """imports a module specified by string `package` which is
    not known until runtime.

    Returns a module/package.

    The purpose of this function is to concentrate a number of uses
    of the "exec" statement in once place where it can be policed
    for security reasons.
    """
    if not PACKAGE_RE.match(package):
        raise ImportError("Invalid dynamic import " + repr(package))

    # six.exec_ doesn't support keyword arguments, even though the
    # documentation says it does.
    six.exec_("import " + package + " as pkg_name", None, locals())

    m = locals()['pkg_name']
    return m


def read_dict(fname):
    ''' read a python dictionary from a file that was written with write_dict. '''
    f=open(fname,'r', encoding="utf-8")
    datastr = f.read()
    f.close()
    # convert DOS file to Unix - otherwise the eval will fail
    datastr = datastr.replace('\r','')
    try :
        datadict = safe_eval(datastr)
    except Exception as e:
        print('EXCEPTION:',e)
        print('cannot eval data in file ',fname)
        raise
    return datadict


def write_dict(datadict, fname, header=None):
    '''write a python dictionary to a file in such a way that it
    can be read with read_dict.
    At present, this means using pprint, but encapsulating it
    here will make it easier if we decide to change the format as the
    project evolves.
    A string, or a list of strings, passed in the optional header keyword
    will be prepended with a # sign before being written to the file,
    which read_dict will ignore. This allows the files to be
    documented.
    '''

    fh=open(fname,'w')
    #First write commented-out header
    if header is not None:
        #Support either a list or a string
        if isinstance(header,list):
            for line in header:
                fh.write("# %s\n"%line)
        else:
            fh.write("# %s\n"%header)
    #Now write the data
    pprint.pprint(datadict,fh)
    fh.close()


# stolen/modified from spidring/metadict.py (see code there as to why this is this way)
def fast_write_dict(adict, fname, comment):
    f = open(fname, 'w')
    f.write("# "+ comment + "\n")
    f.write("{\n")
    for key,val in sorted(adict.items()):
        f.write("\t" + repr(key) + " : " + repr(val) + ",\n")
    f.write("}\n")
    f.close()


# pyetc took from Python-2.7 and enhanced to support np.nan and np.inf
def safe_eval(node_or_string):
    """
    Safely evaluate an expression node or a string containing a Python
    expression.  The string or node provided may only consist of the following
    Python literal structures: strings, numbers, tuples, lists, dicts, booleans,
    and None.
    """
    _safe_names = {'None': None, 'True': True, 'False': False,
                   'nan':nan, 'NaN':np.NaN,
                   'inf':inf, 'Inf':np.Inf}
    if isinstance(node_or_string, six.string_types):
        node_or_string = ast.parse(node_or_string, mode='eval')
    if isinstance(node_or_string, ast.Expression):
        node_or_string = node_or_string.body
    def _convert(node):
        # Old code used to return instances of ast.Dict.
        # New code sometimes retuns dicts disguised as lists
        # of tuples (actuaally, ast.List and ast.Tuple), as
        # the first element of the 'args' attribute of an
        # ast.Call instace.
        if isinstance(node, ast.Call):
            result = OrderedDict()
            if len(node.args) > 0:
                for element in node.args[0].elts:
                    name  = element.elts[0].s
                    v = element.elts[1]
                    if isinstance(v, ast.Str):
                        value = v.s
                    elif isinstance(v, ast.Num):
                        value = v.n
                    else:
                        value = None
                    result[name] = value
            return result
        elif isinstance(node, ast.Str):
            return node.s
        elif isinstance(node, ast.Num):
            return node.n
        elif isinstance(node, ast.Tuple):
            return tuple(map(_convert, node.elts))
        elif isinstance(node, ast.List):
            return list(map(_convert, node.elts))
        elif isinstance(node, ast.Dict):
            return OrderedDict((_convert(k), _convert(v)) for k, v
                         in zip(node.keys, node.values))
        elif isinstance(node, ast.Name):
            if node.id in _safe_names:
                return _safe_names[node.id]
        elif isinstance(node, ast.BinOp) and \
             isinstance(node.op, (ast.Add, ast.Sub)) and \
             isinstance(node.right, ast.Num) and \
             isinstance(node.right.n, complex) and \
             isinstance(node.left, ast.Num) and \
             isinstance(node.left.n, six.integer_types) and \
             isinstance(node.left.n,  float):
            left = node.left.n
            right = node.right.n
            if isinstance(node.op, ast.Add):
                return left + right
            else:
                return left - right

        # when the expression contains a negative number,
        # under python 3 the negative number will be cast
        # as an instance of UnaryOp with first operand of
        # type USub.
        elif isinstance(node, ast.NameConstant):
            if sys.version_info[0] >= 3:
                return node.value
        elif isinstance(node, ast.UnaryOp) and \
             isinstance(node.op, ast.USub) and \
             isinstance(node.operand, ast.Num):
            return -(node.operand.n)

        raise ValueError('malformed string')

    return _convert(node_or_string)


def stable_dict(dictionary=None):
    ''' Returns a dictionary of the appropriate type.

        Under python 2, dictionary iterators are stable
        in between successivea calls within the same process,
        or from process to procees. Under python 3, dictionary
        iterators return elements in random order, different
        from call to call. This function returns an OrderedDict
        instance that works under both 2 and 3, but leaves open
        the possibility to use other kinds of dictionary if the
        need ever arises.

        This function should be used to create dictionary instances
        throughout the pyetc code, instead of direct calls to dict().
    '''
    if dictionary:
        return OrderedDict([item for item in six.iteritems(dictionary)])
    else:
        return OrderedDict()

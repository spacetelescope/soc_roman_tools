import ast
import io
import os
import six
import numpy as np
from numpy import nan, inf, NaN, Inf
from collections import defaultdict, OrderedDict

def read_dict(fname):
    '''read a python dictionary from a file that was written
    with write_dict.
    '''
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

# Taken from Python-2.7 and enhanced to support np.nan and np.inf
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

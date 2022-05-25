"""
Module for setting up logging behavior in other modules.
"""

import logging
import time


def configure_logging(level=logging.INFO):
    """
    Purpose
    -------
    Set up logging messages.

    Inputs
    ------
    level (integer):
        Minimum logging level to display messages. These are technically
        integers, but can use inputs like `logging.INFO` or `logging.DEBUG`.

    Returns
    -------
    None
    """

    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)

    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S %p',
                        level=level)


def timeit(func):
    """
    Purpose
    -------
    A decorator to compute the runtime of a function and include that
    information in the logging.info messages.
    """

    def wrapped(*args, **kwargs):

        # Compute times before function runs
        t1_cpu = time.process_time()
        t1_user = time.time()

        # Run function
        func(*args, **kwargs)

        # Compute times after function completes.
        t2_cpu = time.process_time()
        t2_user = time.time()

        logging.info(f'Elapsed CPU time: {t2_cpu - t1_cpu} seconds')
        logging.info(f'Elapsed real time: {t2_user - t1_user} seconds')

    return wrapped

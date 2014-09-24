"""
vic.py
"""

from __future__ import print_function
import os
import tempfile
import subprocess
import pandas as pd


# -------------------------------------------------------------------- #
class VIC_Error(RuntimeError):
    pass
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
class VIC(object):

    def __init__(self, executable):
        if os.path.isfile(executable) and os.access(executable, os.X_OK):
            self.executable = executable
            self.version = self._get_version()
            self.options = self._get_options()
        else:
            raise VIC_Error('%s is not a valid executable' % executable)

    def _get_version(self):
        """Get the version of VIC from the executable"""
        return self._call_vic('-v')[1]

    def _get_options(self):
        """Get the compile time options of VIC from the executable"""
        return self._call_vic('-o')[1]

    def run(self, global_param):
        """run VIC"""

        if os.path.isfile(global_param):
            global_param_file = global_param
        else:
            # global_param is a string
            f, global_param_file = tempfile.mkstemp(prefix='vic.global.param.',
                                                    suffix='.txt',
                                                    text=True)
            with open(global_param_file, mode='w') as f:
                f.write(global_param)

        stdout = self._call_vic('-g', global_param_file)

        return stdout

    def _call_vic(self, *args):
        vic_args = [self.executable]+[a for a in args]

        proc = subprocess.Popen(' '.join(vic_args),
                                shell=True,
                                stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)
        retvals = proc.communicate()

        stdout = retvals[0]
        stderr = retvals[1]
        returncode = proc.returncode

        if returncode > 0:
            raise VIC_Error(stderr)

        return stdout, stderr
# -------------------------------------------------------------------- #


def read_vic_ascii(filepath, header=True, parse_dates=True,
                   datetime_index=None, names=None, **kwargs):
    """Generic reader function for VIC ASCII output with a standard header
    filepath: path to VIC output file
    header (True or False):  Standard VIC header is present
    parse_dates (True or False): Parse dates from file
    datetime_index (Pandas.tseries.index.DatetimeIndex):  Index to use as datetime index
    names (list like): variable names
    **kwargs: passed to Pandas.read_table

    returns Pandas.DataFrame
    """
    kwargs['header'] = None

    if header is True:
        kwargs['skiprows'] = 6

        # get names
        if names is None:
            with open(filepath) as f:
                # skip lines 0 through 5
                for _ in range(5):
                    next(f)

                # process header
                names = next(f)
                names = names.strip('#').replace('OUT_', '').split()

    kwargs['names'] = names

    if parse_dates:
        time_cols = ['YEAR', 'MONTH', 'DAY']
        if 'HOUR' in names:
            time_cols.append('HOUR')
        kwargs['parse_dates'] = {'datetime': time_cols}
        kwargs['index_col'] = 0

    df = pd.read_table(filepath, **kwargs)

    if datetime_index is not None:
        df.index = datetime_index

    return df

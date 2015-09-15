# stdlib imports
import os
import atexit

# pytables imports
import tables

# rootpy imports
from rootpy.io import root_open, TemporaryFile

# higgstautau imports
from higgstautau import datasets

# local imports
from .. import NTUPLE_PATH, DEFAULT_STUDENT
from . import log; log = log[__name__]


DB = datasets.Database(name='datasets_hh', verbose=False)
FILES = {}
TEMPFILE = TemporaryFile()


def get_file(ntuple_path=NTUPLE_PATH, student=DEFAULT_STUDENT, hdf=False, suffix='', force_reopen=False):
    ext = '.h5' if hdf else '.root'
    filename = student + ext
    if filename in FILES and not force_reopen:
        return FILES[filename]
    file_path = os.path.join(ntuple_path, student + suffix, filename)
    log.info("opening {0} ...".format(file_path))
    if hdf:
        student_file = tables.open_file(file_path)#, driver="H5FD_CORE")
    else:
        student_file = root_open(file_path, 'READ')
    FILES[filename] = student_file
    log.info("got files")
    return student_file


@atexit.register
def cleanup():
    if TEMPFILE:
        TEMPFILE.close()
    for filehandle in FILES.values():
        if filehandle:
            filehandle.close()

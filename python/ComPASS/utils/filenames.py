import os
from mpi4py import MPI

default_case_name = 'compass'

def output_directory(case_name=None):
    tmp, tail = os.path.split(__file__)
    tmp, tail = os.path.split(tmp)
    assert tail=='utils'
    tmp, tail = os.path.split(tmp)
    assert tail=='ComPASS'
    compass_python_directory, tail = os.path.split(tmp)
    if case_name is None:
        case_name = default_case_name
    output = os.path.join(compass_python_directory, 'tests',
                                'output-' + os.path.splitext(os.path.basename(case_name))[0])
    output = os.path.abspath(output)
    # master proc manages directory creation
    comm = MPI.COMM_WORLD
    if comm.rank==0:
        if not os.path.exists(output):
          os.makedirs(output)
    comm.Barrier() # wait for every process to synchronize
    assert os.path.exists(output)
    return output 

def output_directory_and_logfile(case_name=None):
    output = output_directory(case_name)
    if case_name is None:
        case_name = default_case_name
    case_name = os.path.splitext(os.path.basename(case_name))[0]
    logfile = os.path.join(output, case_name + '.log')
    return output, logfile

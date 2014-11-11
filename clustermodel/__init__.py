#
from __future__ import print_function
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

__version__ = "0.13"
from . import simulate
from . import clustermodel
from . import feature
from . import plotting

from clustermodel import clustered_model
from feature import feature_gen, ClusterFeature, cluster_to_dataframe

from multiprocessing import cpu_count
CPUS = min(cpu_count(), 12)

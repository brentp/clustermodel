#
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

__version__ = "0.1"
from . import simulate
from . import clustermodel
from . import feature
from . import plotting

from clustermodel import clustered_model
from feature import feature_gen, ClusterFeature, cluster_to_dataframe

from multiprocessing import cpu_count
CPUS = max(cpu_count(), 12)

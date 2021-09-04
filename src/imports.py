'''Import modules for the project'''

######## user imports ########
from utils import *
from genimport import *
from common import *
##### standard libraries #####
import numpy as np
import scipy
from scipy import linalg as la
import scipy.sparse as sp
from math import floor, ceil
import glob
import time
import zipfile as zp
import gc
from itertools import product
from scipy.special import j0, y0
from scipy.special import hankel1
from math import log10, floor, ceil
from tqdm import tqdm # progress bar
from memory_profiler import memory_usage
import _pickle as cPickle # C language pickle faster than pickle
##### plotting libraries #####
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.lines as mlines
import matplotlib.ticker as mtick
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

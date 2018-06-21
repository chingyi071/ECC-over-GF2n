from BCJR import BCJR
from util import *

p_mat = read_mat('par.csv')
bcjr = BCJR(p_mat)
bcjr.plot_sections(0,5)
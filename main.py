from BCJR import BCJR
from util import *

p_mat = read_mat('file3.csv')
bcjr = BCJR(p_mat, b=2)
bcjr.plot_sections(0,5)

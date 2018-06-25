from BCJR import BCJR
from util import *

p_mat = read_mat('file3.csv')
bcjr1 = BCJR(p_mat, b=2)
bcjr1.plot_sections(0,5)

p_mat = read_mat('file2.csv')
bcjr1 = BCJR(p_mat, b=1)
bcjr1.plot_sections(0,7)

p_mat = read_mat('file1.csv')
bcjr1 = BCJR(p_mat, b=1)
bcjr1.plot_sections(0,5)

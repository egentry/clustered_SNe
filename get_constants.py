from astropy import constants as const 
from astropy import units as u 


print("G = ", str("{0:e}").format(const.G.cgs.value), sep='\t')
print("k_B = ", str("{0:e}").format(const.k_B.cgs.value), sep='\t')
print("m_p = ", str("{0:e}").format(const.m_p.cgs.value), sep='\t')
print("au = ", str("{0:e}").format(const.au.cgs.value), sep='\t')
print("pc = ", str("{0:e}").format(const.pc.cgs.value), sep='\t')
print("M_sun = ", str("{0:e}").format(const.M_sun.cgs.value), sep='\t')
print("R_sun = ", str("{0:e}").format(const.R_sun.cgs.value), sep='\t')
print("yr = ", str("{0:e}").format(u.year.to('s')), sep='\t')

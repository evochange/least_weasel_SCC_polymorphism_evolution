import dadi
import numpy as np

#Define best dadi model:
def model_func(params, ns, pts):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu_1 + (nu11 - nu_1) * (t / t1)
	nu2_func = lambda t: nu_2 * (nu12 / nu_2) ** (t / t1)
	phi = dadi.Integration.two_pops(phi, xx, T=t1, nu1=nu1_func, nu2=nu2_func, m12=m1_12, m21=m1_21)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

#Original data FS:
data = dadi.Spectrum.from_file('./ANC13_NorthCentral_snps.filt3.25K.proj25perc.folded.fs')
pts = [70, 80, 90]
ns = data.sample_sizes

#Best point estimates:
p0 = [2.652037053003643, 0.01, 0.7267855677582798, 0.5307909927243268, 2.924037654441933, 0.44194312627480203, 2.450550638982358]
#lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0]
#upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0]

#Calculate Best model FS:
func_ex = dadi.Numerics.make_extrap_log_func(model_func)

model_fs = func_ex(p0, ns, pts)


#Plot model vs empirical FS comparison:
import matplotlib.pyplot as plt

fig = plt.figure(219033)

fig.clear()

dadi.Plotting.plot_2d_comp_multinom(model_fs, data, vmin=1e-2)

fig.savefig('Plot_best_AIC_ANC13_2D_comparison_multinom.svg')
fig.savefig('Plot_best_AIC_ANC13_2D_comparison_multinom.pdf')
fig.savefig('Plot_best_AIC_ANC13_2D_comparison_multinom.png')


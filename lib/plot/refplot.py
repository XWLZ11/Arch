import matplotlib.pyplot as plt
import numpy as np
import os

path_dir = '/home/tsingularity/Desktop/tianff/202203_MDCarbonicAcid/literature/ref/'

dict_ref = {
        'Aminov et al. 2019': '2019_PNAS_DanielAminov/Fig_1_kelvin.csv',
        'Wang et al. 2010': '2010_JPCA_WangXiaoguang/Sfig_3_kelvin.csv',
        'Pines et al. 2016': '2016_JPCB_DinePines/pka_kelvin.csv',
        'Adamczyk et al. 2009': '2009_Science_KatrinAdamczyk/pka_kelvin.csv',
        }

refnames, dict_data = [], {}
for refname, path_ref in dict_ref.items():
    refnames.append(refname)
    path_file = os.path.join(path_dir, path_ref)
    dict_data[refname] = np.loadtxt(path_file)


plt.style.use('my_scripts/Arch/lib/plot/mplstyle/test1.rc')
fig, ax = plt.subplots()  

ax.errorbar(dict_data[refnames[0]][:, 0], dict_data[refnames[0]][:, 1], 
            yerr=dict_data[refnames[0]][:, 2], fmt=':o', color='#4169E1', 
            ecolor='#4169E1', elinewidth=2.0, capsize=4.0, ms=6.0, 
            label=refnames[0])
ax.errorbar(dict_data[refnames[1]][:, 0], dict_data[refnames[1]][:, 1], 
            yerr=dict_data[refnames[1]][:, 2], fmt=':>', color='#32CD32', 
            ecolor='#32CD32', elinewidth=2.0, capsize=4.0, ms=7.0, 
            label=refnames[1])
ax.errorbar(dict_data[refnames[2]][0], dict_data[refnames[2]][1], 
            yerr=dict_data[refnames[2]][2], fmt=':v', color='#FF8C00', 
            ecolor='#FF8C00', elinewidth=2.0, capsize=4.0, ms=7.0, 
            label=refnames[2])
ax.errorbar(dict_data[refnames[3]][0], dict_data[refnames[3]][1], 
            yerr=dict_data[refnames[3]][2], fmt=':<', color='#BB3F3F', 
            ecolor='#BB3F3F', elinewidth=2.0, capsize=4.0, ms=7.0, 
            label=refnames[3])  

ax.set_xlabel("Temperature (K)")
ax.set_ylabel("pKa")                      
ax.legend(frameon=False)
#plt.savefig("/home/tsingularity/Desktop/test.png")

plt.show()

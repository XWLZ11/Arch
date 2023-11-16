__all__ = ['compute_RDF']

import MDAnalysis as mda


from ..check import check_module, check_func

def main():
    dirpath = '/home/tsingularity/Desktop/scripts_test/target/LAMMPS_lammpstrj/'
    dumpname = 'cp.lammpstrj'

    u = mda.Universe(dirpath+dumpname, format='LAMMPSDUMP')
    s1 = u.select_atoms('type 1')
    s2 = u.select_atoms('type 2')
    [r, g_r] = compute_RDF(u, s1, s2, compute_coord=True)
    coord_number = compute_coord(r, g_r, index_r1=300)
    plot_RDF(r, g_r)

def compute_RDF(dir_init="", dir_target="", file_lammpstrj='dump.lammpstrj', 
                selection1='type 1', selection2='type 2', start=0, stop=-1, step=1):
    os = check_module('os')
    rdf = check_module('MDAnalysis.analysis.rdf')

    u = mda.Universe(os.path.join(dir_init, file_lammpstrj), format='LAMMPSDUMP')
    print()
    s1 = u.select_atoms(selection1)
    s2 = u.select_atoms(selection2)
    value_rdf = rdf.InterRDF(s1, s2, range=(0.0, 15.0), #exclusion_block=(1,7), 
                            nbins=300, norm='rdf', verbose=True)
    value_rdf.run(start, stop, step)    

    return [value_rdf.results.bins, value_rdf.results.rdf]

def compute_RDF1(u, selection1, selection2, compute_coord=False):
    import MDAnalysis.analysis.rdf as rdf
    value_rdf = rdf.InterRDF(selection1, selection2, range=(0.0, 15.0), exclusion_block=(1,7), 
                            nbins=300, norm='rdf', verbose=True)
    value_rdf.run(start=0, stop=-1, step=1)

    return [value_rdf.results.bins, value_rdf.results.rdf]

def compute_coord(r, g_r, index_r1):
    from scipy.integrate import trapezoid
    index_r1 = int(index_r1*300/1500)
    coord_number = trapezoid(g_r[:index_r1]*r[:index_r1]**2, r[:index_r1], axis=0)
    print("Coordination Numbers is %4.2f"%(coord_number))
    return coord_number

def plot_RDF(r, g_r):
    import matplotlib.pyplot as plt
    plt.plot(r, g_r, label=r'$g_{CuNi}$(r)')
    plt.ylabel('g(r)')
    plt.xlim(0, 6)
    plt.legend()
    plt.xlabel(r'Radius/$\AA$')
    plt.title(r'The RDF of CuNi')
    plt.show()

if __name__ == '__main__':
    import warnings
    warnings.simplefilter("ignore", category=UserWarning)
    main()

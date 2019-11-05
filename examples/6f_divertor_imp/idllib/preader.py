import numpy as np
from boututils import save2nc

__version__ = '0.1.0'
__date__ = '11152017'
__author__ = 'J.G. Chen'
__email__ = 'cjgls@pku.edu.cn'

__all__ = ['preader']

def preader(filename):
    """ profile reader for output from CORSICA integrated modeling code.

    Patameters
    ----------
    filename : str
        File output from CORSICA integrated modeling code.

    Returns
    -------
    tuple(data, comments)
    data : dict
        profiles' name and data
    comments : dict
        vars' descriptions

    """

    with open(filename, 'r') as f:
       fl = f.read().strip().splitlines()

    data = {}
    comments = {}
    comments['description'] = fl[0:5]
    size = len(fl)
    ind = 7

    while ind < size:
        if ':' in fl[ind]:
            # new var
            line = fl[ind].split(':')
            var = line[0].strip().replace(' ', '_')
            comments[var] = line[1].strip()
            tmp_data = []
            ind += 1
        elif fl[ind] == 'Additional 0D data':
            ind += 1
            while True:
                try:
                    line = fl[ind].replace(' ', '').split('=')
                    data[line[0]] = float(line[1].replace('D', 'E'))
                    ind += 1
                except IndexError:
                    break
        else:
            # read data line by line
            while ':' not in fl[ind]:
                tmp_data.extend(
                    map(float, fl[ind].replace('D', 'E').split()))
                ind += 1
                if fl[ind] == ' ':
                    ind += 1
                    break

            if len(tmp_data) == 1:
                data[var] = tmp_data[0]
            else:
                data[var] = np.asarray(tmp_data)
    return data, comments


# ------
# run
# ------

savedata = 1
filename = '../CORSICA_files/P_DTSS04N6R3_ITER_MR_01000.TXT'

data, comments = preader(filename)

# normalized poroidal flux
psi = data['Psi']
psin = (psi - psi[0]) / (psi[-1] - psi[0])

# ------
# save psin, Ti, Te, Ni, Ne, er, p0 to netCDF file
# which is required for map_pfile2grid.map_nc2grid()
# ------

if savedata:
    ofile = filename.replace('TXT', 'nc')
    print("save data to {}".format(ofile))
    odata = {}
    odata['psin'] = psin
    odata['er'] = np.zeros_like(psin)
    odata['ti'] = data['Ti']                # eV
    odata['te'] = data['Te']
    odata['ni'] = data['Ni'] / 1.e20        # 1e20 / m^3
    odata['ne'] = data['Ne'] / 1.e20
    #odata['p0'] = data['Pressure_Thermal']  # Pa
    odata['p0'] = data['Ni'] * (data['Ti'] + data['Te']) * 1.6021765e-19
    save2nc(ofile, 'w', **odata)

# -----
# save all data to netCDF
# WARNING: some vars' names contain special chars, failed to save
# -----
if savedata == 'all':
    ofile = filename.replace('.TXT', '_all.nc')
    print("save all data to {}".format(ofile))
    data['psin'] = psin
    save2nc(ofile, 'w', **data)


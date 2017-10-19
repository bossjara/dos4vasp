import numpy as np
import matplotlib.pyplot as plt


def main():
    tdos, pdoses, natom, title, efermi, npoints = read_doscar()
    tdos = np.array(tdos, dtype=np.float)
    brtdos = np.zeros([npoints, 3], dtype=np.float)
    tdos[:, 0] = tdos[:, 0] - efermi
    write_dos(tdos, dosfile='tdos.dat')
    doses = np.array(pdoses)
    Mn_s = np.zeros([npoints, 3], dtype=np.float)
    Mn_p = np.zeros([npoints, 3], dtype=np.float)
    Mn_d = np.zeros([npoints, 3], dtype=np.float)
    Ge_p = np.zeros([npoints, 3], dtype=np.float)
    Mn_s[:, 0] = Mn_p[:, 0] = Mn_p[:, 0] = Mn_d[:, 0] = tdos[:, 0]
    Ge_p[:, 0] = tdos[:, 0]
    for i in range(doses.shape[1]):
        Mn_s[i, 1] = doses[0, i, 1] + doses[1, i, 1]
        Mn_s[i, 2] = doses[0, i, 2] + doses[1, i, 2]
        Mn_p[i, 1] = doses[0, i, 3] + doses[0, i, 4] + doses[0, i, 5] +\
                     doses[1, i, 3] + doses[1, i, 4] + doses[1, i, 5]
        Mn_p[i, 2] = doses[0, i, 6] + doses[0, i, 7] + doses[0, i, 8] +\
                     doses[1, i, 6] + doses[1, i, 7] + doses[1, i, 8]

        Mn_d[i, 1] = doses[0, i, 9] + doses[0, i, 10] + doses[0, i, 11] +\
                     doses[0, i, 12] + doses[0, i, 13] + doses[1, i, 9] + \
                     doses[1, i, 10] + doses[1, i, 11] + doses[1, i, 12] + \
                     doses[1, i, 13]
        Mn_d[i, 2] = doses[0, i, 14] + doses[0, i, 15] + doses[0, i, 16] + \
                     doses[0, i, 17] + doses[0, i, 18] + doses[1, i, 14] + \
                     doses[1, i, 15] + doses[1, i, 16] + doses[1, i, 17] + \
                     doses[1, i, 18]
        Ge_p[i, 1] = doses[19, i, 6] + doses[19, i, 7] + doses[19, i, 8] +\
                     doses[20, i, 6] + doses[20, i, 7] + doses[20, i, 8] + \
                     doses[21, i, 6] + doses[21, i, 7] + doses[21, i, 8]
        Ge_p[i, 2] = doses[19, i, 9] + doses[19, i, 10] + doses[19, i, 11] +\
                     doses[20, i, 9] + doses[20, i, 10] + doses[20, i, 11] +\
                     doses[21, i, 9] + doses[21, i, 10] + doses[21, i, 11]
    brtdos[:, 0] = tdos[:, 0]
    brtdos[:, 1] = calc_spectra(tdos[:, [0, 1]])
    brtdos[:, 2] = calc_spectra(tdos[:, [0, 2]])
    write_dos(brtdos, dosfile='brtdos.dat')
    brtdos[:, 1] = calc_spectra(Mn_s[:, [0, 1]])
    brtdos[:, 2] = calc_spectra(Mn_s[:, [0, 2]])
    write_dos(brtdos, dosfile='Mn12-s.dat')
    brtdos[:, 1] = calc_spectra(Mn_p[:, [0, 1]])
    brtdos[:, 2] = calc_spectra(Mn_p[:, [0, 2]])
    write_dos(brtdos, dosfile='Mn12-p.dat')
    brtdos[:, 1] = calc_spectra(Mn_d[:, [0, 1]])
    brtdos[:, 2] = calc_spectra(Mn_d[:, [0, 2]])
    write_dos(brtdos, dosfile='Mn12-d.dat')
    brtdos[:, 1] = calc_spectra(Ge_p[:, [0, 1]])
    brtdos[:, 2] = calc_spectra(Ge_p[:, [0, 2]])
    write_dos(brtdos, dosfile='Ge123-p.dat')


def read_doscar(doscar='DOSCAR'):
    with open(doscar, 'r') as infile:
        natom = int(infile.readline().split()[0])
        print(natom)
        infile.readline()
        infile.readline()
        infile.readline()
        title = infile.readline().strip('\n')
        tmp = infile.readline().split()
        npoints = int(tmp[2])
        efermi = float(tmp[3])
        print(title)
        print(npoints)
        tdos = []
        for _ in range(npoints):
            tdos.append([float(qq) for qq in infile.readline().split()])
        pdoses = []
        for _ in range(natom):
            npoints = int(infile.readline().split()[2])
            tmp = []
            for _ in range(npoints):
                tmp.append([float(qq) for qq in infile.readline().split()])
            pdoses.append(tmp)
    return tdos, pdoses, natom, title, efermi, npoints


def write_dos(dos, dosfile='dos.dat'):
    with open(dosfile, 'w') as f:
        if dos.shape[1] == 2:
            for i in range(dos.shape[0]):
                f.write(f'{dos[i, 0]: .4f} {dos[i, 1]: .4f}\n')
        else:
            for i in range(dos.shape[0]):
                f.write(f'{dos[i, 0]: .4f} {dos[i, 1]: .4f} {-dos[i, 2]:.4f}\n')


def calc_spectra(dos):
    for i in range(dos.shape[0]):
        if dos[i, 0] >= 0:
            dos[i, 1] = np.float(0)
    return convolute(dos[:, 0], dos[:, 1], lorentzian)


def lorentzian(e, gamma):
    return gamma/(e * e + gamma * gamma)


def convolute(x, f, fbr):
    fsmear = np.zeros(len(x), dtype=np.float)
    for i in range(len(x)):
        for j in range(len(x)):
            fsmear[i] += f[j] * fbr(x[j]-x[i], 0.1)
    return fsmear*np.pi/len(x)


class Atom:

    def __init__(self):
        pass


class EfAtom(Atom):

    def __init__(self, *args):
        pass


if __name__ == "__main__":
    main()

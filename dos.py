
def main():
    with open('DOSCAR', 'r') as infile:
        natom = int(infile.readline().split()[0])
        print(natom)
        infile.readline()
        infile.readline()
        infile.readline()
        title = infile.readline().strip('\n')
        npoints = int(infile.readline().split()[2])
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
        print(len(pdoses), len(pdoses[0]))


class Atom:

    def __init__(self):
        pass


class EfAtom(Atom):

    def __init__(self, *args):
        pass


if __name__ == "__main__":
    main()

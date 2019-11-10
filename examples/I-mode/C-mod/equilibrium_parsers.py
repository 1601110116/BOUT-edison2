import re
from typing import List


class PfileParser(object):
    """This class parses a p-file, returning the profiles in Lists. Originally designed
    using Alcator C-mod p1120907032.01010. Other pfiles might require modifications of
    the parser. The pfile is supposed to be in the current directory.

    Parameters
    ----------
    filename : The file name of the pfile in the current directory

    """

    def __init__(self, filename):
        self.filename = filename

    def get_psinorm(self):
        """
        Returns
        -------
        psinorm : normalized psi of all profiles in the pfile.

        """
        with open(self.filename, 'r') as f:
            f.readline()  # skip the first line
            psinorm: List[float] = []
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse pfile: %s' % self.filename)
                if re.match(r'.*psinorm.*', line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                psinorm.append(float(values[0]))
        return psinorm

    def get_ne(self):
        """
        Returns
        -------
        ne : List of the electron density profile
        unit : unit of the values in the ne List
        """
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse pfile: %s' % self.filename)
                groups = re.match(r'^(.*ne\(10\^)(\d+)(/)(.*)(\).*)$', line)
                if groups:
                    order = float(groups.group(2))
                    unit = groups.group(4)
                    break
            ne: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                ne.append(10 ** order * float(values[1]))
        return ne, unit

    def get_ni(self):
        """
        Returns:
        --------
        ni : List of the ion density profile
        unit : unit of the values in the ni List
        """
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse pfile: %s' % self.filename)
                groups = re.match(r'^(.*ni\(10\^)(\d+)(/)(.*)(\).*)$', line)
                if groups:
                    order = float(groups.group(2))
                    unit = groups.group(4)
                    break
            ni: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                ni.append(10 ** order * float(values[1]))
        return ni, unit

    def get_te(self):
        """
        Returns
        -------
        te : List of the electron temperature profile
        unit : unit of the values in the te List
        """
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse pfile: %s' % self.filename)
                groups = re.match(r'^(.*te\()(\w+)(\).*)$', line)
                if groups:
                    unit = groups.group(2)
                    break
            te: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                te.append(float(values[1]))
        return te, unit

    def get_ti(self):
        """
        Returns
        -------
        ti : List of the ion temperature profile
        unit :  unit of the values in the ti List
        """
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse pfile: %s' % self.filename)
                groups = re.match(r'^(.*ti\()(\w+)(\).*)$', line)
                if groups:
                    unit = groups.group(2)
                    break
            ti: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                ti.append(float(values[1]))
        return ti, unit


class CXRSParser(object):
    """This class parses a CXRS file, returning the profiles in Lists. Originally designed
    using Alcator C-mod cxrs1120907032.01010_v20140623.txt. Other CXRS files might require
    modifications of the parser. CXRS file is supposed to be in the current directory.

    Parameters
    ----------
    filename : The file name of the CXRS file in the current directory

    """

    def __init__(self, filename):
        self.filename = filename

    def get_tB(self):
        """
        Returns
        -------
        psinorm : normalized psi of the profile
        tB : List of the Boron temperature profile
        tB_err : error of the Boron temperature profile
        unit : unit of tB and tB_err
        """
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse CXRS file: %s' % self.filename)
                groups = re.match(r'^(.*T_.*\()(.+)(\).*)$', line)
                if groups:
                    unit = groups.group(2)
                    break
            psinorm: List[float] = []
            tB: List[float] = []
            tB_err: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                if len(values) == 3:
                    psinorm.append(float(values[0]))
                    tB.append(float(values[1]))
                    tB_err.append(float(values[2]))
        return psinorm, tB, tB_err, unit

    def get_nB(self):
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse CXRS file: %s' % self.filename)
                groups = re.match(r'^(.*n_.*\(10\^)(\d+)(\s+)(.+)(\).*)$', line)
                if groups:
                    order = float(groups.group(2))
                    unit = groups.group(4)
                    break
            psinorm: List[float] = []
            nB: List[float] = []
            nB_err: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                if len(values) == 3:
                    psinorm.append(float(values[0]))
                    nB.append(10 ** order * float(values[1]))
                    nB_err.append(10 ** order * float(values[2]))
        return psinorm, nB, nB_err, unit

    def get_vpolB(self):
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse CXRS file: %s' % self.filename)
                groups = re.match(r'^(.*vpol_.*\()(.+)(\).*)$', line)
                if groups:
                    unit = groups.group(2)
                    break
            psinorm: List[float] = []
            vpolB: List[float] = []
            vpolB_err: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                if len(values) == 3:
                    psinorm.append(float(values[0]))
                    vpolB.append(float(values[1]))
                    vpolB_err.append(float(values[2]))
        return psinorm, vpolB, vpolB_err, unit

    def get_vtorB(self):
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse CXRS file: %s' % self.filename)
                groups = re.match(r'^(.*vtor_.*\()(.+)(\).*)$', line)
                if groups:
                    unit = groups.group(2)
                    break
            psinorm: List[float] = []
            vtorB: List[float] = []
            vtorB_err: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                if len(values) == 3:
                    psinorm.append(float(values[0]))
                    vtorB.append(float(values[1]))
                    vtorB_err.append(float(values[2]))
        return psinorm, vtorB, vtorB_err, unit

    def get_Er(self):
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    raise EOFError('Failed to parse CXRS file: %s' % self.filename)
                groups = re.match(r'^(.*Er.*\()(.+)(\).*)$', line)
                if groups:
                    unit = groups.group(2)
                    break
            psinorm: List[float] = []
            Er: List[float] = []
            Er_err: List[float] = []
            while True:
                line = f.readline()
                if re.match(r'.*psinorm.*', line) or (not line):
                    break
                line = line.strip()
                values = re.split(r'\s+', line)
                if len(values) == 3:
                    psinorm.append(float(values[0]))
                    Er.append(float(values[1]))
                    Er_err.append(float(values[2]))
        return psinorm, Er, Er_err, unit

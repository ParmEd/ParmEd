import numpy as np

class Symmetry(object):
    def __init__(self, lines=None):
        if lines is not None:
            self.data = self.parse(lines)
        else:
            self.data = []

    def parse(self, lines):
        data = []
        for line in lines:
            if line.startswith('REMARK 290   SMTRY'):
                data.append(line.split()[4:])

        return np.asarray(data, dtype='f8')

    def write_string(self):
        # rename to write?
        fmt = '%d%4d%10.6f%10.6f%10.6f%15.5f'

        lines = []

        for index, arr in enumerate(self.data):
            arr_list = [1 + index % 3, 1 + int(index/3)] + arr.tolist()
            line = "REMARK 290   SMTRY" + fmt % tuple(arr_list)
            lines.append(line)

        output = '\n'.join(lines)
        return output

if __name__ == '__main__':
    lines = """
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000       36.30027
REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000       59.50256
REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       46.45545
REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       59.50256
REMARK 290   SMTRY1   4  1.000000  0.000000  0.000000       36.30027
REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       46.45545
REMARK 290   SMTRY3   4  0.000000  0.000000 -1.000000        0.00000
""".split('\n')

    symm = Symmetry(lines)
    data = symm.data

    symm_string = symm.write_string()
    print(symm_string)

    for line in lines:
        assert line.strip() in symm_string

    print(symm.data)

import re
import numpy as np
import h5py


def extract(inp, key, shape, output):
    f = open(inp, 'r')

    pattern = re.compile(r'(\d+)\s+(\d+)\s+(\d+)\s+([\d.Ee+\-]+)\s+([\d.Ee+\-]+)')

    count = shape[0] * shape[1] * shape[2]

    data = np.array([]).reshape((0, shape[0], shape[1], shape[2], 2))
    step = 1
    c = 0
    while True:
        line = f.readline()
        c += 1
        if not line:
            break
        if key in line:
            f.readline()
            c += 1
            values = np.zeros(shape=(1, shape[0], shape[1], shape[2], 2), dtype=float)
            for i in range(count):
                line = f.readline()
                c += 1
                m = pattern.search(line)
                if m:
                    values[0, int(m.group(1)) - 1, int(m.group(2)) - 1, int(m.group(3)) - 1, 0] = float(m.group(4))
                    values[0, int(m.group(1)) - 1, int(m.group(2)) - 1, int(m.group(3)) - 1, 1] = float(m.group(5))
                else:
                    print("not parsed: {}".format(line.strip()))
            data = np.concatenate([data, values], axis=0)
            print("Step {} finished, line count: {}".format(step, c))
            step += 1

    f.close()

    with h5py.File(output, 'w') as f:
        f.create_dataset(name='tally', data=data)


if __name__ == "__main__":
    extract("problem8_inp.RMC.Tally", "ID = 5,  Type = power, Number of mesh grids  = 3186225",
            (255, 255, 49), "tally.h5")


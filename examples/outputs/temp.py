import pickle
from pathlib import Path


if __name__ == "__main__":
    dir = Path.cwd()
    path = dir / "msa/MSA/POE-0.1-in-50-years/Record1nodeDisp.pickle"
    path1 = dir / "msa/MSA_old/POE-0.1-in-50-years/Record1nodeDisp.pickle"

    with open(path, "rb") as f:
        data = pickle.load(f)

    with open(path1, "rb") as f:
        data1 = pickle.load(f)

    for i in range(len((data[0][0][0]))):
        print(data[0][0][1][i], data1[0][0][1][i])

from __future__ import print_function
import os
import uproot as ur

def main():
    for file in os.scandir("/sphenix/tg/tg01/hf/zhji/data/sphenix/histos"):
        if (file.name.startswith("training-") and
            file.name.endswith(".root") and
            file.is_file()):
            try:
                ur.dask(file.path + ":T")
            except ur.exceptions.KeyInFileError:
                print("bad " + file.name)
            else:
                os.system("cp " + file.path + " /phenix/spin/phnxsp01/zji/data/sphenix/histos")

if __name__ == '__main__':
    main()

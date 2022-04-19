from __future__ import print_function
import os
import uproot as ur

def main():
    for file in os.scandir("/phenix/spin/phnxsp01/zji/data/sphenix/histos"):
        if (file.name.startswith("training-") and
            file.name.endswith(".root") and
            file.is_file()):
            try:
                ur.lazy(file.path + ":T")
            except ur.exceptions.KeyInFileError:
                print("rm " + file.name)
                os.system("rm " + file.path)

if __name__ == '__main__':
    main()

import os
import glob

npp = glob.glob("data/S*")



indexFile = open("sampleIndex.csv", "w")

for pp in npp:
    a1 = os.path.basename(pp)
    a2 = os.path.relpath(pp)
    indexFile.write(a1)
    indexFile.write(",")
    indexFile.write(a2)
    indexFile.write("\n")

indexFile.close()

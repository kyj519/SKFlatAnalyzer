#!/usr/bin/env python3

fin = open("reshaping_deepJet_v3.csv")
fout = open("reshaping_deepJet_v3_simplified.csv", "w")

while True:
    line = fin.readline()
    if not line: break
        
    if "_jes" not in line:
        fout.write(line)
    else:
        if "_jes," in line:
            fout.write(line)
fin.close()
fout.close()

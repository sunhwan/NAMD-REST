import os
for i in range(32):
    for line in open('output/%d/dist.job0.%d.history' % (i,i)):
        entries = line.strip().split()
        if entries[0] == '154000':
            rep = int(entries[1])
            os.system('~/catdcd -stype psf -s input/separation.psf -otype pdb -o input/input.%d.pdb -first 77 -last 77 -dcd output/%d/dist.job0.%d.dcd' % (rep, i, i))

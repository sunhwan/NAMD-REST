import os

nrep = 32
job = 1
images = [None for i in range(nrep)]
for i in range(nrep):
    history = 'output/%d/dist.job%d.%d.history' % (i, job, i)
    time, rep = open(history).readlines()[-1].strip().split()[:2]
    images[int(rep)] = i

for i in range(nrep):
    rep = images[i]
    print """
mol new input/solvated.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile output/%(rep)d/dist.job%(job)d.%(rep)d.coor type namdbin first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top [format "Image %(i)d"]
mol delrep 0 top""" % locals()

    os.system('~/catdcd -s input/solvated.psf -stype psf -otype pdb -o pdbs/image.%(i)02d.pdb -namdbin output/%(rep)d/dist.job%(job)d.%(rep)d.coor' % locals())


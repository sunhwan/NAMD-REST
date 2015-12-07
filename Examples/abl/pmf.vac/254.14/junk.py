import os
import glob
import commands

nx = 12
ny = 4
nset = 10
for j in range(ny):
    #os.system('mkdir -p data/%d' % j)
    #os.system('rm -rf data/%d/*' % j)

    for i in range(nx):
        ix = j * nx + i
        for traj in glob.glob('output/%d/dist.job*.sort.colvars.traj' % ix):
            print traj, 'data/%d/%d' % (j, i)
            #os.system("""awk '{if ($1!="#") print $1, $4}' %s >> data/%d/%d.traj""" % (traj, j, i))
            #os.system("cat %s >> data/%d/%d.sort.colvars.traj""" % (traj, j, i))

            hist = traj[:-12] + 'history'
            #os.system("cat %s >> data/%d/%d.sort.history""" % (hist, j, i))

        nl = int(commands.getoutput('wc -l data/%d/%d.traj' % (j, i)).split()[0])
        #os.system("""awk 'FNR==NR {t[NR]=$4; off=NR; next} {i=NR-off; if (i%%2==0) print $1,$2,t[i/2]}' data/%d/%d.sort.history data/%d/%d.traj >> data/%d/%d.temp""" % (j, i, j, i, j, i))

        for k in range(1, nset):
            #os.system('mkdir -p data/%d/%d' % (j, k))
            #open('data/%d/%d/%d.traj' % (j, k, i), 'w').write("".join(open('data/%d/%d.traj' % (j, i)).readlines()[int(float(k)/nset*nl):]))
            #open('data/%d/%d/%d.sort.history' % (j, k, i), 'w').write("".join(open('data/%d/%d.sort.history' % (j, i)).readlines()[int(float(k)/nset*nl):]))
            open('data/%d/%d/%d.temp' % (j, k, i), 'w').write("".join(open('data/%d/%d.temp' % (j, i)).readlines()[int(float(k)/nset*nl/2):]))

        #nl = int(commands.getoutput('wc -l data/%d/%d.traj' % (j, i)).split()[0])
        #os.system('mkdir -p data/%d/last/' % (j))
        #open('data/%d/last/%d.traj' % (j, i), 'w').write("".join(open('data/%d/%d.traj' % (j, i)).readlines()[int(nl*0.5):]))

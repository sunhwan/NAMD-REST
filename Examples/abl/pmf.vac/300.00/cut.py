import os
for i in range(15):
    print i
    os.system('cut -d " " -f 2 output/%d/dist.job*.%d.history > history/%d.history' % (i,i,i))

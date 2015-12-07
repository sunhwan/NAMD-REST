import glob
total = 0
count = 0
runs = 0
nx = 12
ny = 4
replica_x = [0 for i in range(nx*ny)]
replica_y = [0 for i in range(nx*ny)]

for i in range(nx*ny):
    prev = None
    _runs = 0
    for history in sorted(glob.glob('output/*/*.%d.history'%i)):
        print history
        for line in open(history):
            cnt = int(line.split()[1])

            if prev is None:
                prev = cnt
                continue
            if prev != cnt:
                if abs(cnt - prev) == nx: count = replica_y
                else: count = replica_x
                count[prev] += 1
                count[cnt] += 1
                prev = cnt
            _runs += 1
if runs == 0: runs = _runs

for i in range(nx*ny):
    print "%4d %8.2f %8.2f" % (i, replica_x[i] / float(runs) / 2 * 100, replica_y[i] /float(runs) / 2 * 100)

import sys

prefix = sys.argv[1]
step = int(sys.argv[2])
num_replica = int(sys.argv[3])
final_step = -1
if len(sys.argv) > 4:
    final_step = int(sys.argv[4])

history_fps = [open('output/%d/%s.job%d.%d.history' % (i, prefix, step, i)) for i in range(num_replica)]
colvars_fps = [open('output/%d/%s.job%d.%d.colvars.traj' % (i, prefix, step, i)) for i in range(num_replica)]

sorted_history_fps = [open('output/%d/%s.job%d.%d.sort.history' % (i, prefix, step, i), 'w') for i in range(num_replica)]
sorted_colvars_fps = [open('output/%d/%s.job%d.%d.sort.colvars.traj' % (i, prefix, step, i), 'w') for i in range(num_replica)]

cnt = 1
while 1:
    for i in range(num_replica):
        try:
            line = history_fps[i].readline()
            time, rep = map(int, line.split()[:2])
        except:
            final_step = time
            break
        sorted_history_fps[rep].write(line)

        while 1:
            line = colvars_fps[i].readline()
            if line.strip().split()[0] != str(time): continue
            sorted_colvars_fps[rep].write(line)
            break

    if final_step > 0 and time >= final_step: break
    print time


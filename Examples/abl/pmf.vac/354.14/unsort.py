
fps = []
for i in range(48):
    fps.append(open('../%d/dist.job0.%d.sort.history' % (i, i)))

f = open('1.history', 'w')
count = 0
cnt = 1
reps = [0 for i in range(48)]

while 1:
    lines = []
    for i in range(48):
        line = fps[i].readline()
        if not line: sys.exit()
        lines.append(line)

    rep = int(lines[cnt].strip().split()[1])
    if rep != cnt:
        cnt = reps[rep]

    entries = lines[cnt].strip().split()
    entries[1] = str(cnt)
    print " ".join(entries)
    f.write(" ".join(entries) + "\n")
    count += 1

    for i in range(48):
        rep = int(lines[i].strip().split()[1])
        reps[rep] = i


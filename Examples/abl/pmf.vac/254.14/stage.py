import commands
import sys
icnt = int(sys.argv[1])
states = commands.getoutput("grep -c 'stage' output/*/dist.job%d.*colvars.state | grep -v :2" % icnt)
for f in states.splitlines():
    f = f.split(':')[0]
    state = []
    flag = False
    for line in open(f).readlines():
        if line.startswith('harmonic'):
            flag = True
        if flag and line.strip().startswith('stage'):
            flag = False
        if flag and line.strip().startswith('}'):
            state.append('stage 0')
            flag = False
        
        state.append(line.rstrip())

    open(f, 'w').write("\n".join(state))


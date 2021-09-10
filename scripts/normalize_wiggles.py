import sys


read_total = float(sys.argv[1])/1000000

for line in sys.stdin:
    line = line.replace('\n','').split(' ')
    line[3] = str(int(line[3])/read_total)
    print('\t'.join(line))



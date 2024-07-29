s = []
x = [[] for i in range(4)]

def readrec(line):
    #if line.startswith('a'):
    # the header line could be anything, assume we're at the beginning when called
    s.append(line)
    line=file.readline()
    while not(line.startswith('b')):
        x[irec].append(int(line))
        line=file.readline()

irec=0
with open('dummydata.txt','r') as file:
    line=file.readline()
    while not(line.startswith('c')):
        readrec(line) # repeatedly reads records
        irec+=1
        line=file.readline()

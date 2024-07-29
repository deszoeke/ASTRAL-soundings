#import string as str

lines = [' ggg', 'hhh ', ' ', '  ', '\n', '', 'iii', 'last']

i = 0
line = str.strip( lines[i] )
while not(line.startswith('last')):
    if (line != ''):
        print(line)
    i = i + 1
    line = str.strip( lines[i] )

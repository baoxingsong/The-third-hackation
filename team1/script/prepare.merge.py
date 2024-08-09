import pandas as pd

block_id = []
with open('synteny.txt', mode='r') as f:
    for line in f:
        line = line.strip()
        if len(line.split(' ')) < 6:
            continue
        else:
            first = line.split(' ')[0]
            block_id.append(first.split(':')[0])

with open('sequenceColor.txt', mode='r') as f1,open('mgr_macro.txt', mode='w') as f2:
    for line in f1:
        line = line.strip()
        n = 0
        first_id = -1
        for i in line.split(' '):
            if n % 2 == 0:
                n += 1
                continue
            elif str(i) not in block_id:
                n += 1
                continue
            elif i == first_id:
                n += 1
                continue
            else:
                n += 1
                first_id = i
                f2.write(f'{i} ')
        f2.write('\n')


def readTextFile(filePath):
    with open(filePath, 'r') as f:
        return "".join([l.strip() for l in f.readlines()])

def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq +'\n')

def readFASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAfile = [l.strip() for l in f.readlines()]

    FASTAdict = {}
    FASTAlabel = ""

    for line in FASTAfile:
        if '>' in line:
            FASTAlabel = line
            FASTAdict[FASTAlabel] = ""
        else:
            FASTAdict[FASTAlabel] += line

    return FASTAdict
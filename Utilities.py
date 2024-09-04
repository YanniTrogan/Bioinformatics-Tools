#Document for python utilities
def colored(seq):
    """Colorizes DNA nucleotide bases"""
    bcolors = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset': '\033[0;0m'
    }

    tmpStr = ""

    for nuc in seq: 
        if nuc in bcolors:
            tmpStr += bcolors[nuc] + nuc
        else:
            tmpStr += bcolors['reset'] + nuc
    
    return tmpStr + '\033[0;0m'

def Fibonacci_Loop(months, offspring):
    parent = 1
    child = 1
    for int in range(months - 1):
        tmpVal = child
        child = parent
        parent = parent + tmpVal
    return child

def Fibonacci_Loop_Python(months, offspring):
    parent, child = 1, 1
    for int in range(months - 1):
        child, parent = parent, parent + (child * offspring)
    return child

def readTextFile(filepath):
    """Reading a file and returning a list of lines"""
    with open(filepath, 'r') as f:
        return "".join([l.strip() for l in f.readlines()])

def writeTextFile(filePath, seq, mode ='w'):
    with open(filePath, mode) as f:
        f.write(seq + '\n')

def read_Fasta(filePath):
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]
    FASTADict = {}
    FASTALabel = ""
    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line
    return FASTADict
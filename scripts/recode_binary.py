import sys, os

def recode_binary(fl):
    newfl = open("BIN_"+fl,"w")
    aln = open(fl,"r")
    h = aln.readline()
    newfl.write(h)
    for i in aln:
        if len(i.strip()) == 0:
            continue
        spls = i.strip().split()
        nm = spls[0]
        seq = spls[1]
        newseq = ""
        for site in seq:
            if site == "G" or site == "A":
                newseq += "0"
            elif site == "C" or site == "T":
                newseq += "1"
        newfl.write(nm+"\t"+newseq+"\n")

def recode_4_state(fl):
    newfl = open("recode_"+fl,"w")
    aln = open(fl,"r")
    h = aln.readline()
    newfl.write(h)
    for i in aln:
        if len(i.strip()) == 0:
            continue
        spls = i.strip().split()
        nm = spls[0]
        seq = spls[1]
        newseq = ""
        for site in seq:
            if site == "G":
                newseq += "0"
            elif site == "C":
                newseq += "1"
            elif site == "A":
                newseq+= "2"
            elif site == "T":
                newseq += "3"
            elif site == "-" or site == "?":
                newseq += site
        newfl.write(nm+"\t"+newseq+"\n")

def single_fl(fl,states):
    if states == "2":
        recode_binary(fl)
    elif states == "4":
        recode_4_state(fl)

def multi_fl(directory,states):
    if states == "2":
        for i in os.listdir(directory):
            recode_binary(directory+i)
    elif states == "4":
        for i in os.listdir(directory):
            recode_4_state(directory+i)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "usage: " + sys.argv[0]+ " <alignment(s)> <nstates to recode to>"
        sys.exit(0)

    single_fl(sys.argv[1],str(sys.argv[2]))
    #multi_fl(sys.argv[1],str(sys.argv[2]))

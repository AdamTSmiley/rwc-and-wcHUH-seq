TotalDataExport = open("rwcHUHSeqData.txt", "w")
Mn = open("rwcHUHSeqMnCl2.txt", "w")
Mg = open("rwcHUHSeqMgCl2.txt", "w")

from Bio import SeqIO

def ForwardReads(file):
    Reads = []
    open(file)
    for record in SeqIO.parse(file, "fastq"):
        temp = str(record.seq)
        Reads.append(temp)
    
    TrimmedData = []
    
    for i in range(len(Reads)):
        if len(Reads[i]) == 85:
            for z in range(len(Reads[i])):
                if Reads[i][z:z+3] == "TAT":
                    if Reads[i][z+4:z+7] == "ATT":
                        if Reads[i][z+7] != "N":
                            if Reads[i][z+8:z+10] == "CC":
                                TrimmedData.append(Reads[i][z:z+10])
                                
    Seqs = {}

    for i in range(len(TrimmedData)):
        if TrimmedData[i] in Seqs:
            Seqs[TrimmedData[i]] += 1
        else:
            Seqs[TrimmedData[i]] = 1
    return(Seqs)

def DataExport(D1,D2,D3,D4,D5,D6):
    print("Sequence" + "\t" + "MgCl2_1" + "\t" + "MgCl2_2" + "\t" + "MgCl2_3" + "\t" + "MnCl2_1" + "\t" + "MnCl2_2" + "\t" + "MnCl2_3", file = TotalDataExport)
    for i in sorted(D1):
        print(i + "\t" + str(D1[i]) + "\t" + str(D2[i]) + "\t" + str(D3[i]) + "\t" + str(D4[i]) + "\t" + str(D5[i]) + "\t" + str(D6[i]), file = TotalDataExport)

def Prop(D1,D2,D3):
    counts = {}
    for i in sorted(D1):
        counts[i] = D1[i] + D2[i] + D3[i]
    
    total = 0
    
    for i in counts:
        total += counts[i]
    
    props = {}
    
    for i in counts:
        props[i] = ((counts[i]/total)*100)
        
    return props

###############################################################################

Mg1 = ForwardReads("MgCl21_R1_001.fastq")
Mg2 = ForwardReads("MgCl22_R1_001.fastq")
Mg3 = ForwardReads("MgCl23_R1_001.fastq")

Mn1 = ForwardReads("MnCl21_R1_001.fastq")
Mn2 = ForwardReads("MnCl22_R1_001.fastq")
Mn3 = ForwardReads("MnCl23_R1_001.fastq")

###############################################################################

MnProp = Prop(Mn1, Mn2, Mn3)
MgProp = Prop(Mg1, Mg2, Mg3)

out = DataExport(Mg1, Mg2, Mg3, Mn1, Mn2, Mn3)

###############################################################################

print("Sequence" + "\t" + "Percent Enrichment", file = Mg)
for i in MgProp:
    print(i + "\t" + str(MgProp[i]), file = Mg)
    
print("Sequence" + "\t" + "Percent Enrichment", file = Mn)
for i in MnProp:
    print(i + "\t" + str(MnProp[i]), file = Mn)

###############################################################################

TotalDataExport.close()
Mn.close()
Mg.close()
TotalDataExport = open("wcHUHSeqCleavageData.txt", "w")
Mn = open("wcHUHSeqMnCl2PercentReductions.txt", "w")
Mg = open("wcHUHSeqMgCl2PercentReductions.txt", "w")

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

def DataExport(D1,D2,D3,D4,D5,D6,D7,D8,D9):
    print("Sequence" + "\t" + "Ref1" + "\t" + "Ref2" + "\t" + "Ref3" + "\t" + "MgCl2_1" + "\t" + "MgCl2_2" + "\t" + "MgCl2_3" + "\t" + "MnCl2_1" + "\t" + "MnCl2_2" + "\t" + "MnCl2_3", file = TotalDataExport)
    for i in sorted(D1):
        print(i + "\t" + str(D1[i]) + "\t" + str(D2[i]) + "\t" + str(D3[i]) + "\t" + str(D4[i]) + "\t" + str(D5[i]) + "\t" + str(D6[i]) + "\t" + str(D7[i]) + "\t" + str(D8[i]) + "\t" + str(D9[i]), file = TotalDataExport)

def PercentReductions(R1,R2,R3,D1,D2,D3):
    RefCounts = {}
    for i in sorted(R1):
        temp = 0
        temp = R1[i]+R2[i]+R3[i]
        RefCounts[i] = temp
    
    ConditionCounts = {}
    for i in sorted(D1):
        temp = 0
        temp = D1[i]+D2[i]+D3[i]
        ConditionCounts[i] = temp
    
    PercentReductions = []
    for i in sorted(RefCounts):
        temp = []
        temp.append(i)
        if ((RefCounts[i] - ConditionCounts[i])/RefCounts[i]) > 0:
            temp.append(round((RefCounts[i] - ConditionCounts[i])/RefCounts[i],3)*100)
            PercentReductions.append(temp)
        else:
            temp.append(0)
            PercentReductions.append(temp)
    return(PercentReductions) 

###############################################################################

MgCl2_1 = ForwardReads("MgCl2-1_R1_001.fastq")
MgCl2_2 = ForwardReads("MgCl2-2_R1_001.fastq")
MgCl2_3 = ForwardReads("MgCl2-3_R1_001.fastq")
MnCl2_1 = ForwardReads("MnCl2-1_R1_001.fastq")
MnCl2_2 = ForwardReads("MnCl2-2_R1_001.fastq")
MnCl2_3 = ForwardReads("MnCl2-3_R1_001.fastq")
Ref_1 = ForwardReads("Ref-1_R1_001.fastq")
Ref_2 = ForwardReads("Ref-2_R1_001.fastq")
Ref_3 = ForwardReads("Ref-3_R1_001.fastq")

###############################################################################

DataExport(Ref_1,Ref_2,Ref_3,MgCl2_1,MgCl2_2,MgCl2_3,MnCl2_1,MnCl2_2,MnCl2_3)

###############################################################################

MgCl2PercentReductions = PercentReductions(Ref_1,Ref_2,Ref_3,MgCl2_1,MgCl2_2,MgCl2_3)
MnCl2PercentReductions = PercentReductions(Ref_1,Ref_2,Ref_3,MnCl2_1,MnCl2_2,MnCl2_3)

###############################################################################

print("Sequence" + "\t" + "Percent Reductions", file = Mg)
for i in MgCl2PercentReductions:
    print(i[0] + "\t" + str(i[1]), file = Mg)
    
print("Sequence" + "\t" + "Percent Reductions", file = Mn)
for i in MnCl2PercentReductions:
    print(i[0] + "\t" + str(i[1]), file = Mn)

###############################################################################    

TotalDataExport.close()
Mn.close()
Mg.close()
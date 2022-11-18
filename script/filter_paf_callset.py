import sys,os
dd = {}
fout = open(sys.argv[3],'w')
with open(sys.argv[1],'r') as fp: #het variants from Deepvariant
    for line in fp:
        if line.startswith('#'):continue
        else:
            line_temp = line.strip().split('\t')
            dd.setdefault(line_temp[0],{})
            dd[line_temp[0]][int(line_temp[1])] = None
with open(sys.argv[2],'r') as fp: #paf vcv
    for line in fp:
        if line.startswith('#'):
            fout.write(line)
        else:
            line_temp = line.strip().split('\t')
            if line_temp[0] not in dd:continue
            if int(line_temp[1]) not in dd[line_temp[0]]:continue
            fout.write(line)
fout.close()


'''
chr1_RagTag     124757557       .       G       T       63.8    PASS    .       GT:GQ:DP:AD:VAF:PL      0/1:63:40:20,20:0.5:63,0,69
chr1_RagTag     124759716       .       G       GTT     0       RefCall .       GT:GQ:DP:AD:VAF:PL      0/0:36:37:20,5:0.135135:0,35,55
chr1_RagTag     124760047       .       A       G       59.9    PASS    .       GT:GQ:DP:AD:VAF:PL      0/1:59:39:20,18:0.461538:59,0,65
'''

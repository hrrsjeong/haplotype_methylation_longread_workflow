import sys,os,gzip
fout = open(sys.argv[2],'w')
with gzip.open(sys.argv[1],'rt') as fp:
    for line in fp:
        if line.startswith('#'):
            fout.write(line)
        else:
            line_temp = line.strip().split('\t')
            if line_temp[6] != "PASS":continue
            infos = line_temp[-1].split(':')
            if infos[0] != "0/1":continue
            if float(infos[1]) < 30:continue
            if int(infos[2]) < 20:continue
            if ',' in line_temp[4]:continue
            ref,alt = [int(x) for x in infos[3].split(',')]
            if ref < 10 or alt < 10:continue
            fout.write(line)
fout.close()


'''
chr1_RagTag     124757557       .       G       T       63.8    PASS    .       GT:GQ:DP:AD:VAF:PL      0/1:63:40:20,20:0.5:63,0,69
chr1_RagTag     124759716       .       G       GTT     0       RefCall .       GT:GQ:DP:AD:VAF:PL      0/0:36:37:20,5:0.135135:0,35,55
chr1_RagTag     124760047       .       A       G       59.9    PASS    .       GT:GQ:DP:AD:VAF:PL      0/1:59:39:20,18:0.461538:59,0,65
'''

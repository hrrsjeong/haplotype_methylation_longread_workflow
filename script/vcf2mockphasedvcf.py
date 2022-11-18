import sys
fout = open(sys.argv[1].replace(".filtered.vcf",".phased.vcf"),'w')
with open(sys.argv[1],'r') as fp:
    for line in fp:
        if line.startswith("#"):
            fout.write(line)
        else:
            line_temp = line.strip().split('\t')
            if line_temp[-1] != "1/1":continue
            fout.write('\t'.join(line_temp[0:7])+'\t.\tGT\t0|1\n')
fout.close()


'''
h2tg000001l	1254592	.	A	C	60	.	QNAME=h1tg000067l;QSTART=857904;QSTRAND=+	GT	1/1
h2tg000001l	1255105	.	G	C	60	.	QNAME=h1tg000067l;QSTART=858417;QSTRAND=+	GT	1/1
h2tg000001l	1255413	.	C	CT	60	.	QNAME=h1tg000067l;QSTART=858726;QSTRAND=+	GT	1/1
'''
    


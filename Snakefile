import pandas as pd
import Bio.SeqIO
import Bio.Seq
import collections
import gzip
import numpy as np
import os
import pandas as pd
import pysam
import re
import shutil

SDIR=os.path.dirname(workflow.snakefile)
print (SDIR)
#sys.exit()
CWD=os.getcwd()
#shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")
shell.prefix(f"source ~/.bashrc")
configfile: "config.yaml"
#deepvar_image = '/net/eichler/vol26/7200/software/pipelines/deepvariant/images/deepvariant_1.4.0-gpu.sif'
deepvar_image = config["deepvariant"] #'/net/eichler/vol26/7200/software/pipelines/deepvariant/images/deepvariant_1.4.0.sif'
sample_df = pd.read_csv('manifest.tab', sep='\t', index_col=['SAMPLE'])

#it requires hifi bam with kinetic info
dict_movie_bam   = {}
dict_tandem_pbsv = {}

for sample_idx in sample_df.index:
	dict_movie_bam.setdefault(sample_idx,{})
	dict_tandem_pbsv.setdefault(sample_idx,{})
	with open(sample_df.at[sample_idx,"BAM_fofn"],'r') as fp:
		for movie_file in fp:
			dict_movie_bam[sample_idx][movie_file.strip().split('/')[-1].split('.')[0]] = movie_file.strip()
	dict_tandem_pbsv[sample_idx]["hap1"] = sample_df.at[sample_idx,"hap1_TandemRepeat"]
	dict_tandem_pbsv[sample_idx]["hap2"] = sample_df.at[sample_idx,"hap2_TandemRepeat"]
print(dict_movie_bam)
#sys.exit()

localrules: run_deepvariant
localrules: run_pbcpg
localrules: run_pbcpg_hap

def gpu_check(deepvar_image):
	if "gpu" in deepvar_image:
		return "--nv"
	else:
		return ""

def find_asm(wildcards):
	return sample_df.at[wildcards.sample, wildcards.hap]

def movie_bam(wildcards):
	return dict_movie_bam[wildcards.sample][wildcards.movie]
	
def alt_hap_reads(wildcards):
	if wildcards.hap == "hap1":
		return '{sample}/align/readid/{sample}.movie_{movie}.haplotype_hap2.ref_hap2.readid.txt',
	elif wildcards.hap == "hap2":
		return '{sample}/align/readid/{sample}.movie_{movie}.haplotype_hap1.ref_hap1.readid.txt',
	else:
		print ("haplotype index conflict")
		sys.exit()

def tandem_repeat(wildcards):
	return dict_tandem_pbsv[wildcards.sample][wildcards.hap]

def whatshap_vcf(wildcards):
	if wildcards.hap == "hap1":
		return '{sample}/paf/{sample}.wga.hap2_to_hap1.phased.vcf.gz'
	elif wildcards.hap == "hap2":
		return '{sample}/paf/{sample}.wga.hap1_to_hap2.phased.vcf.gz'
	else:
		print ("haplotype index conflict")
		sys.exit()

def pbcpg_true_hap1_check(wildcards):
	if wildcards.hap == "hap1":
		return '{sample}/pb_cpg/{sample}.haplotagged.ref_hap1_pbcpg.hap1.reference.bed'
	elif wildcards.hap == "hap2":
		return '{sample}/pb_cpg/{sample}.haplotagged.ref_hap2_pbcpg.hap2.reference.bed'
	else:
		print ("haplotype pbcpg error")

def pbcpg_true_hap2_check(wildcards):
	if wildcards.hap == "hap1":
		return '{sample}/pb_cpg/{sample}.haplotagged.ref_hap1_pbcpg.hap2.reference.bed'
	elif wildcards.hap == "hap2":
		return '{sample}/pb_cpg/{sample}.haplotagged.ref_hap2_pbcpg.hap1.reference.bed'
	else:
		print ("haplotype pbcpg error")

def pbcpg_denovo_true_hap1_check(wildcards):
	if wildcards.hap == "hap1":
		return '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_hap1_pbcpg.hap1.denovo.bed'
	elif wildcards.hap == "hap2":
		return '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_hap2_pbcpg.hap2.denovo.bed'
	else:
		print ("haplotype pbcpg error")

def pbcpg_denovo_true_hap2_check(wildcards):
	if wildcards.hap == "hap1":
		return '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_hap1_pbcpg.hap2.denovo.bed'
	elif wildcards.hap == "hap2":
		return '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_hap2_pbcpg.hap1.denovo.bed'
	else:
		print ("haplotype pbcpg error")

alt_hap = {"hap1":"hap2","hap2":"hap1"}


rule all:
	input:
		#expand('{sample}/paf/{sample}.wga.hap2_to_hap1.vcf',sample=sample_df.index.get_level_values('SAMPLE')),
		#expand('{sample}/paf/{sample}.wga.hap1_to_hap2.vcf',sample=sample_df.index.get_level_values('SAMPLE')),
		#expand('{sample}/align/{sample}.mCG.ref_{hap}.bam.bai', sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/align/{sample}.haplotagged.ref_{hap}.bam.bai', sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/align/{sample}.readid_{hap}.haplotagged.ref_{hap}.bam.bai', sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/pbsv/{sample}.haplotagged.ref_{hap}.vcf',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/pbsv/{sample}.readid_{hap}.haplotagged.ref_{hap}.vcf',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/pb_cpg/{sample}.haplotagged.ref_{hap}_pbcpg-aligned_bam_to_cpg_scores.log',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/pb_cpg/{sample}.readid_{hap}.haplotagged.ref_{hap}_pbcpg-aligned_bam_to_cpg_scores.log',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		expand('{sample}/pb_cpg/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		expand('{sample}/pb_cpg/renamed/{sample}.readid_{hap}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		expand('{sample}/pb_cpg/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.{haplo}.reference.bed',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2'],haplo=['hap1','hap2']),
		#expand('{sample}/pb_cpg_denovo/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.combined.denovo.bed',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2']),
		#expand('{sample}/pb_cpg_denovo/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.{haplo}.denovo.bed',sample=sample_df.index.get_level_values('SAMPLE'), hap=['hap1','hap2'],haplo=['hap1','hap2']),

rule remove_short_reads:
	input:
		bam = movie_bam,
	output:
		outbam = temp('{sample}/bam/{movie}.mCG.bam')
	resources:
		mem = 12,
		hrs = 70,
	params:
		min_len = r'length(seq)>10000'
	threads: 4
	shell:
		'''
		source activate methylation_hifi
		samtools view -@{threads} -e '{params.min_len}' -o {output.outbam} {input.bam}
		'''

rule pbmm2_align:
	input:
		fa = find_asm,
		bam = '{sample}/bam/{movie}.mCG.bam'
	output:
		outbam = '{sample}/align/{sample}.movie_{movie}.mCG.ref_{hap}.bam'
	resources:
		mem = 12,
		hrs = 70,
	params:
		align_t = 4,
		sort_t = 2,
	threads: 6
	shell:
		'''
		source activate methylation_hifi
		pbmm2 align --preset HiFi --sort -j {params.align_t} -J {params.sort_t} {input.fa} {input.bam} {output.outbam}
		'''

rule minimap_asm:
	input:
		ref_fa =  lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		query_fa = lambda wildcards: sample_df.at[wildcards.sample, wildcards.alt]
	output:
		paf = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.paf'
	resources:
		mem = 12,
		hrs = 90,
	threads: 8
	shell:
		'''
		source activate methylation_hifi
		minimap2 -cx asm20 --cs --eqx --secondary=no -s 10000 -t {threads} {input.ref_fa} {input.query_fa} > {output.paf}
		'''

rule paf_call:
	input:
		paf = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.paf',
		ref = lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
	output:
		vcf = rules.minimap_asm.output.paf.replace(".paf",".vcf")
	resources:
		mem = 8,
		hrs = 12,
	params:
		min_aligned_len_cov = 50000,
		min_aligned_len_call = 100000,
	threads: 2
	shell:
		'''
		source activate methylation_hifi
		sort -k6,6 -k8,8n {input.paf} | paftools.js call -l {params.min_aligned_len_cov} -L {params.min_aligned_len_call} -f {input.ref} -s {wildcards.sample} - > {output.vcf}
		'''

rule merge_bam_deepvariant:
	input:
		bam = lambda wildcards: expand('{{sample}}/align/{{sample}}.movie_{movie}.mCG.ref_{{hap}}.bam',movie=dict_movie_bam[wildcards.sample].keys())
	output:
		outbam = temp('{sample}/align/{sample}.mCG.ref_{hap}.bam'),
		outbai = temp('{sample}/align/{sample}.mCG.ref_{hap}.bam.bai')
	resources:
		mem = 12,
		hrs = 70,
	threads: 4
	shell:
		'''
		source activate methylation_hifi
		samtools merge -@{threads} {output.outbam} {input.bam}
		samtools index {output.outbam}
		'''

rule run_deepvariant:
	input:
		fa = find_asm,
		bam = rules.merge_bam_deepvariant.output.outbam,
		bai = rules.merge_bam_deepvariant.output.outbai,

	output:
		vcf = '{sample}/deepvariant/{sample}.deepvar.ref_{hap}.vcf.gz',
		tbi = '{sample}/deepvariant/{sample}.deepvar.ref_{hap}.vcf.gz.tbi',
	resources:
		#mem = 16,
		#hrs = 72,
		load = 100,
	params:
		gpu = gpu_check(deepvar_image),
	threads: 50
	shell:
		'''
		module load singularity/3.6.0
		singularity exec {params.gpu} --bind /net/:/net {deepvar_image} run_deepvariant --model_type PACBIO --ref {input.fa} --reads {input.bam} --output_vcf {output.vcf} --num_shards {threads}
		'''

rule paf_call_het:
	input:
		deepvar_vcf = '{sample}/deepvariant/{sample}.deepvar.ref_{hap}.vcf.gz',
	output:
		tmp_vcf = '{sample}/deepvariant/{sample}.deepvar.ref_{hap}.het.vcf'
	resources:
		mem = 4,
		hrs = 8,
	threads: 1
	shell:
		'''
		python {SDIR}/script/het_deepvar.py {input.deepvar_vcf} {output.tmp_vcf}
		'''

rule filter_paf_call:
	input:
		het_vcf = '{sample}/deepvariant/{sample}.deepvar.ref_{hap}.het.vcf',
		paf_vcf = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.vcf',
	output:
		vcf = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.filtered.vcf'
	resources:
		mem = 4,
		hrs = 8,
	threads: 1
	shell:
		'''
		python {SDIR}/script/filter_paf_callset.py {input.het_vcf} {input.paf_vcf} {output.vcf}
		'''

rule make_mockphasedvcf:
	input:
		filtered_vcf = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.filtered.vcf',
	output:
		phased_vcf = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.phased.vcf.gz',
		phased_vcf_tbi = '{sample}/paf/{sample}.wga.{alt}_to_{hap}.phased.vcf.gz.tbi', 
	resources:
		mem = 8,
		hrs = 8,
	threads: 1
	shell:
		'''
		python {SDIR}/script/vcf2mockphasedvcf.py {input.filtered_vcf}
		bgzip -c {wildcards.sample}/paf/{wildcards.sample}.wga.{wildcards.alt}_to_{wildcards.hap}.phased.vcf > {output.phased_vcf}
		tabix -p vcf {output.phased_vcf}
		'''

rule samtools_reformat:
	input:
		aligned_bam = '{sample}/align/{sample}.movie_{movie}.mCG.ref_{hap}.bam',
	output:
		reheader_bam = temp('{sample}/align/{sample}.movie_{movie}.mCG.ref_{hap}.reheader.bam'),
		reheader_bai = '{sample}/align/{sample}.movie_{movie}.mCG.ref_{hap}.reheader.bam.bai',
	resources:
		mem = 12,
		hrs = 32,
	params:
		sed = r"s/SM:[^\t]*/SM:{sample}/g",
	threads: 4
	shell:
		'''
		samtools view -H {input.aligned_bam} | sed "{params.sed}" | samtools reheader - {input.aligned_bam} > {output.reheader_bam}
		samtools index {output.reheader_bam}
		'''

rule run_whatshap:
	input:
		ref = lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		vcf = whatshap_vcf,
		bam = '{sample}/align/{sample}.movie_{movie}.mCG.ref_{hap}.reheader.bam',
	output:
		haplotag_bam = temp('{sample}/align/{sample}.movie_{movie}.haplotagged.ref_{hap}.bam'),
		readid = temp('{sample}/align/readid/{sample}.movie_{movie}.all.ref_{hap}.readid.txt'),
	resources:
		mem = 12,
		hrs = 32,
	threads: 1
	shell:
		'''
		source activate methylation_hifi
		whatshap haplotag --ignore-read-groups --output-haplotag-list {output.readid} -o {output.haplotag_bam} --reference {input.ref} {input.vcf} {input.bam}
		'''

rule partitioning_hap_reads:
	input:
		all_reads = '{sample}/align/readid/{sample}.movie_{movie}.all.ref_{hap}.readid.txt',
	output:
		ref_reads = '{sample}/align/readid/{sample}.movie_{movie}.haplotype_{hap}.ref_{hap}.readid.txt',
	resources:
		mem = 4,
		hrs = 2,
	params:
		alt_reads = lambda wildcards: '{sample}/align/readid/{sample}.movie_{movie}.haplotype_{alt}.ref_{hap}.readid.txt'.format(sample=wildcards.sample,movie=wildcards.movie,hap=wildcards.hap,alt = alt_hap[wildcards.hap])
	threads: 1
	shell:
		'''
		cat {input.all_reads} | awk -v OFS='\\t' '{{if ($2 == "H2") print $1}}' > {params.alt_reads}
		cat {input.all_reads} | awk -v OFS='\\t' '{{if ($2 == "H1") print $1}}' > {output.ref_reads}
		'''

rule exclude_conflict_hap_reads:
	input:
		hap_ref_reads = '{sample}/align/readid/{sample}.movie_{movie}.haplotype_{hap}.ref_{hap}.readid.txt',
		hap_alt_reads = alt_hap_reads,
	output:
		curated_hap_reads = '{sample}/align/readid/{sample}.movie_{movie}.haplotype_{hap}.haplotype.readid.txt',
	resources:
		mem = 4,
		hrs = 2,
	threads: 1
	shell:
		'''
		comm -23 <(sort {input.hap_ref_reads}) <(sort {input.hap_alt_reads}) > {output.curated_hap_reads}
		'''

rule filter_bam_readid_haps:
	input:
		hap_reads = '{sample}/align/readid/{sample}.movie_{movie}.haplotype_{hap}.haplotype.readid.txt',
		haplo_bam   = '{sample}/align/{sample}.movie_{movie}.haplotagged.ref_{hap}.bam',
	output:
		haplo_hap_bam = temp('{sample}/align/{sample}.movie_{movie}.readid_{hap}.haplotagged.ref_{hap}.bam'),
	resources:
		mem = 12,
		hrs = 32,
	threads: 4
	shell:
		'''
		samtools view -N {input.hap_reads} -@{threads} -o {output.haplo_hap_bam} {input.haplo_bam}
		'''

rule merge_bam_haplotag:
	input:
		bam = lambda wildcards: expand('{{sample}}/align/{{sample}}.movie_{movie}.haplotagged.ref_{{hap}}.bam',movie=dict_movie_bam[wildcards.sample].keys())
	output:
		outbam = '{sample}/align/{sample}.haplotagged.ref_{hap}.bam',
		outbai = '{sample}/align/{sample}.haplotagged.ref_{hap}.bam.bai'
	resources:
		mem = 12,
		hrs = 30,
	threads: 4
	shell:
		'''
		source activate methylation_hifi
		samtools merge -@{threads} {output.outbam} {input.bam}
		samtools index {output.outbam}
		'''

rule pbsv_discover:
	input:
		repeat = tandem_repeat,
		bam = '{sample}/align/{sample}.haplotagged.ref_{hap}.bam',
	output:
		svsig = '{sample}/pbsv/{sample}.haplotagged.ref_{hap}.svsig.gz'
	resources:
		mem = 12,
		hrs = 50,
	threads: 1
	shell:
		'''
		source activate methylation_hifi
		pbsv discover --hifi -b {input.repeat} {input.bam} {output.svsig}
		'''

rule pbsv_call:
	input:
		ref_fa =  lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		svsig = rules.pbsv_discover.output.svsig,
	output:
		vcf = '{sample}/pbsv/{sample}.haplotagged.ref_{hap}.vcf',
	resources:
		mem = 4,
		hrs = 50,
	threads: 12
	shell:
		'''
		source activate methylation_hifi
		pbsv call -j {threads} {input.ref_fa} {input.svsig} {output.vcf}
		'''

rule run_pbcpg:
	input:
		ref_fa =  lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		bam = '{sample}/align/{sample}.haplotagged.ref_{hap}.bam',
	output:
		bed = '{sample}/pb_cpg/{sample}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',
		h1 = '{sample}/pb_cpg/{sample}.haplotagged.ref_{hap}_pbcpg.hap1.reference.bed',
		h2 = '{sample}/pb_cpg/{sample}.haplotagged.ref_{hap}_pbcpg.hap2.reference.bed',
		out = '{sample}/pb_cpg/{sample}.haplotagged.ref_{hap}_pbcpg-aligned_bam_to_cpg_scores.log',
	#resources:
	#	pbcpg = 50,
	#	mem = 24,
	#	hrs = 120,
	threads: 48
	shell:
		'''
		source activate methylation_hifi
		python /net/eichler/vol26/home/hsjeong/nobackups/programs/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b {input.bam} -f {input.ref_fa} -o {wildcards.sample}/pb_cpg/{wildcards.sample}.haplotagged.ref_{wildcards.hap}_pbcpg -m reference -t {threads}
		'''
rule rename_pbcpg:
	input:
		true_hap1 = pbcpg_true_hap1_check,
		true_hap2 = pbcpg_true_hap2_check,
		true_combined = '{sample}/pb_cpg/{sample}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',
	output:
		bed = '{sample}/pb_cpg/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',
		bed_h1 = '{sample}/pb_cpg/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.hap1.reference.bed',
		bed_h2 = '{sample}/pb_cpg/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.hap2.reference.bed',
	resources:
		mem = 4,
		hrs = 5,
	threads: 1
	shell:
		'''
		cp {input.true_hap1} {output.bed_h1}
		cp {input.true_hap2} {output.bed_h2}
		cp {input.true_combined} {output.bed}
		'''
rule run_pbcpg_denovo:
	input:
		ref_fa =  lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		bam = '{sample}/align/{sample}.haplotagged.ref_{hap}.bam',
	output:
		bed = '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_{hap}_pbcpg.combined.denovo.bed',
		h1 = '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_{hap}_pbcpg.hap1.denovo.bed',
		h2 = '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_{hap}_pbcpg.hap2.denovo.bed',
		out = '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_{hap}_pbcpg-aligned_bam_to_cpg_scores.log',
	resources:
		pbcpg = 50,
		mem = 32,
		hrs = 120,
	threads: 8
	shell:
		'''
		source activate methylation_hifi
		python /net/eichler/vol26/home/hsjeong/nobackups/programs/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b {input.bam} -f {input.ref_fa} -o {wildcards.sample}/pb_cpg_denovo/{wildcards.sample}.haplotagged.ref_{wildcards.hap}_pbcpg -t {threads}
		'''
rule rename_pbcpg_denovo:
	input:
		true_hap1 = pbcpg_denovo_true_hap1_check,
		true_hap2 = pbcpg_denovo_true_hap2_check,
		true_combined = '{sample}/pb_cpg_denovo/{sample}.haplotagged.ref_{hap}_pbcpg.combined.denovo.bed',
	output:
		bed = '{sample}/pb_cpg_denovo/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.combined.denovo.bed',
		bed_h1 = '{sample}/pb_cpg_denovo/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.hap1.denovo.bed',
		bed_h2 = '{sample}/pb_cpg_denovo/renamed/{sample}.haplotagged.ref_{hap}_pbcpg.hap2.denovo.bed',
	resources:
		mem = 4,
		hrs = 5,
	threads: 1
	shell:
		'''
		cp {input.true_hap1} {output.bed_h1}
		cp {input.true_hap2} {output.bed_h2}
		cp {input.true_combined} {output.bed}
		'''

rule merge_bam_haplotag_hap:
	input:
		bam = lambda wildcards: expand('{{sample}}/align/{{sample}}.movie_{movie}.readid_{{hap}}.haplotagged.ref_{{hap}}.bam',movie=dict_movie_bam[wildcards.sample].keys())
	output:
		outbam = '{sample}/align/{sample}.readid_{hap}.haplotagged.ref_{hap}.bam',
		outbai = '{sample}/align/{sample}.readid_{hap}.haplotagged.ref_{hap}.bam.bai'
	resources:
		mem = 12,
		hrs = 30,
	threads: 4
	shell:
		'''
		source activate methylation_hifi
		samtools merge -@{threads} {output.outbam} {input.bam}
		samtools index {output.outbam}
		'''

rule pbsv_discover_hap:
	input:
		repeat = tandem_repeat,
		bam = '{sample}/align/{sample}.readid_{hap}.haplotagged.ref_{hap}.bam',
	output:
		svsig = '{sample}/pbsv/{sample}.readid_{hap}.haplotagged.ref_{hap}.svsig.gz'
	resources:
		mem = 12,
		hrs = 50,
	threads: 1
	shell:
		'''
		source activate methylation_hifi
		pbsv discover --hifi -b {input.repeat} {input.bam} {output.svsig}
		'''

rule pbsv_call_hap:
	input:
		ref_fa =  lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		svsig = rules.pbsv_discover_hap.output.svsig,
	output:
		vcf = '{sample}/pbsv/{sample}.readid_{hap}.haplotagged.ref_{hap}.vcf',
	resources:
		mem = 4,
		hrs = 50,
	threads: 8
	shell:
		'''
		source activate methylation_hifi
		pbsv call -j {threads} {input.ref_fa} {input.svsig} {output.vcf}
		'''

rule run_pbcpg_hap:
	input:
		ref_fa =  lambda wildcards: sample_df.at[wildcards.sample, wildcards.hap],
		bam = '{sample}/align/{sample}.readid_{hap}.haplotagged.ref_{hap}.bam',
	output:
		bed = '{sample}/pb_cpg/{sample}.readid_{hap}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',
		out = '{sample}/pb_cpg/{sample}.readid_{hap}.haplotagged.ref_{hap}_pbcpg-aligned_bam_to_cpg_scores.log',
	#resources:
	#	mem = 16,
	#	hrs = 32,
	threads: 16
	shell:
		'''
		source activate methylation_hifi
		python /net/eichler/vol26/home/hsjeong/nobackups/programs/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b {input.bam} -f {input.ref_fa} -o {wildcards.sample}/pb_cpg/{wildcards.sample}.readid_{wildcards.hap}.haplotagged.ref_{wildcards.hap}_pbcpg -m reference -t {threads}
		'''
rule rename_pbcpg_hap:
	input:
		true_combined = '{sample}/pb_cpg/{sample}.readid_{hap}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',
	output:
		bed = '{sample}/pb_cpg/renamed/{sample}.readid_{hap}.haplotagged.ref_{hap}_pbcpg.combined.reference.bed',
	resources:
		mem = 4,
		hrs = 5,
	threads: 1
	shell:
		'''
		cp {input.true_combined} {output.bed}
		'''




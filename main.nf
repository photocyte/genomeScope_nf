nextflow.preview.dsl=2

process jellyfishCount {
cpus 10 
memory "50GB"
conda "jellyfish"
input:
 val theK
 path reads
output:
 path "reads.histo"
shell:
'''
##Decompress files in parallel using named pipes.
for f in !{reads}
do
 mkfifo ${f%.gz}
 zless $f > ${f%.gz} &
 echo ${f%.gz} >> fifos.txt
done

jellyfish count -C -m !{theK} -s 1000000000 -t !{task.cpus} $(cat fifos.txt | tr '\n' ' ') -o reads.jf
jellyfish histo -t !{task.cpus} reads.jf > reads.histo

##Cleanup temporary files. The named pipes won't take up any space.
rm -f reads.jf

'''
}

process genomeScope2 {
conda "r"
input:
 val theK
 path kmerCounts
shell:
'''
git clone https://github.com/tbenavi1/genomescope2.0.git
cd genomescope2.0/
sed -i "s^\"~/R_libs/\"^\"$CONDA_PREFIX/lib/R_libs\"^g" install.R
Rscript install.R
genomescope.R -i !{kmerCounts} -o output_dir -k !{theK}
'''
}

workflow countAndPlot_wf {
 take: kmerK ; reads
 main:
  jellyfishCount(kmerK,reads)
  genomeScope2(kmerK,jellyfishCount.out)
}

workflow {
 log.info """\
 G E N O M E S C O P E - 2 - N E X T F L O W
 =============================================
 reads        : ${params.reads}
 (If 'null', you have to pass a --reads parameter)
 """
 countAndPlot_wf(Channel.value("21"),Channel.fromPath(params.reads).collect())
}

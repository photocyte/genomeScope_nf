nextflow.preview.dsl=2

//// Process definitions ////

process jellyfishCount {
cpus 10 
memory "50GB"
conda "jellyfish"
input:
 val theK
 val maxKmerCov
 path reads
output:
 val theK
 val maxKmerCov
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
jellyfish histo -t !{task.cpus} --high !{maxKmerCov} reads.jf > reads.histo

##Cleanup temporary files. The named pipes won't take up any space.
rm -f reads.jf
'''
}

process genomeScope2 {
conda "r r-minpack.lm r-argparse"
input:
 val theK
 val maxKmerCov
 path kmerCounts
shell:
'''
git clone https://github.com/tbenavi1/genomescope2.0.git
cd genomescope2.0/

#mkdir "$CONDA_PREFIX/lib/R_libs"
#sed -i "s^\"~/R_libs/\"^\"$CONDA_PREFIX/lib/R_libs\"^g" install.R
#Rscript install.R
#echo "$(head -n 1 install.R)" >> $CONDA_PREFIX/lib/R/etc/Renviron

./genomescope.R -i !{kmerCounts} -m !{maxKmerCov} -o output_dir -k !{theK}
'''
}

//// Workflow definitions ////

workflow countAndPlot_wf {
 take: kmerK ; maxCov ; reads
 main:
  jellyfishCount(kmerK,maxCov,reads) | genomeScope2
}

workflow {
 log.info """\
  G E N O M E S C O P E - 2 - N E X T F L O W
 =============================================
 reads        : ${params.reads}
 (If 'null', you have to pass a --reads parameter)
 """
 theK = Channel.value("21")
 maxKmerCov = Channel.value("1000000")
 reads = Channel.fromPath(params.reads).collect()

 countAndPlot_wf(theK,maxKmerCov,reads)
}

nextflow.enable.dsl=2

//// Process definitions ////

process jellyfishCount {
cpus 10 
memory "100GB"
conda "jellyfish"
input:
 val theK
 val maxKmerCov
 path reads
output:
 path "*.jf", emit:jellyfishDbs
tag "${reads[0]}"
shell:
'''
##Decompress files in parallel using named pipes.
for f in !{reads}
do
 mkfifo ${f%.gz}
 zless $f > ${f%.gz} &
 echo ${f%.gz} >> fifos.txt
done

jellyfish count -C -m !{theK} -s 10000000000 -t !{task.cpus} --upper-count=!{maxKmerCov} $(cat fifos.txt | tr '\n' ' ') -o !{reads[0]}.jf

##Cleanup temporary files. The named pipes won't take up any space.
##rm -f reads.jf
'''
}

process jellyfishMerge {
executor "local"
cpus 1
memory "100GB"
conda "jellyfish"
input:
 val theK
 val maxKmerCov
 path jellyfishDbs
output:
 val theK
 val maxKmerCov
 path "merged.jf"
shell:
'''
jellyfish merge --upper-count=!{maxKmerCov} !{jellyfishDbs} -o merged.jf
'''
}

process jellyfishHisto {
cpus 10 
memory "100GB"
conda "jellyfish"
input:
 val theK
 val maxKmerCov
 path jellyfishDbs
output:
 val theK
 val maxKmerCov
 path "merged.histo"
shell:
'''
jellyfish histo -t !{task.cpus} --high=!{maxKmerCov} !{jellyfishDbs} -o merged.histo
'''
}

process genomeScope2 {
executor "local"
publishDir './results/' , mode:'link'
conda "r r-minpack.lm r-argparse"
input:
 val theK
 val maxKmerCov
 path kmerCounts
output:
 path "output_dir/*"
shell:
'''
git clone https://github.com/tbenavi1/genomescope2.0.git

if [[ ! -d "$CONDA_PREFIX/lib/R_libs" ]]
then
 mkdir "$CONDA_PREFIX/lib/R_libs"
fi

sed -i "s^\"~/R_libs/\"^\"$CONDA_PREFIX/lib/R_libs\"^g" ./genomescope2.0/install.R
sed -i "s^, lib=local_lib_path^^g" ./genomescope2.0/install.R
cat ./genomescope2.0/install.R | grep -v "minpack.lm" | grep -v "argparse" > tmp.R ##installs handled by bioconda
mv -f tmp.R ./genomescope2.0/install.R
cd ./genomescope2.0/
Rscript install.R
cd ../

##echo "$(head -n 1 install.R)" >> $CONDA_PREFIX/lib/R/etc/Renviron

echo "Now running GenomeScope 2.0"
./genomescope2.0/genomescope.R -i !{kmerCounts} -m !{maxKmerCov} -o output_dir -k !{theK}
'''
}

//// Workflow definitions ////

workflow countAndPlot_wf {
 take: kmerK ; maxCov ; reads
 main:
  jellyfishCount(kmerK,maxCov,reads)
  jellyfishMerge(kmerK,maxCov,jellyfishCount.out.jellyfishDbs.collect()) | jellyfishHisto | genomeScope2
  
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
 reads = Channel.fromPath(params.reads)

 countAndPlot_wf(theK,maxKmerCov,reads)
}

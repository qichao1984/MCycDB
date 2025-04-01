# MCycDB
Alternative links containing unzipped database files: https://www.alipan.com/s/Gs62Gxax7ii 

# MCycDB: a curated database for comprehensively profiling methane cycling processes of environmental microbiomes
Methane is a critical greenhouse gas with significant impacts on global change. Shotgun metagenome sequencing has greatly facilitated profiling methane cycle gene families and associated microbial taxa. However, there are many limitations for public orthology databases in accurately profiling methane cycle process, such as inefficient database searching, difficulties in distinguishing homologous genes, and low coverage of methane cycle genes and/or gene (sub)families. Therefore, it is necessary to build a manually curated database for accurately and rapidly profiling functions and taxonomy of methane cycling microbial communities in metagenome sequencing data.
The developed MCycDB contains 298 methane cycling gene families covering 10 methane metabolism pathways with 610,208 representative sequences, and 20,616 homologous orthology groups were also included to reduce false positive sequence assignments.

Four files are included in MCycDB:
1. MCycDB_2021.zip: fasta format representative sequences obtained by clustering curated sequences at 100% sequence identity. This file can be used for "BLAST" searching methane cycle genes in shotgun metagenomes.
2. id2gene.map: a mapping file that maps sequence IDs to gene names, only sequences belonging to methane cycle gene families are included. Sequences for methane cycle gene homologs are not included. This file is used to generate methane cycle profiles from BLAST-like results against MCycDB.
3. MCycDB_FunctionProfiler.PL: a perl script for functional profiling of methane cycle genes.
4. MCycDB_TaxonomyProfiler.PL: a perl script for taxonomical profiling of methane cycle microbial communities

<b>DOWNLOAD/INSTALLATION</b>
<p>git clone https://github.com/qichao1984/MCycDB.git</p>

<b>Dependencies and Tools</b>
<p>Perl modules that can be easily installed via cpan:</p>
<p>List::Util</p>
<p>Getopt::Long<p>

<b>Dependencies for MCycDB_FunctionProfiler.PL, currently supported database searching tools are:</b>
<p>usearch: https://www.drive5.com/usearch/download.html
<p>diamond: https://github.com/bbuchfink/diamond/releases
<p>blast: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz

<b>Dependencies for MCycDB_TaxonomyProfiler.PL:</b>
<p>seqtk: https://github.com/lh3/seqtk.git
<p>kraken2:  https://github.com/DerrickWood/kraken2.git

<b>USAGE</b>
Before getting started, please modify both scripts (MCycDB_FunctionProfiler.PL, MCycDB_TaxonomyProfiler.PL) at lines 6-18 to specify the locations of third party tools and their parameters. If the tools are already in the system path, no revision is needed. By default, basic parameters are used for these tools. Users are encouraged to make revisions in cases of short reads and/or expecting more strict/relaxed results. We also encourage users to develop useful implementations based on MCycDB.

Note: Kraken2 database could be downloaded from https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads, or built locally.

<b>Example for using MCycDB_FunctionProfiler.PL:</b>

perl MCycDB_FunctionProfiler.PL -d <workdir> -m <diamond|usearch|blast> -f <filetype> -s <seqtype> -si <sample size info file> -rs <random sampling size> -o <outfile>

<b>Detailed explanations:</b>

-d: specify the directory where your fasta/fastq (or gzipped) files are located.

-m: specify the database searching program you plan to use, currently diamond, usearch and blast are supported.

-f: specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s: sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have, note that no file extensions should be included here.

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

-o: the output file for methane cycle gene profiles.

<b>Example for using MCycDB_TaxonomyProfiler.PL:</b>

perl MCycDB_TaxonomyProfiler.PL -d <workdir> -m <diamond|usearch|blast> -f <filetype> -s <seqtype> -si <sample size info file> -rs <random sampling size>

<b>Detailed explanations:</b>

-d: specify the directory where your fasta/fastq (or gzipped) files are located.

-m: specify the database searching program you plan to use, currently diamond, usearch and blast are supported.

-f: specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s: sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have, note that no file extensions should be included here.

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

<b>Example for sampleinfo.txt:</b>

sample_name1 sequence_number

sample_name2 sequence_number

sample_name3 sequence_number

Reference tools:seqkit https://bioinf.shenwei.me/seqkit/usage/#seqkit

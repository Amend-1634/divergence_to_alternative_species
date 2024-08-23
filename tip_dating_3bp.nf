// Enable DSL2 functionality
nextflow.enable.dsl = 2
include { REMAP_ANC_READS } from './remapping_ancient_reads_2.nf'

// Workflow definition
workflow {

	// 
	Channel.fromPath(params.all_input).splitCsv(sep: "\n").collect().set{ all_files }
        //.map{ it -> it.trim() }
        //.filter{ it != null }

	//all_files.view()

	//first line as ancient bam file
	all_files.map{ it[0] }.set{ anc_bam }

	//second line as modern reference fasta file
	all_files.map{ it[1] }.set{ mod_fa }

	//since the third line, it's all modern samples other than the modern reference for the bam
	all_files.subList(2, all_files.size()).set{ other_mod_fa }

	//n_files.view()

	// turn the ancient sample bam into fastq file for later remapping to other modern accessions
	BAM_TO_FASTQ( anc_bam )

	// remap the reads in the bam file to other modern references
	REMAP_ANC_READS( BAM_TO_FASTQ.out, other_mod_fa)

	// output the mapped bam and its corresponding modern reference
	//REMAP_ANC_READS.out.view()

	//collect all bam file

        //(1)ancient reads remmaped to be bam, with its corresponding other modern fasta file
        anc_bam_other_mod_fa = REMAP_ANC_READS.out
        //anc_bam_other_mod_fa.view()

	//(2)ancient bam with the modern reference
	anc_bam_mod_fa = anc_bam.merge(mod_fa)
	//anc_bam_mod_fa.view()

	paired_bam_fa = anc_bam_mod_fa.concat( anc_bam_other_mod_fa )
	paired_bam_fa.view()

	BAM_TO_FASTA( paired_bam_fa )
	//BAM_TO_FASTA.out.view()

	COMPUTE_MUTATION( BAM_TO_FASTA.out )
	//COMPUTE_MUTATION.out.view()

	SUMMARY_MUTATION(  COMPUTE_MUTATION.out.collect(), params.label )
	SUMMARY_MUTATION.out.view()
}

process BAM_TO_FASTQ {
	cpus params.threads

        input:
                path(bam)

        output:
                path("*.fastq")

        script:
        """
	samtools bam2fq $bam > \$(echo $bam | sed 's/.bam//').fastq
        """

}

process BAM_TO_FASTA {

    //publishDir params.results_dir, mode: "copy"

    input:
        tuple path(bam), path(fasta)

    output:
        path("*.fa")
	path("cp.*") // the copied corresponding modern fasta

    script:
    """ 
	#sort the file
	bam_name=\$(echo $bam | sed 's/.bam//')
	sorted_bam=\${bam_name}_sorted.bam
	samtools sort -o \$sorted_bam \${bam_name}.bam

	#(1)turn the ancient reads in bam into a fasta file
        fa_name=\$(echo $bam | sed 's/.bam//').anc

        #extract fasta file from a bam--using angsd
        angsd -i \$sorted_bam -doFasta 1 -doCounts 1 -out \$fa_name
                # -doFasta      Generate a fasta for a BAM file
                # -doCounts     Calculate various counts statistics
        pigz -d *.fa.gz

	#(2)also output the corresponding modern reference fasta file
	cp $fasta cp.\$(basename $fasta)
    """
}

process COMPUTE_MUTATION {

    publishDir params.results_dir, mode: "copy"

    input:
        path(anc_fa)
	path(mod_fa)
    output:
        path("*.stat")

    script:
    """ 
        #estimate divergence of this contig
        rates=\$(est_mutations_3bp_transition.py $anc_fa $mod_fa)

        mutation=\$(echo \$rates | cut -f 1 -d' ')
        transition=\$(echo \$rates | cut -f 2 -d' ')
        transversion=\$(echo \$rates | cut -f 3 -d' ')
        total_covered_sites=\$(echo \$rates | cut -f 4 -d' ')

        echo -e "$mod_fa\t\$mutation\t\$transition\t\$transversion\t\$total_covered_sites" > \$(basename ${anc_fa} | sed 's/.fa//' | sed 's/.fna//' | sed 's/.fasta//')-\$(basename ${mod_fa} | sed 's/.fa//' | sed 's/.fna//' | sed 's/.fasta//').stat

    """
}

process SUMMARY_MUTATION {

        publishDir params.results_dir, mode: "copy"

        cpus params.threads

        input:
                path(stat)
                val(label)

        output:
                path("*.all.stat")

        script:
        """
        #combine all stat
        cat *.stat >> ${label}.all.stat

        """
}

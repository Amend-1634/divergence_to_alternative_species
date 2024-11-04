// Enable DSL2 functionality
nextflow.enable.dsl = 2

// Workflow definition
workflow REMAP_ANC_READS{

	take:
		anc_fastq
		mod_fasta

	main:

		INDEXING( mod_fasta, params.label, params.threads )

		MAPPING( INDEXING.out, anc_fastq )

		REMOVE_C_TO_T_FORWARD_G_TO_A_REVERSE( MAPPING.out )

	emit:
		REMOVE_C_TO_T_FORWARD_G_TO_A_REVERSE.out // pair the bam with its corresponding fasta file of other modern samples

}

process INDEXING {
        cpus params.threads

        input:
                each(fasta)
                val(label)
                val(threads)

        output:
                tuple path("*.bt2l"), path("cp.*")

        script:
        """
        #build the mapping index from the concacenated fasta #direct the error to output
        bowtie2-build --large-index --threads $threads $fasta $label

	#output the fasta (other modern reference) again
	cp $fasta cp.\$(basename $fasta)
        """
}

process MAPPING {
        //publishDir params.mapping, mode: "copy"

        input:
                tuple path(index), path(other_mod_fa)
                path(reads)

        output:
                tuple path("*.bam"), path("cp.cp.*")
		

        script:
        """
        #obtain the name of the index files
        ind_name=\$(ls *.bt2* | sed -n 1p | cut -f 1 -d'.' )

        bowtie2 --very-sensitive -p ${params.threads} -x \$ind_name -U $reads | \
  samtools view -@ ${params.threads} -Sb -q 1 - > \$(basename $reads)_\$(basename $other_mod_fa).bam


        #output the corresponding fasta file
	cp $other_mod_fa cp.\$(basename $other_mod_fa)
        """
}

process REMOVE_C_TO_T_FORWARD_G_TO_A_REVERSE {

	//publishDir params.mapping, mode: "copy"
	
	container 'docker://quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0'

        input:
                tuple path(bam), path(fa)

        output:
                tuple path("*.rm_damage.bam"), path("cp.cp.cp.*")


        script:
	"""
	# Sort the BAM file
	samtools sort -o \$(basename $bam .bam).sorted.bam $bam

	# Index the sorted BAM file
	samtools index -c \$(basename $bam .bam).sorted.bam

	# Run the Python script on the sorted and indexed BAM file
	adna_sslib_damage_removal.py -b \$(basename $bam .bam).sorted.bam -o \$(basename $bam .bam).rm_damage.bam

        #make sure the fasta file is output with cp.cp.cp.
	if [[ "$fa" != cp.cp.* ]]; then
	    # Add "cp.cp." as the prefix
	    cp $fa cp.cp.$fa
	    fa="cp.cp.$fa"
	else
	    fa="$fa"
	fi

        cp \$fa cp.\$(basename \$fa)
        """


}


process BAM_TO_FASTA {

    publishDir params.mapping, mode: "copy"

    input:
        path(bam)

    output:
        path("*.fa")

    script:
    """
        fa_name=\$(echo $bam | sed 's/.bam//').anc
        
        #extract fasta file from a bam--using angsd
        angsd -i $bam -doFasta 1 -doCounts 1 -out \$fa_name
                # -doFasta      Generate a fasta for a BAM file
                # -doCounts     Calculate various counts statistics
        pigz -d *.fa.gz
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
        rates=\$(est_mutations2.py $anc_fa $mod_fa)

        mutation=\$(echo \$rates | cut -f 1 -d' ')
        transition=\$(echo \$rates | cut -f 2 -d' ')
        transversion=\$(echo \$rates | cut -f 3 -d' ')
        total_covered_sites=\$(echo \$rates | cut -f 4 -d' ')

        echo -e "$mod_fa\t\$mutation\t\$transition\t\$transversion\t\$total_covered_sites" >> \$(basename ${anc_fa} | sed 's/.fa//' | sed 's/.fna//' | sed 's/.fasta//')-\$(basename ${mod_fa} | sed 's/.fa//' | sed 's/.fna//' | sed 's/.fasta//').stat

    """
}

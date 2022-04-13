#!/bin/bash


#
# This is primary bioinformatics pipeline intended for processing, mapping, assembling, and annotation of bulk Nanopore RNA-seq data
#


echo -e ""

echo -e "** Wellcome to the transcripts assembling pipeline!\n"

echo -e "/** This pipeline is intended for processing, aligning, transcript assembling, and annotation of Nanopore RNA-seq data! You will make use of tools like Pychopper, minimap2, Samtools, StringTie, gffcompare and gffread to obtain assembled and annotated full-length RNA transcripts. **/\n"


read -p "Would you like to proceed: (yes/no)? " USER_ANSWER



if [[ $USER_ANSWER = yes ]]

then

    echo -e "\n/** In order to go through the pipeline from the very beginning to the end successfully, you will have to specify 3 paths: PATH to single fastq file containing all reads of interest (cDNA or dRNA); PATH to reference genome in FASTA format; PATH to reference annotation in GTF/GFF format. **/\n"

else

    echo -e "\nBye!"
    exit 0

fi


echo -e "[Step 1/6] The first step in this pipeline is trimming of sequencing adapters if you analyse cDNA. You have to skip this step if your data comes from direct RNA-seq experiment or if you already have your adapters trimmed.\n"

read -p "Would you like to run Pychopper to trim your adapters: (yes/no)? " RUNNING_PYCHOPPER


if [[ $RUNNING_PYCHOPPER = yes ]]

then

    echo -e "\n/** You chose to run Pychopper to trim and orient full-length Nanopore cDNA reads. **/\n"

    read -p "Now you have to specify full path to the directory that contains the input fastq file. NOTE: you must provide single fastq file! Path to file: " PYCHOPPER_INPUT_FASTQ

    if [[ $PYCHOPPER_INPUT_FASTQ = "" ]]

    then

        echo -e  "\n* Invalid input! You have to supply one positional argument - path to input fastq file. Either 0, or more than 1 argumnets have been passed."

        exit 1

    else
	
        echo -e ""

        mkdir -p pychopper_results

        PYCHOPPER_OUT_DIR="./pychopper_results"
        
        PYCHOPPER_OUT_FILE="./pychopper_results/trimmed_full_length_cdna_reads.fq"

        echo -e "* Running Pychopper...\n"

        cdna_classifier.py -t 8 -r "$PYCHOPPER_OUT_DIR"/report.pdf \
            -A "$PYCHOPPER_OUT_DIR"/aln_hits.bed \
            -S "$PYCHOPPER_OUT_DIR"/statistics.tsv \
            -u "$PYCHOPPER_OUT_DIR"/unclassified.fq \
            -w "$PYCHOPPER_OUT_DIR"/rescued.fq \
	    $PYCHOPPER_INPUT_FASTQ $PYCHOPPER_OUT_FILE && \ 
	    cecho success "✔ Pychopper run successfully!" || cecho error "Pychopper failed to run!"

    fi

fi


if [[ $RUNNING_PYCHOPPER = no ]]

then

    echo -e  "\n* Skipping running pychopper."

fi


echo -e "\n[Step 2/6] Step 2 in this pipeline is building reference genome index with minimap2.\n"

echo -e "/** You can skip this step if you already have indexed genome. It's recommended building index as it will speed up the downstream analysis. **/\n"


read -p "Would you like to index your reference genome: (yes/no)? " BUILDING_GENOME_INDEX


if [[ $BUILDING_GENOME_INDEX = yes ]] 

then 

    echo -e ""

    read -p "Enter path to reference genome: " PATH_TO_REFERENCE_GENOME

    if [[ $PATH_TO_REFERENCE_GENOME = "" ]]

    then 

        echo -e "\nInvalid input! You have not specified path to reference genome.\n"

	exit 2

    else 
	
        mkdir -p minimap_genome_index
     	
        MINIMAP_INDEX="./minimap_genome_index/genome_index.mmi"
	
        THREADS=8
	
        echo -e "\n* Building minimap index...\n"

        minimap2 -t $THREADS -k14 -I 1000G -d $MINIMAP_INDEX "$PATH_TO_REFERENCE_GENOME" && \
		cecho success "✔ Build index!" || cecho error "Minimap2 failed to build index!"

    fi

else

    echo -e "\n* Skipping building genome index."

fi


echo -e "\n[Step 3/6] The third step in this pipeline is aligning your reads to reference (indexed) genome with minimap2.\n"

echo -e "/** You can skip this step if you alredy have aligned and sorted reads. **/\n"


read -p "Would you like to align your reads with minimap2 in splice-aware mode: (yes/no)? " ALIGN_READS


if [[ $ALIGN_READS = yes ]] 

then

    echo -e ""

    read -p "Have you run pychopper in the first stage: (yes/no)? " RUN_PYCHOPPER

    echo -e ""

    read -p "Have you built genome index: (yes/no)? " BUILD_GEN_INDEX

    
    if [[ $RUN_PYCHOPPER = yes ]] && [[ $BUILD_GEN_INDEX = yes ]]
    
    then
	
        echo -e ""

        read -p "How many cores do you want to use for the alignment: " NUM_CORES

        mkdir -p alignment_dir
	
        echo -e "\n* Aligning reads...\n"

        minimap2 -t $NUM_CORES -ax splice -uf $MINIMAP_INDEX \
	        "./pychopper_results/trimmed_full_length_cdna_reads.fq" > ./alignment_dir/aligned_reads.sam \
	       	&& cecho success "✔ Reads aligned!" || cecho error "Minimap2 failed to align reads!" && \
		samtools view -q 40 -F2304 -S -b \
                ./alignment_dir/aligned_reads.sam > ./alignment_dir/aligned_reads.bam && \
                samtools sort -o ./alignment_dir/sorted_aligned_reads.bam ./alignment_dir/aligned_reads.bam && \
                samtools index ./alignment_dir/sorted_aligned_reads.bam && \
                cecho success "✔ Reads sorted and indexed!" || \
                cecho error "Samtools failed to sort reads!"


    elif [[ $RUN_PYCHOPPER = yes ]]  && [[ $BUILD_GEN_INDEX = no ]]

    then 
	
        echo -e ""
	
        read -p "How many cores do you want to use for the alignment: " NUM_THREADS
        
        echo -e ""        

        read -p "Enter path to indexed genome: " INDEXED_GENOME_PATH
	
        mkdir -p alignment_dir

        echo -e "\n* Aligning reads...\n"

        minimap2 -t $NUM_THREADS -ax splice -uf $INDEXED_GENOME_PATH \
                "./pychopper_results/trimmed_full_length_cdna_reads.fq" > ./alignment_dir/aligned_reads.sam && \
                cecho success "✔ Reads aligned!" || cecho error "Minimap2 failed to align reads!" && \
		samtools view -q 40 -F2304 -S -b \
                ./alignment_dir/aligned_reads.sam > ./alignment_dir/aligned_reads.bam && \
                samtools sort -o ./alignment_dir/sorted_aligned_reads.bam ./alignment_dir/aligned_reads.bam && \
                samtools index ./alignment_dir/sorted_aligned_reads.bam && \
                cecho success "✔ Reads sorted and indexed!" || \
                cecho error "Samtools failed to sort reads!"
    

    elif [[ $RUN_PYCHOPPER = no ]]  && [[ $BUILD_GEN_INDEX = yes ]]

    then

        echo -e ""

        read -p "How many cores do you want to use for the alignment: " NUM_THREADS
	
        echo -e "" 

        read -p "Enter path to fastq file that will be aligned: " READS

        mkdir -p alignment_dir

        echo -e "\n* Aligning reads...\n"

        minimap2 -t $NUM_THREADS -ax splice -uf $MINIMAP_INDEX \
                $READS > ./alignment_dir/aligned_reads.sam && \
                cecho success "✔ Reads aligned!" || cecho error "Minimap2 failed to align reads!" && \
                samtools view -q 20 -F2304 -S -b \
                ./alignment_dir/aligned_reads.sam > ./alignment_dir/aligned_reads.bam && \
                samtools sort -o ./alignment_dir/sorted_aligned_reads.bam ./alignment_dir/aligned_reads.bam && \
                samtools index ./alignment_dir/sorted_aligned_reads.bam && \
                cecho success "✔ Reads sorted and indexed!" || \
                cecho error "Samtools failed to sort reads!"


    elif [[ $RUN_PYCHOPPER = no ]]  && [[ $BUILD_GEN_INDEX = no ]]

    then
	
        echo -e ""

        read -p "How many cores do you want to use for the alignment: " NUMBER_THREADS
	
        echo -e "" 

        read -p "Enter path to fastq file that will be aligned: " READS_TO_ALIGN
	
        echo -e "" 

        read -p "Enter path to indexed genome: " PATH_INDEXED_GENOME
        
        mkdir -p alignment_dir        

        echo -e "\n* Aligning reads...\n"

        minimap2 -t $NUMBER_THREADS -ax splice -uf $PATH_INDEXED_GENOME \
                $READS_TO_ALIGN > ./alignment_dir/aligned_reads.sam && \
                cecho success "✔ Reads aligned!" || \
                cecho error "Minimap2 failed to align reads!" && \
		samtools view -q 40 -F2304 -S -b \
	       	./alignment_dir/aligned_reads.sam > ./alignment_dir/aligned_reads.bam && \
		samtools sort -o ./alignment_dir/sorted_aligned_reads.bam ./alignment_dir/aligned_reads.bam && \
		samtools index ./alignment_dir/sorted_aligned_reads.bam && \
	       	cecho success "✔ Reads sorted and indexed!" || \
		cecho error "Samtools failed to sort reads!"
    
    fi

else 

    echo -e "\n* Skipping aligning reads."

fi


echo -e "\n[Step 4/6] Now you can run StringTie to assemble the RNA-Seq alignments obtainded in the previous step into potential transcripts.\n"

read -p "Would you like to run StringTie to assemble RNA transcripts: (yes/no)? " RUN_STRINGTIE


if [[ $RUN_STRINGTIE = yes ]]

then 

    echo -e "\n/** StringTie takes as input a binary SAM (BAM) file sorted by reference position. This means that if run minimap2 in the previous step you already have this file and you can run StringTie directly. Otherwise, you will have to specify path to sorted BAM file. **/\n"

    read -p "Have you run minimap2 in the previous step: (yes/no)? " RUN_MINIMAP_AND_SAMTOOLS

    if [[ $RUN_MINIMAP_AND_SAMTOOLS = yes ]] 

    then 

        echo -e ""

        read -p "Enter path to reference annotation to use for guiding the assembly process (GTF/GFF): " PATH_TO_ANNOTATION
	
        echo -e ""

        read -p "On how many cores do you want to run StringTie: " NUM_CORES_STRINGTIE
	
        mkdir -p ./stringtie_results

        STRINGTIE_OUT_FILE="./stringtie_results/stringtie_long_reads.out.gff"

        STRINGTIE_IN_FILE="./alignment_dir/sorted_aligned_reads.bam"

        echo -e "\n* Running StringTie...\n"

        stringtie --rf -G $PATH_TO_ANNOTATION -L -p $NUM_CORES_STRINGTIE \
		--conservative -o $STRINGTIE_OUT_FILE $STRINGTIE_IN_FILE && \
		cecho success "✔ Trnascripts assembled!" || cecho error "StringTie failed to assemble transcripts!"

    else 
	
        echo -e ""

        read -p "Then enter path to sorted BAM file: " PATH_TO_BAM_FILE

        echo -e ""

        read -p "Enter path to reference annotation to use for guiding the assembly process (GTF/GFF): " PATH_TO_ANNOTATION

        echo -e ""

        read -p "On how many cores do you want to run StringTie: " NUM_CORES_STRINGTIE
	
        mkdir -p ./stringtie_results

        STRINGTIE_OUT_FILE="./stringtie_results/stringtie_long_reads.out.gff"

        echo -e "\n* Running StringTie...\n"

        stringtie --rf -G $PATH_TO_ANNOTATION -L -p $NUM_CORES_STRINGTIE \
		--conservative -o $STRINGTIE_OUT_FILE $PATH_TO_BAM_FILE && \
		cecho success "✔ Trnascripts assembled!" || cecho error "StringTie failed to assemble transcripts!"

    fi

else

    echo -e "\n* Skipping running StringTie."

fi


echo -e "\n[Step 5/6] The assembled transcripts generated in the previous step can be annotated by using gffcompare. One potential result from this analysis can be the discovery of novel transcripts not present in the reference database.\n"

read -p "Would you like to annotate the assembled transcriptome: (yes/no)? " ANNOTATE_TRANSCRIPTOME


if [[ $ANNOTATE_TRANSCRIPTOME = yes ]] 

then 

    echo -e "\n/** NOTE: if you skipped running StringTie in the previous step you will have to provide path to GTF or GFF file generated from StringTie. **/\n"

    read -p "Have you run StringTie in the previous step: (yes/no)? " RUN_STRINGTIE_IN_STEP_4

    if [[ $RUN_STRINGTIE_IN_STEP_4 = yes ]]

    then 

        mkdir -p ./gffcompare_results

        OUT_FILE_GFFCOMPARE="./gffcompare_results/stringtie"

        echo -e "\n* Running gffcompare...\n"

        gffcompare -o $OUT_FILE_GFFCOMPARE -r $PATH_TO_ANNOTATION -R $STRINGTIE_OUT_FILE && \
		cecho success "✓ Completed gffcompare"  || cecho error "✗ Failed gffcompare"


    else

	mkdir -p ./gffcompare_results

        echo -e "\n/** You will have to enter path to reference annotation in GFF or GTF format as well as path to GTF or GFF file generated from Stringtie. **/\n"
	
        read -p "Enter path to reference annotation: " PATH_TO_REF_ANNOTATION

        echo -e ""

        read -p "Enter path to GFF or GTF file generated by StringTie: " PATH_TO_STRINGTIE_GFF_FILE

        OUT_FILE_GFFCOMPARE="./gffcompare_results/stringtie"

        echo -e "\n* Running gffcompare...\n"

        gffcompare -o $OUT_FILE_GFFCOMPARE -r $PATH_TO_REF_ANNOTATION -R $PATH_TO_STRINGTIE_GFF_FILE && \
		cecho success "✓ Completed gffcompare"  || cecho error "✗ Failed gffcompare"

    fi

else 

    echo -e "\n* Skipping running gffcompare."

fi


echo -e "\n[Step 6/6] The last step in this pipeline is generation of reference transcriptome in FASTA format. For that purpose, we can use gffread to combine the assembled and annotated transcripts with genomic reference sequence.\n"

read -p "Would you like to create reference transcriptome: (yes/no)? " GEN_REF_TRANSCRIPTOME


if [[ $GEN_REF_TRANSCRIPTOME = yes ]]

then 

    echo -e "\n/** To create reference transcriptome with gffread you need reference genome in FASTA format and GFF/GTF file obtainded from StringTie and annotated with gffcompare. **/\n"
	
    read -p "Have you built genome index in step 2: (yes/no)? " HAS_GENOME_INDEX

    if [[ $HAS_GENOME_INDEX = yes ]]

    then 
	
	echo -e ""

	read -p "Have you run gffcompare in step 5: (yes/no)? " RUN_GFFCOMPARE_IN_STEP_5

	
	if [[ $RUN_GFFCOMPARE_IN_STEP_5 = yes ]] 

	then

	    mkdir -p gffread_results 
	    
	    GFFREAD_OUT="./gffread_results/gffread_transcriptome.fasta"

	    GFFREAD_IN="./gffcompare_results/stringtie.annotated.gtf"

	    echo -e "\n* Running gffread...\n"

	    gffread -F -g $PATH_TO_REFERENCE_GENOME -w $GFFREAD_OUT $GFFREAD_IN && \
		    cecho success "✓ Reference transcriptome created!" || \
		    cecho error "✗ Gffread failed to create reference transcriptome!"

	else 
	
	    echo -e ""

	    read -p "Enter path to GTF/GFF file generated from StringTie and annotated by gffcompare: " PATH_TO_STRINGTIE_GFF

	    mkdir -p gffread_results

        GFFREAD_OUT="./gffread_results/stringtie_transcriptome.fasta"

        echo -e "\n* Running gffread...\n"

        gffread -F -g $PATH_TO_REFERENCE_GENOME -w $GFFREAD_OUT $PATH_TO_STRINGTIE_GFF && \
            cecho success "✓ Reference transcriptome created!" || \
		    cecho error "✗ Gffread failed to create reference transcriptome!"

	fi

    else 

        echo -e ""

        read -p "Have you run gffcompare in step 5: (yes/no)? " RUN_GFFCOMPARE_IN_STEP_5


	if [[ $RUN_GFFCOMPARE_IN_STEP_5 = yes ]]

	then
	
        echo -e ""

        read -p "Enter path to genome reference in FASTA format: " PATH_TO_GENOME_REF

	    mkdir -p gffread_results

        GFFREAD_OUT="./gffread_results/stringtie_transcriptome.fasta"

	    GFFREAD_IN="./gffcompare_results/stringtie.annotated.gtf"

        echo -e "\n* Running gffread...\n"
		
        gffread -F -g $PATH_TO_GENOME_REF -w $GFFREAD_OUT $GFFREAD_IN && \
            cecho success "✓ Reference transcriptome created!"  || \
		    cecho error "✗ Gffread failed to create reference transcriptome!"
	    
	else 

	    echo -e ""

	    read -p "Enter path to GTF/GFF file generated from StringTie and annotated by gffcompare: " PATH_TO_STRINGTIE_GTF_FILE

	    echo -e ""

	    read -p "Enter path to genome reference in FASTA format: " PATH_TO_GENOME_REF

	    mkdir -p gffread_results

	    GFFREAD_OUT="./gffread_results/stringtie_transcriptome.fasta"

	    echo -e "\n* Running gffread...\n"

	    gffread -F -g $PATH_TO_GENOME_REF -w $GFFREAD_OUT $PATH_TO_STRINGTIE_GTF_FILE && \
            cecho success "✓ Reference transcriptome created!"  || \
		    cecho error "✗ Gffread failed to create reference transcriptome!"

    	fi

    fi
    
else 

    echo -e "\n* Skipping running gffread."

fi


echo -e "\nWell done! You went through the transcripts assembling pipeline successsfully!\n"

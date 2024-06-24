version 1.0

workflow CDR_Detection_Workflow {
	input {
		File input_fasta
		File censat_bed
		File mC_bam
		String sample_id
	}

	call extract_h1l {
		input:
			censat_bed = censat_bed,
			sample_id = sample_id
	}

	call samtools_index {
		input:
			input_bam = mC_bam,
			sample_id = sample_id
	}

	call modkit_pileup {
		input:
			input_fasta = input_fasta,
			indexed_bam = samtools_index.bam_file,
			index_bai = samtools_index.bam_index,
			sample_id = sample_id
	}

	scatter (contig in extract_h1l.contig_names) {
		call contig_extract_and_h1l_intersect {
			input:
				modkit_pileup_bed = modkit_pileup.modkit_pileup_bed,
				censat_h1l_bed = extract_h1l.censat_h1l_bed,
				contig_name = contig,
				sample_id = sample_id
		}
  
		call strict_detect {
			input:
				pileup_bed = contig_extract_and_h1l_intersect.pileup_h1l_intersect_bed,
				censat_h1l_bed = extract_h1l.censat_h1l_bed,
				contig_name = contig,
				sample_id = sample_id
		}
		
		call hmm_detect {
			input:
				pileup_bed = contig_extract_and_h1l_intersect.pileup_h1l_intersect_bed,
				strict_bed = strict_detect.strict_cdr_bed,
				contig_name = contig,
				sample_id = sample_id
		}
		
		call validate_cdrs {
			input:
				pileup_bed = contig_extract_and_h1l_intersect.pileup_h1l_intersect_bed,
				strict_bed = strict_detect.strict_cdr_bed,
				hmm_bed = hmm_detect.viterbi_cdr_bed,
				emissionMatrix = hmm_detect.emissionMatrix,
				transitionMatrix = hmm_detect.transitionMatrix,
				contig_name = contig,
				sample_id = sample_id
		}
	}
	
	call gather_outputs {
		input:
			hmm_cdrs = hmm_detect.viterbi_cdr_bed,
			pileup_beds = contig_extract_and_h1l_intersect.pileup_h1l_intersect_bed,
			strict_cdrs_beds = strict_detect.strict_cdr_bed,
			pileup_histograms = validate_cdrs.histogram,
			emission_heatmaps = validate_cdrs.emission_heatmap,
			transition_heatmaps = validate_cdrs.transition_heatmap,
			emissionBoundaries = hmm_detect.emissionBoundaries,
			emissionMatrix = hmm_detect.emissionMatrix,
			transitionMatrix = hmm_detect.transitionMatrix,
			sample_id = sample_id
	}

	output {
		File hmm_cdrs_bed = gather_outputs.hmm_cdrs_bed
		File pileup_bed = modkit_pileup.modkit_pileup_bed
		
		File igv_zip = gather_outputs.igv_zip
		File summary_zip = gather_outputs.summary_zip	
	}

	meta {
		author: "Julian Menendez"
		email: "jmmenend@ucsc.edu"
		description: "Automated CDR Detection, using windowing then an HMM."
	}
}

task extract_h1l {
	input {
		File censat_bed
		String sample_id

		Int memSizeGB = 16
		Int threadCount = 4
		Int preempts = 1
	}

	# Estimate disk size required
	Int final_disk_dize = 16

	# set outputs as wdl variables
	String censat_h1l_bed_output = "~{sample_id}.CenSat.H1L.bed"
	String contig_names_output = "contig_names.json"

	command <<<
		set -eux -o pipefail

		## extract only rows containing H1L(alpha-sat) in the fourth column
		awk '$4 ~ /H1L/' ~{censat_bed} > ~{censat_h1l_bed_output}

		## Declare an empty array
		declare -a contig_array=()

		## Use awk to extract unique entries from the first column and add them to the array
		mapfile -t contig_array < <(cut -f1 "~{censat_h1l_bed_output}" | sort | uniq)

		## save the contig names to a json file
		printf '%s\n' "${contig_array[@]}" | jq -R . | jq -s . > ~{contig_names_output}
	>>>

	output {
		File censat_h1l_bed = censat_h1l_bed_output
		Array[String] contig_names = read_json(contig_names_output)
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/jq:v1.0" 
		preemptible: preempts
	}
}

task samtools_index {
	input {
		File input_bam
		String sample_id

		Int memSizeGB = 256
		Int threadCount = 32
		Int preempts = 1
	}

	# Estimate disk size required
	Int input_bam_size = ceil(size(input_bam, "GB"))       
	Int final_disk_dize = input_bam_size * 2

	String bam_output = "~{sample_id}.bam"
	String index_output = "~{sample_id}.bam.bai"

	command <<<
		set -eux -o pipefail

		samtools sort -@ ~{threadCount} \
			-o ~{bam_output} \
			~{input_bam}
		samtools index -@ ~{threadCount} \
			-o ~{index_output} \
			~{bam_output} \
	>>>

	output {
		File bam_file = bam_output
		File bam_index = index_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/samtools:1.20"
		preemptible: preempts
	}
}

task modkit_pileup {
	input {
		File input_fasta 
		File indexed_bam
		File index_bai
		String sample_id

		Int memSizeGB = 256
		Int threadCount = 32
		Int preempts = 1
	}

	# Estimate disk size required
	Int input_fasta_size = ceil(size(input_fasta, "GB"))
	Int input_bam_size = ceil(size(indexed_bam, "GB"))       
	Int final_disk_dize = input_fasta_size + input_bam_size * 8

	String modkit_pileup_bed_output = "~{sample_id}.5mCpileup.bed"

	command <<<
		set -eux -o pipefail

		## run modkit pileup!
		modkit pileup \
			~{indexed_bam} ~{modkit_pileup_bed_output} \
			-t ~{threadCount} \
			--filter-percentile 0.66 \
			--ignore h \
			--force-allow-implicit \
			--cpg \
			--ref ~{input_fasta} \
			--combine-strands
	>>>

	output {
		File modkit_pileup_bed = modkit_pileup_bed_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/modkit:0.3.0"
		preemptible: preempts
	}
}

task contig_extract_and_h1l_intersect {
	input {
		File modkit_pileup_bed
		File censat_h1l_bed
		String contig_name
		String sample_id

		Int memSizeGB = 64
		Int threadCount = 4
		Int preempts = 1
	}

	# Estimate disk size required
	Int input_bed_size = ceil(size(modkit_pileup_bed, "GB"))    
	Int final_disk_dize = input_bed_size * 16

	String pileup_h1l_intersect_bed_output = "~{contig_name}_~{sample_id}_5mCpileup.H1L.bed"

	command <<<
		set -eux -o pipefail

		## grabs the 4 important columns for CDR predictions, and bedtools intersects with alpha-sat array
		awk -v contig="~{contig_name}" 'BEGIN {OFS="\t"} ($1 ~ contig && $4 ~ /m/) {print $1, $2, $3, $11}' ~{modkit_pileup_bed} | \
			bedtools intersect -a - -b ~{censat_h1l_bed} -wa > ~{pileup_h1l_intersect_bed_output}
	>>>

	output {
		File pileup_h1l_intersect_bed  = pileup_h1l_intersect_bed_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "biocontainers/bedtools:v2.28.0_cv2"
		preemptible: preempts
	}
}

task strict_detect {
	input {
		File pileup_bed
		File censat_h1l_bed
		String contig_name
		String sample_id

		Int memSizeGB = 16
		Int threadCount = 4
		Int diskSizeGB = 128
		Int preempts = 1
	}

	String strict_cdr_bed_output = "~{contig_name}_~{sample_id}.strictCDRs.bed"

	command <<<
		set -eux -o pipefail

		bash /opt/strictCDRDetection.sh \
			-i ~{pileup_bed} \
			-r ~{censat_h1l_bed} \
			-o "~{contig_name}_~{sample_id}"

		# for handling cases in which H1L doesn't have CDR(mostly broken contigs)
		touch ~{strict_cdr_bed_output}
	>>>

	output {
		File strict_cdr_bed = strict_cdr_bed_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "jmmenend/cdr_detect:0.4.1"
		preemptible: preempts
	}
}

task hmm_detect {
	input {
		File pileup_bed
		File strict_bed
		String contig_name
		String sample_id

		Int memSizeGB = 64
		Int threadCount = 4
		Int diskSizeGB = 128
		Int preempts = 1
	}

	String viterbi_cdr_bed_output = "~{contig_name}_~{sample_id}.hmmCDR.bed"

	String emissionBoundaries_output = "~{contig_name}_~{sample_id}.hmmCDR.emission_boundaries.csv"
	String emissionMatrix_output = "~{contig_name}_~{sample_id}.hmmCDR.emission_matrix.csv"
	String transitionMatrix_output = "~{contig_name}_~{sample_id}.hmmCDR.transition_matrix.csv"

	command <<<
		set -eux -o pipefail
		
		if [ -s ~{strict_bed} ]; then
			# run the HMM with a VERY low learning rate maybe 1e-7?
			python3 /opt/HMMCDRDetection.py \
				-l 0.0000001 \
				-p ~{pileup_bed} \
				-s ~{strict_bed} \
				-o ~{viterbi_cdr_bed_output}

			# sort the viterbi cdr bed file by position
			sort -k 1,1 -k2,2n -o ~{viterbi_cdr_bed_output} ~{viterbi_cdr_bed_output}
		else
			echo "No Strict CDRs on Contig Skipping HMM Detection!"
			touch ~{viterbi_cdr_bed_output}
			touch ~{emissionBoundaries_output}
			touch ~{emissionMatrix_output}
			touch ~{transitionMatrix_output}
		fi
	>>>

	output {
		File viterbi_cdr_bed = viterbi_cdr_bed_output
		
		File emissionBoundaries = emissionBoundaries_output
		File emissionMatrix = emissionMatrix_output
		File transitionMatrix = transitionMatrix_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "jmmenend/cdr_detect:0.4.1"
		preemptible: preempts
	}
}

task validate_cdrs {
	input {
		File pileup_bed
		File strict_bed
		File hmm_bed
		File emissionMatrix
		File transitionMatrix
		String contig_name
		String sample_id

		Int memSizeGB = 64
		Int threadCount = 4
		Int diskSizeGB = 128
		Int preempts = 1
	}
	
	String histogram_output = "~{contig_name}_~{sample_id}.hmmCDR.histogram.png"
	String emission_heatmap_output = "~{contig_name}_~{sample_id}.hmmCDR.emissionMatrixHeatmap.png"
	String transition_heatmap_output = "~{contig_name}_~{sample_id}.hmmCDR.transitionMatrixHeatmap.png"
	
	command <<<
		set -eux -o pipefail
		
		# need to fix each of these python scripts... :/
		if [ -s ~{strict_bed} ]; then
			if [ -s ~{hmm_bed} ]; then
				cat ~{strict_bed}
				cat ~{hmm_bed}
				# if strict_cdrs is NOT empty
				
				# generate graphs for faster CDR validation
				python3 /opt/cdr_histogram.py \
					-i ~{pileup_bed} \
					-r ~{hmm_bed} \
					-o "~{histogram_output}"

				python3 /opt/hmm_heatmaps.py \
					-e ~{emissionMatrix} \
					-t ~{transitionMatrix} \
					-o "~{contig_name}_~{sample_id}.hmmCDR"
			else
				# if hmm_cdrs IS empty
				echo "Contig has NO CDRs!"
			fi	
		else
			# if hmm_cdrs IS empty
			echo "Contig has NO CDRs!"
		fi

		touch ~{histogram_output}
		touch ~{emission_heatmap_output}
		touch ~{transition_heatmap_output}
	>>>

	output {
		File histogram = histogram_output
		File emission_heatmap = emission_heatmap_output
		File transition_heatmap = transition_heatmap_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "jmmenend/cdr_detect:0.4.1"
		preemptible: preempts
	}
}

task gather_outputs {
	input {
		Array[File] hmm_cdrs
		
		Array[File] pileup_beds
		Array[File] strict_cdrs_beds

		Array[File] pileup_histograms
		Array[File] emission_heatmaps
		Array[File] transition_heatmaps
		Array[File] emissionBoundaries
		Array[File] emissionMatrix
		Array[File] transitionMatrix

		String sample_id

		Int memSizeGB = 64
		Int threadCount = 1
		Int diskSizeGB = 128
		Int preempts = 1
	}

	String hmm_cdrs_bed_output = "~{sample_id}.hmmCDR.bed"
	
	String igv_zip_output = "~{sample_id}.igv.tar.gz"
	String summary_zip_output = "~{sample_id}.summary.tar.gz"

	command <<<
		set -eux -o pipefail

		# combine HMM CDR predictions 
		cat ~{sep=" " hmm_cdrs} | sort -k 1,1 -k2,2n > ~{hmm_cdrs_bed_output}
		
		# combine H1L pileups
		cat ~{sep=" " pileup_beds} | sort -k 1,1 -k2,2n > "~{sample_id}.pileup.H1L.bedgraph"
		
		# combine strict CDRs and Transitions
		cat ~{sep=" " strict_cdrs_beds} | sort -k 1,1 -k2,2n > "~{sample_id}.strictCDR.bed"

		# zip igv tracks together
		mkdir -p igv_tracks
		cp "~{sample_id}.pileup.H1L.bedgraph" igv_tracks
		cp "~{sample_id}.strictCDR.bed" igv_tracks
		cd igv_tracks && \
			tar -czf ../~{igv_zip_output} . && \
			cd ..
		
		# Zip summary files together
		mkdir -p summary_zip
		for file in ~{sep=" " pileup_histograms}; do
			file_basename=$(basename $file)
			find . -type f -name $file_basename -size +0 -exec cp {} summary_zip \;
		done

		for file in ~{sep=" " emission_heatmaps}; do
			file_basename=$(basename $file)
			find . -type f -name $file_basename -size +0 -exec cp {} summary_zip \;
		done

		for file in ~{sep=" " transition_heatmaps}; do
			file_basename=$(basename $file)
			find . -type f -name $file_basename -size +0 -exec cp {} summary_zip \;
		done

		for file in ~{sep=" " emissionBoundaries}; do
			file_basename=$(basename $file)
			find . -type f -name $file_basename -size +0 -exec cp {} summary_zip \;
		done

		for file in ~{sep=" " emissionMatrix}; do
			file_basename=$(basename $file)
			find . -type f -name $file_basename -size +0 -exec cp {} summary_zip \;
		done
		
		for file in ~{sep=" " transitionMatrix}; do
			file_basename=$(basename $file)
			find . -type f -name $file_basename -size +0 -exec cp {} summary_zip \;
		done

		cd summary_zip && \
			tar -czf ../~{summary_zip_output} . && \
			cd ..
	>>>

	output {
		File hmm_cdrs_bed = hmm_cdrs_bed_output
		
		File igv_zip = igv_zip_output
		File summary_zip = summary_zip_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "jmmenend/jq:v1.0"
		preemptible: preempts
	}
}
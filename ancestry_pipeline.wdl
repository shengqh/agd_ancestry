version 1.0 

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/genotype/Utils.wdl" as http_GenotypeUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/agd/AgdUtils.wdl" as http_AgdUtils


# WORKFLOW

workflow agd_ancestry_workflow{
    input{
        # workflow choices 

        Boolean external_spike_in = true 
        Boolean scope_supervised = true
        Boolean run_pca = true
        Boolean run_scope= true

        # required inputs original data as array of chromosomes 

        Array[File] source_pgen_files
        Array[File] source_pvar_files
        Array[File] source_psam_files
        Array[String] chromosomes
        String target_prefix

        File id_map_file

        # optional outputs for exporting 
        String? project_id
        String target_gcp_folder

        # optional inputs for PCA - required if running PCA 

        File? pca_variants_extract_file
        File? pca_loadings_file
        File? pca_af_file


        # optional inputs for spike in data - required if merging spike in data 

        Array[File]? source_bed_files
        Array[File]? source_bim_files
        Array[File]? source_fam_files
        File? spike_in_pgen_file
        File? spike_in_pvar_file
        File? spike_in_psam_file
        File? spike_in_relatives_exclude 
    
        # optional inputs for supervised scope    - required if running supervised scope 

        File? supervised_scope_reference_pgen_file
        File? supervised_scope_reference_pvar_file
        File? supervised_scope_reference_psam_file
        File? supervised_scope_reference_superpop_file   
        File? supervised_scope_reference_relatives_exclude
        File? supervised_scope_reference_freq

        #optional inputs for unsupervised scope - required if running unsupervised scope
        String? scope_plink2_maf_filter = "--maf 0.01"

        String? scope_plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
        File? scope_long_range_ld_file
        Int? K = 4
        Int? seed = 1234
    }

    # If the user chose to use the supervised scope and there is no precalculated reference allele frequency provided, allele frequency can be calculated from provided plink files instead
    if(scope_supervised && !defined(supervised_scope_reference_freq)){
        File required_supervised_scope_reference_pgen_file = select_first([supervised_scope_reference_pgen_file])
        File required_supervised_scope_reference_pvar_file = select_first([supervised_scope_reference_pvar_file])
        File required_supervised_scope_reference_psam_file = select_first([supervised_scope_reference_psam_file])
        File required_supervised_scope_reference_superpop_file = select_first([supervised_scope_reference_superpop_file])

        call CalculateFreq{
        input: 
            pgen_file = required_supervised_scope_reference_pgen_file,
            pvar_file = required_supervised_scope_reference_pvar_file,
            psam_file = required_supervised_scope_reference_psam_file,
            superpop_file = required_supervised_scope_reference_superpop_file,
            relatives_exclude = supervised_scope_reference_relatives_exclude
        }
        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFileFreq{
                input:
                    source_file = CalculateFreq.freq_file,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }
    }

    # If the user chose to use an external spike in, then first the spike-in data must be merged with the original data, and then the pipeline can proceed 
    if(external_spike_in){
        if(!defined(source_bed_files)){
            scatter (idx in range(length(chromosomes))) {
                String chromosome_for_spike_in_conversion = chromosomes[idx]
                File pgen_file_for_conversion = source_pgen_files[idx]
                File pvar_file_for_conversion = source_pvar_files[idx]
                File psam_file_for_conversion = source_psam_files[idx]

                call ConvertPgenToBed as ConvertPgenToBedForSpikeIn{
                    input:
                        pgen = pgen_file_for_conversion,
                        pvar = pvar_file_for_conversion,
                        psam = psam_file_for_conversion, 
                        out_prefix = chromosome_for_spike_in
                }
            }
         } 
        Array[File] source_bed_files_required=select_first([source_bed_files, ConvertPgenToBedForSpikeIn.convert_Pgen_out_bed])
        Array[File] source_bim_files_required=select_first([source_bim_files, ConvertPgenToBedForSpikeIn.convert_Pgen_out_bim])
        Array[File] source_fam_files_required=select_first([source_fam_files, ConvertPgenToBedForSpikeIn.convert_Pgen_out_bim])

         scatter (idx in range(length(chromosomes))) {
            String chromosome_for_spike_in = chromosomes[idx]
            File bed_file_for_spike_in = source_bed_files_required[idx]
            File bim_file_for_spike_in = source_bim_files_required[idx]
            File fam_file_for_spike_in = source_fam_files_required[idx]

            call SubsetChromosomeTGP{
                input: 
                    pgen_file = spike_in_pgen_file,
                    pvar_file = spike_in_pvar_file,
                    psam_file = spike_in_psam_file,
                    chromosome = chromosome_for_spike_in,
                    relatives_exclude = spike_in_relatives_exclude
            }

            call Merge1000genomesAGD {
                input:
                    agd_bed_file = bed_file_for_spike_in,
                    agd_bim_file = bim_file_for_spike_in,
                    agd_fam_file = fam_file_for_spike_in,
                    TGP_bed_file = SubsetChromosomeTGP.subset_reference_out_bed_file,
                    TGP_bim_file = SubsetChromosomeTGP.subset_reference_out_bim_file,
                    TGP_fam_file = SubsetChromosomeTGP.subset_reference_out_fam_file,
                    chromosome = chromosome_for_spike_in
            }
        }
    }

    #now we can proceed with the pipeline, using select_first to get the first non empty array value to select either the merged or original input files 

    Array[File] my_pgen_files=select_first([Merge1000genomesAGD.merge_spike_out_pgen_file, source_pgen_files])
    Array[File] my_pvar_files=select_first([Merge1000genomesAGD.merge_spike_out_pvar_file, source_pvar_files])
    Array[File] my_psam_files=select_first([Merge1000genomesAGD.merge_spike_out_psam_file, source_psam_files])

    if(run_pca){
        scatter (idx in range(length(chromosomes))) {
            String chromosome_for_pca = chromosomes[idx]
            File pgen_file_for_pca = my_pgen_files[idx]
            File pvar_file_for_pca = my_pvar_files[idx]
            File psam_file_for_pca = my_psam_files[idx]
            String replaced_sample_name_for_pca = "~{chromosome_for_pca}.psam"

            #I think I need this to get the IDs correctly as GRIDS

            call http_AgdUtils.ReplaceICAIdWithGrid as ReplaceICAIdWithGridForPCA {
                input:
                    input_psam = psam_file_for_pca,
                    id_map_file = id_map_file,
                    output_psam = replaced_sample_name_for_pca
            }

            call ExtractVariants as ExtractVariants{
                input:
                    pgen_file = pgen_file_for_pca,
                    pvar_file = pvar_file_for_pca,
                    psam_file = ReplaceICAIdWithGridForPCA.output_psam,
                    chromosome = chromosome_for_pca,
                    variants_extract_file = pca_variants_extract_file
            }
        }

        call http_GenotypeUtils.MergePgenFiles as MergePgenFilesForPCA{
            input:
                pgen_files = ExtractVariants.extract_variants_output_pgen_file,
                pvar_files = ExtractVariants.extract_variants_output_pvar_file,
                psam_files = ExtractVariants.extract_variants_output_psam_file,
                output_prefix = target_prefix
        }

        call ProjectPCA{
            input: 
                pgen_file = MergePgenFilesForPCA.output_pgen, 
                pvar_file = MergePgenFilesForPCA.output_pvar,
                psam_file = MergePgenFilesForPCA.output_psam,  
                PCA_loadings = pca_loadings_file,
                PCA_AF = pca_af_file,
                OUTNAME = target_prefix
        }

        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFile_PCAone {
                input:
                    source_file = ProjectPCA.output_pca_file,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }

        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFile_PCAtwo {
                input:
                    source_file = ProjectPCA.output_pca_variants,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }
    }

    if(run_scope){
        scatter (idx in range(length(chromosomes))) {
            String chromosome_for_scope = chromosomes[idx]
            File pgen_file_for_scope = my_pgen_files[idx]
            File pvar_file_for_scope = my_pvar_files[idx]
            File psam_file_for_scope = my_psam_files[idx]
            String replaced_sample_name_for_scope = "~{chromosome_for_scope}.psam"

            #I think I need this to get the IDs correctly as GRIDS

            call http_AgdUtils.ReplaceICAIdWithGrid as ReplaceICAIdWithGridForScope {
                input:
                    input_psam = psam_file_for_scope,
                    id_map_file = id_map_file,
                    output_psam = replaced_sample_name_for_scope
            }
            call PreparePlink as PreparePlink{
                input:
                    pgen_file = pgen_file_for_scope,
                    pvar_file = pvar_file_for_scope,
                    psam_file = ReplaceICAIdWithGridForScope.output_psam,
                    long_range_ld_file = scope_long_range_ld_file,
                    plink2_maf_filter = scope_plink2_maf_filter,
                    plink2_LD_filter_option = scope_plink2_LD_filter_option,
                    chromosome = chromosome_for_scope 
            }
        }
        call http_GenotypeUtils.MergePgenFiles as MergePgenFilesForScope{
            input:
                pgen_files = PreparePlink.prepare_plink_unsupervised_output_pgen_file,
                pvar_files = PreparePlink.prepare_plink_unsupervised_output_pvar_file,
                psam_files = PreparePlink.prepare_plink_unsupervised_output_psam_file,
                output_prefix = target_prefix
        }

        call ConvertPgenToBed as ConvertPgenToBedForScope{
            input: 
                pgen = MergePgenFilesForScope.output_pgen, 
                pvar = MergePgenFilesForScope.output_pvar,
                psam = MergePgenFilesForScope.output_psam, 
        }

        call RunScopeUnsupervised{    
            input:
                bed_file = ConvertPgenToBedForScope.convert_Pgen_out_bed,
                bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
                fam_file = ConvertPgenToBedForScope.convert_Pgen_out_fam,
                K = K,
                output_string = target_prefix,
                seed = seed
        }

        if(scope_supervised){
            call QCAllelesBim{
                input:
                    bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
                    freq_file = select_first([CalculateFreq.freq_file, supervised_scope_reference_freq])
            }

            call PreparePlinkSupervised{
                input:
                    bed_file = ConvertPgenToBedForScope.convert_Pgen_out_bed,
                    bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
                    fam_file = ConvertPgenToBedForScope.convert_Pgen_out_fam,
                    variant_list = QCAllelesBim.out_variants
            }

            call RunScopeSupervised{
                input:
                    bed_file = PreparePlinkSupervised.prepare_plink_supervised_out_bed,
                    bim_file = PreparePlinkSupervised.prepare_plink_supervised_out_bim,
                    fam_file = PreparePlinkSupervised.prepare_plink_supervised_out_fam,
                    K = K,
                    output_string = target_prefix,
                    seed = seed,
                    topmed_freq = QCAllelesBim.out_frq
            }
        }
        
        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyThreeFiles as CopyFiles_one{
                input:
                    source_file1 = RunScopeUnsupervised.outP,
                    source_file2 = RunScopeUnsupervised.outQ,
                    source_file3 = RunScopeUnsupervised.outV,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
            if(scope_supervised){
                call http_GcpUtils.MoveOrCopyThreeFiles as CopyFiles_two {
                    input:
                        source_file1 = select_first([RunScopeSupervised.outP]),
                        source_file2 = select_first([RunScopeSupervised.outQ]),
                        source_file3 = select_first([RunScopeSupervised.outV]),
                        is_move_file = false,
                        project_id = project_id,
                        target_gcp_folder = select_first([target_gcp_folder])
                }
            }
        }
    }

 output {
        Array[File] ancestry_outputs = select_all([CopyFileFreq.output_file, CopyFile_PCAone.output_file, CopyFile_PCAtwo.output_file, CopyFiles_one.output_file1, CopyFiles_one.output_file2, CopyFiles_one.output_file3, CopyFiles_two.output_file1, CopyFiles_two.output_file2, CopyFiles_two.output_file3])
    }
}


# TASKS

## for frequency calculation 
task CalculateFreq{
    input{ 
        File pgen_file
        File pvar_file
        File psam_file
        File superpop_file
        File? relatives_exclude

        String? plink2_maf_filter = "--maf 0.001"

        Int memory_gb = 20
        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 4) + 20
    String out_name =  basename(pgen_file, ".pgen") + "_within_superpop_freqs"
    String out_file = out_name + ".frq.strat"

    command <<<
        # take the TGP data, remove duplicates, restrict to biallelic SNPs (necessary since cannot calculate frq within in plink2 in the same way as in plink1), and calculate within super populations 
        # include a MAF filter (default 0.001) to reduce size of output file
        plink2 --pgen ~{pgen_file} \
            --pvar ~{pvar_file} \
            --psam ~{psam_file} \
            --allow-extra-chr \
            --chr 1-22, X, Y \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 10 truncate \
            --max-alleles 2 \
            ~{plink2_maf_filter} \
            --rm-dup 'exclude-all' \
            --remove ~{relatives_exclude} \
            --make-bed \
            --out tgp_nodup  

        plink --bfile tgp_nodup \
            --allow-extra-chr \
            --freq --within ~{superpop_file} \
            --out ~{out_name}
    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output{
        File freq_file = out_file
    }
}

## for merging with spike-in data task 
task SubsetChromosomeTGP {
    input {
        File? pgen_file
        File? pvar_file
        File? psam_file
        String chromosome
        File? relatives_exclude
        Int? memory_gb = 20
        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }
    
    String out_string = "TGP_" + chromosome
    String in_chromosome = sub(chromosome, "chr", "")
    Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB") * 2) * 2 + 20
    
    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }
    
    command {
         plink2 \
             --pgen ~{pgen_file} --pvar ~{pvar_file} --psam ~{psam_file} \
             --allow-extra-chr \
             --chr ~{in_chromosome} \
             --remove ~{relatives_exclude} \
             --make-bed \
             --out ~{out_string}
    }
    
    output {
        File subset_reference_out_bed_file = out_string + ".bed"
        File subset_reference_out_bim_file = out_string + ".bim"
        File subset_reference_out_fam_file = out_string + ".fam"
    }
}


task ConvertPgenToBed{
    input {
        File pgen 
        File pvar 
        File psam 

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

        String? out_prefix

        Int? memory_gb = 20
    }

    Int disk_size = ceil(size([pgen, pvar, psam], "GB")  * 2)*2 + 20


    String out_string = if defined(out_prefix) then out_prefix else basename(pgen, ".pgen")

    command {
        plink2 \
            --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --make-bed \
            --out ~{out_string}
    }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File convert_Pgen_out_bed = "${out_string}.bed"
        File convert_Pgen_out_bim = "${out_string}.bim"
        File convert_Pgen_out_fam = "${out_string}.fam"
    }
}

task Merge1000genomesAGD{
    input{
        File agd_bed_file
        File agd_bim_file
        File agd_fam_file

        File TGP_bed_file
        File TGP_bim_file
        File TGP_fam_file

        String? chromosome

        Int? memory_gb = 20


        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([agd_bed_file, TGP_bed_file], "GB")  * 2)*2 + 20


    String out_string = "AGD_TGP_" + chromosome
    String agd_prefix = basename(agd_bed_file, ".bed")
    String TGP_prefix = basename(TGP_bed_file, ".bed")
    String agd_prefix_rename= agd_prefix + "_renamed"
    String agd_prefix_2 = agd_prefix + "_2"
    String TGP_prefix_2 = TGP_prefix + "_2"

    String relocated_bed = agd_prefix + ".bed"
    String relocated_bim = agd_prefix + ".bim"
    String relocated_fam = agd_prefix + ".fam"

    String relocated_tgp_bed = TGP_prefix + ".bed"
    String relocated_tgp_bim = TGP_prefix + ".bim"
    String relocated_tgp_fam = TGP_prefix + ".fam"

    String agd_snp_list = agd_prefix + "_renamed.snplist"
    String TGP_snp_list = TGP_prefix + ".snplist"


    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    command{

        ln ~{agd_bed_file} ./~{relocated_bed}
        ln ~{agd_bim_file} ./~{relocated_bim}
        ln ~{agd_fam_file} ./~{relocated_fam}

        ln ~{TGP_bed_file} ./~{relocated_tgp_bed}
        ln ~{TGP_bim_file} ./~{relocated_tgp_bim}
        ln ~{TGP_fam_file} ./~{relocated_tgp_fam}

        plink2 \
            --bfile ~{agd_prefix} \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 10 truncate \
            --make-bed \
            --out ~{agd_prefix_rename}

        plink2 \
            --bfile ~{agd_prefix_rename} \
            --write-snplist \
            --out ~{agd_prefix_rename} 
        
        plink2 \
            --bfile ~{TGP_prefix} \
            --write-snplist \
            --out ~{TGP_prefix}

        plink \
            --bfile ~{agd_prefix_rename} \
            --bmerge ~{TGP_prefix} \
            --make-bed \
            --out merged_beds_files

        plink2 \
            --bfile ~{agd_prefix_rename} \
            --exclude merged_beds_files-merge.missnp \
            --extract-intersect ~{agd_snp_list} ~{TGP_snp_list} \
            --make-bed \
            --out ~{agd_prefix_2}

        plink2 \
            --bfile ~{TGP_prefix} \
            --exclude merged_beds_files-merge.missnp \
            --extract-intersect ~{agd_snp_list} ~{TGP_snp_list} \
            --make-bed \
            --out ~{TGP_prefix_2}

        plink \
            --bfile ~{agd_prefix_2} \
            --bmerge ~{TGP_prefix_2} \
            --make-bed \
            --out merged_beds_files2
        
        plink2 \
            --bfile merged_beds_files2 \
            --make-pgen \
            --out ~{out_string}
    }

    output{
        File merge_spike_out_pgen_file = out_string + ".pgen"
        File merge_spike_out_pvar_file = out_string + ".pvar"
        File merge_spike_out_psam_file = out_string + ".psam"
    }
}

## for PCA 
task ExtractVariants{
  input {
    File? pgen_file
    File? pvar_file
    File? psam_file 

    String? chromosome

    File? variants_extract_file

    Int? memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"
  String intermediate_pgen = chromosome + "_varids.pgen"
  String intermediate_pvar = chromosome + "_varids.pvar"
  String intermediate_psam = chromosome + "_varids.psam"

  command {
    plink2 \
      --pgen ~{pgen_file} \
      --pvar ~{pvar_file} \
      --psam ~{psam_file} \
      --snps-only \
      --set-all-var-ids chr@:#:\$r:\$a \
      --new-id-max-allele-len 1000 \
      --make-pgen \
      --out ~{chromosome}_varids
    
    plink2 \
      --pgen ~{intermediate_pgen} \
      --pvar ~{intermediate_pvar} \
      --psam ~{intermediate_psam} \
      --extract ~{variants_extract_file} \
      --make-pgen \
      --out ~{chromosome}
  }

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }

  output {
    File extract_variants_output_pgen_file = new_pgen
    File extract_variants_output_pvar_file = new_pvar
    File extract_variants_output_psam_file = new_psam
  }

}

task ProjectPCA{
  input{
    File? pgen_file
    File? pvar_file
    File? psam_file
    File? PCA_loadings
    File? PCA_AF
    String? OUTNAME

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

  }

  Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB")  * 2) + 20

  String pca_file = OUTNAME + ".genotype.pca.sscore"
  String pca_variants = OUTNAME + "genotype.sscore.vars"

  command {
    plink2 --pgen ~{pgen_file} --pvar ~{pvar_file} --psam ~{psam_file} --score ~{PCA_loadings} \
    variance-standardize \
    cols=-scoreavgs,+scoresums \
    list-variants \
    header-read \
    --score-col-nums 3-22 \
    --read-freq ~{PCA_AF} \
    --out ~{OUTNAME}

    cp ~{OUTNAME}.sscore ~{pca_file}
    cp ~{OUTNAME}.sscore.vars ~{pca_variants}
    }

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  
  output {
    File output_pca_file = "~{pca_file}"
    File output_pca_variants="~{pca_variants}"
  }
}

## for scope
task PreparePlink{
  input {
    File? pgen_file
    File? pvar_file
    File? psam_file 

    String? chromosome

    String? plink2_maf_filter = "--maf 0.01"
    String? plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
    File? long_range_ld_file


    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"
  String out_prefix = chromosome 


  command {
    plink2 \
      --pgen ~{pgen_file} \
      --pvar ~{pvar_file} \
      --psam ~{psam_file} \
      ~{plink2_maf_filter} \
      --snps-only \
      --const-fid \
      --set-all-var-ids chr@:#:\$r:\$a \
      --new-id-max-allele-len 1000 \
      --make-pgen \
      --out maf_filtered
      
    plink2 \
        --pgen maf_filtered.pgen \
        --pvar maf_filtered.pvar \
        --psam maf_filtered.psam \
        --exclude range ~{long_range_ld_file} \
        --make-pgen \
        --out maf_filtered_longrange
    
    plink2 \
      --pgen maf_filtered_longrange.pgen \
      --pvar maf_filtered_longrange.pvar \
      --psam maf_filtered_longrange.psam \
      ~{plink2_LD_filter_option}

    plink2 \
        --pgen maf_filtered_longrange.pgen \
        --pvar maf_filtered_longrange.pvar \
        --psam maf_filtered_longrange.psam \
        --extract plink2.prune.in \
        --make-pgen \
        --out ~{out_prefix}
  }

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }

  output {
    File prepare_plink_unsupervised_output_pgen_file = new_pgen
    File prepare_plink_unsupervised_output_pvar_file = new_pvar
    File prepare_plink_unsupervised_output_psam_file = new_psam
  }
}

task QCAllelesBim{
    input {
        File? bim_file
        File? freq_file

        String docker = "blosteinf/r_utils_terra:0.1"
        Int memory_gb = 20
    }

    Int disk_size = ceil(size([bim_file, freq_file], "GB")  * 2) + 20

    command {
        ls /home/r-environment/
        Rscript /home/r-environment/allele_qc.R --in_freq ~{freq_file} --in_bim ~{bim_file} 
    }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File out_frq = "corrected_freq.frq"
        File out_variants = "variants_to_extract.txt"
    }
}

task PreparePlinkSupervised{
    input { 
        File? bed_file
        File? bim_file
        File? fam_file 

        File? variant_list 
        String? out_string = "variant_filtered"

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

  Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20
  Int memory_gb = 20

  String new_bed = out_string + ".bed"
  String new_bim = out_string + ".bim"
  String new_fam= out_string + ".fam"

  command { 
    plink2 \
        --bed ~{bed_file} \
        --bim ~{bim_file} \
        --fam ~{fam_file} \
        --extract ~{variant_list} \
        --make-bed \
        --out ~{out_string}
  }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File prepare_plink_supervised_out_bed = new_bed
        File prepare_plink_supervised_out_bim = new_bim
        File prepare_plink_supervised_out_fam = new_fam
        String out_prefix = out_string
    }

}

task RunScopeUnsupervised{
    input{

        File bed_file
        File bim_file
        File fam_file

        Int? K
        String output_string
        Int? seed

        Int memory_gb = 60
        String docker = "blosteinf/scope:0.1"
    }

    String plink_binary_prefix =  basename(bed_file, ".bed")
    String relocated_bed = plink_binary_prefix + ".bed"
    String relocated_bim = plink_binary_prefix + ".bim"
    String relocated_fam = plink_binary_prefix + ".fam"

    String unsup_output = output_string + "_unsupervised_" 
    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<
        ln ~{bed_file} ./~{relocated_bed}
        ln ~{bim_file} ./~{relocated_bim}
        ln ~{fam_file} ./~{relocated_fam}
        scope -g ~{plink_binary_prefix} -k ~{K} -seed ~{seed} -o ~{unsup_output}
        awk '{ for (i=1; i<=NF; i++) { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++) { str=str" "a[i,j]; } print str } }' ~{unsup_output}Qhat.txt > transposed_Qhat.txt
        cut -f2 ./~{relocated_fam} | paste - transposed_Qhat.txt > ~{unsup_output}Qhat.txt
    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
  }

    output {
        File outP= "${unsup_output}Phat.txt"
        File outQ= "${unsup_output}Qhat.txt"
        File outV= "${unsup_output}V.txt"
    }
}

task RunScopeSupervised{
    input{
       
        File bed_file
        File bim_file
        File fam_file

        Int? K
        String output_string
        Int? seed

        File? topmed_freq

        Int memory_gb = 60

        String docker = "blosteinf/scope:0.1"
    }

    String plink_binary_prefix = basename(bed_file, ".bed")
    String relocated_bed= plink_binary_prefix + ".bed"
    String relocated_bim= plink_binary_prefix + ".bim"
    String relocated_fam= plink_binary_prefix + ".fam"

    String sup_output = output_string + "_supervised_"

    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<
        ln ~{bed_file} ./~{relocated_bed}
        ln ~{bim_file} ./~{relocated_bim}
        ln ~{fam_file} ./~{relocated_fam}
        scope -g ~{plink_binary_prefix} -freq ~{topmed_freq} -k ~{K} -seed ~{seed} -o ~{sup_output}
        ls
        awk '{ for (i=1; i<=NF; i++) { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++) { str=str" "a[i,j]; } print str } }' ~{sup_output}Qhat.txt > transposed_Qhat.txt
        cut -f2 ./~{relocated_fam} | paste - transposed_Qhat.txt > ~{sup_output}Qhat.txt
    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
  }

    output {
        File outP= "${sup_output}Phat.txt"
        File outQ= "${sup_output}Qhat.txt"
        File outV= "${sup_output}V.txt"
    }
}


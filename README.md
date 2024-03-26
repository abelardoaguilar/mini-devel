# mini-devel

```mermaid
flowchart TD
    subgraph ide1 [notation]
    B(folder)
    C[[script]]
    D[(database)]
    end


```
Development repo for the miniphages project

```mermaid
flowchart TD

    classDef src stroke:#f00
    classDef bin stroke:#0f0
    classDef results stroke:#00f
    classDef data stroke:#f96

    A[mini-devel] --- S(src):::src
    A --- B(bin):::bin
    A --- R(results):::results
    A --- D(data):::data
    



    B --- B1[[1_split_global_DTR_fna_by_biomes.ipynb]]:::bin -.- B2[[2_execute_VIBRANT_in_site.sh]]:::bin -.- B3[[3_filter_VIBRANT_output.ipynb]]:::bin

    B1 -.-> R1
    D1 ---> R1

    D --- 1_environmental_seqs_from_JGI_IMGvir[1_environmental_seqs_from_JGI_IMGvir]:::data

    subgraph ide1 [files shared by JGI]
    1_environmental_seqs_from_JGI_IMGvir --- D1[(DTRs_20kb.fna\nDTRs_20kb.csv)]:::data
    1_environmental_seqs_from_JGI_IMGvir --- D3[(biome1.tsv\nbiome2.tsv\nbiome3.tsv)]:::data
    end

    R --- R1(1_biome_fasta_files):::results

    subgraph ide2 [output for: 1_split_global_DTR_fna_by_biomes.ipynb]
    R1 --- biome0[(DTRs_20kb)]:::results
    R1 --- biome1[(biome1)]:::results
    R1 --- biome2[(biome2)]:::results
    R1 --- biome3[(biome3)]:::results
    end

    subgraph ide3 [VIBRANT output folders]
    biome0 ---> vibrant0[(VIBRANT OUTPUTS\n\nDTRs_20kb\nbiome1\nbiome2\nbiome3)]:::results
    biome1 ---> vibrant0
    biome2 ---> vibrant0
    biome3 ---> vibrant0
    B2 -.-> vibrant0 
    end

    S --- auxiliary_functions(auxiliary_functions):::src
    auxiliary_functions --- S1[[__init__.py]]:::src --- S2[[VIBRANT_filtering_functions.py]]:::src 

    vibrant0 ---> 2_filtered_VIBRANT_output(2_filtered_VIBRANT_output\n*one per VIBRANT output):::results


    subgraph ide4 [filtered VIBRANT outputs]
    B3 -.-> 2_filtered_VIBRANT_output
    2_filtered_VIBRANT_output --- filters_by_mcp[(FILTERS BY MCP\n\nProteins_identified_as_MCPs.faa\nProteins_identified_as_not_MCPs.faa\n\n_contigs_with/without_mcp.fna .faa\n\n_contigs_w_mcp_.tsv)]:::results
    2_filtered_VIBRANT_output --- filters_by_quality[(FILTERS BY QUALITY\n\n_above_quality_genomes_.fna .faa\n_bellow_quality_genomes_.fna .faa\n\n_quality_genomes_.csv)]:::results
    2_filtered_VIBRANT_output --- filters_by_length[(FILTERS BY CONTIG LENGTH\n\n_contigs_within_range.fna .faa\n\n_metadata_filtered_by_size.csv)]:::results

    filters_by_mcp --- multifilter[(\nCOMBINATIONS OF FILTERS\n+\nphages_meeting_all_criteria.fna .faa)]:::results
    filters_by_quality --- multifilter
    filters_by_length --- multifilter
    end
    
    vibrant0 --- |output 23|2_global_VIBRANT_output23[3_All_phages_combined]:::results

    subgraph ide5 [All encoded proteins found by VIBRANT]
    2_global_VIBRANT_output23 --- output_23[(All_phages_combined.fna\n.faa\n.ffn\n.gbk\n.txt)]:::results
    end

    



    D --- efams[2_Zayed_efam_2020]:::data
    
    subgraph ide6 [efam Viral Protein Families]
    efams --- cluster41798[(cluster41798\n\n.csv\n.faa\n_aligned.faa\n.hmm)]:::data
    end
    
    subgraph ide7 [viral-protein-function-plm]
    output_23 --- |faa_prediction.sh|viral_pllm[(*_embeddings_dict.pkl
*_functional_probabilities.csv
*_function_predictions.csv)]:::results
    end

    subgraph ide8 [Search of new MCP detected by Zachary 2023]
    output_23 --- |hmmsearch|hmm[(All_phages_combined.hmmsearch*)]:::results
    cluster41798 -.-> |hmmbuild|hmm
    end
    
    

```
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
    A[mini-devel] --- S(src)
    A --- B(bin)
    A --- R(results)
    A --- D(data)


    B --- B1[[1_split_global_DTR_fna_by_biomes.ipynb]] -.- B2[[2_execute_VIBRANT_in_site.sh]] -.- B3[[3_filter_VIBRANT_output.ipynb]]

    B1 -.-> |splits DTRs_20kb.fna\ninto biomes| R1

    D1 ---> |input| R1

    subgraph ide1 [files shared by JGI]
    D --- D1[(DTRs_20kb.fna)]
    D --- D2[(DTRs_20kb.csv)]
    D --- D3[(biome1.tsv)]
    D --- D4[(biome2.tsv)]
    D --- D5[(biome3.tsv)]
    end

    R --- R1(1_biome_fasta_files)

    subgraph ide2 [output for: 1_split_global_DTR_fna_by_biomes.ipynb]
    R1 --- biome0[(biome0)]
    R1 --- biome1[(biome1)]
    R1 --- biome2[(biome2)]
    R1 --- biome3[(biome3)]
    end

    subgraph ide3 [VIBRANT output folders]
    biome0 ---> vibrant0[(VIBRANT_DTRs_20kb)]
    B2 -.-> vibrant0 
    biome1 --- vibrant_all[(VIBRANT_biome ...)]
    biome2 --- vibrant_all
    biome3 --- vibrant_all
    end

    S --- auxiliary_functions(auxiliary_functions)
    auxiliary_functions --- S1[[__init__.py]] --- S2[[VIBRANT_filtering_functions.py]] 

    vibrant_all ---> 2_filtered_VIBRANT_output(2_filtered_VIBRANT_output\n*one per VIBRANT output)
    vibrant_all ---> 2_filtered_VIBRANT_output
    vibrant_all ---> 2_filtered_VIBRANT_output
    vibrant0 ---> 2_filtered_VIBRANT_output


    subgraph ide4 [filtered VIBRANT outputs]
    B3 -.-> 2_filtered_VIBRANT_output
    2_filtered_VIBRANT_output --- _MCPs[(_MCPs.faa)]
    2_filtered_VIBRANT_output --- _phages_with_MCP[(_phages_with_MCP.faa)]
    2_filtered_VIBRANT_output --- complete_phages[(_complete_phages.faa)]
    2_filtered_VIBRANT_output --- complete_phages_with_MCP[(_complete_phages_with_MCP.faa)]
    end


```
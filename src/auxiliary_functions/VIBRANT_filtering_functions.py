from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import os


def read_VIBRANT_custom_fasta_23(file_path):
    """
    Reads a custom FASTA file in the VIBRANT format and returns a list of records.

    Args:
        file_path (str): The path to the input file.

    Returns:
        list: A list of dictionaries, where each dictionary represents a record in the FASTA file.
              Each record dictionary contains the following keys:
              - 'id': The identifier of the record.
              - 'range': The range of the record.
              - 'strand': The strand of the record.
              - 'annotation_code': The annotation code of the record.
              - 'description': The description of the record.
              - 'sequence': The sequence of the record.

    """
    records = []
    with open(file_path, 'r') as file:
        record = {}
        for line in file:
            if line.startswith('>'):  # Header line indicates a new record
                if record:  # If there's an existing record, save it before starting a new one
                    records.append(record)
                parts = line[1:].strip().split('\t')  # Split header line into components
                record = {
                    'id': parts[0],
                    'range': parts[1],
                    'strand': parts[2],
                    'annotation_code': parts[3],
                    'description': parts[4],
                    'sequence': ''
                }
            else:  # Sequence line, might span multiple lines
                record['sequence'] += line.strip()
        if record:  # Don't forget to save the last record
            records.append(record)
    return records


def search_identifiers_in_VIBRANT_phages_proteins(path,VIBRANT_phages_proteins):
    
    def read_fasta_file(path):
        identifiers = []
        for record in SeqIO.parse(path, "fasta"):
            identifiers.append(record.id)
        return identifiers
    
    
    faa_output_path = path.replace('.fna', '.faa')
    genome_identifiers = read_fasta_file(path)
    VIBRANT_phages_proteins_ids = [record['id'] for record in VIBRANT_phages_proteins]
    modified_protein_ids = ['_'.join(protein_id.split('_')[:2]) for protein_id in VIBRANT_phages_proteins_ids]
    df_ids = pd.DataFrame({'id':modified_protein_ids, 'VIBRANT_phages_proteins_ids':VIBRANT_phages_proteins_ids})
    matching_proteins_ids = df_ids[df_ids['id'].isin(genome_identifiers)]['VIBRANT_phages_proteins_ids']
    matching_proteins_ids = matching_proteins_ids.tolist()    

    matching_protein_sequences = [SeqRecord(Seq(record['sequence']), id=record['id'], description='') for record in VIBRANT_phages_proteins if record['id'] in matching_proteins_ids]
    SeqIO.write(matching_protein_sequences, faa_output_path, "fasta")
    
    print('\nFrom the file:', path)
    print('\tNumber of genomes:', len(genome_identifiers), 'Number of matching proteins:', len(matching_proteins_ids))
    print('\tMatching proteins written to:', faa_output_path)
    
    return matching_protein_sequences



def filter_VIBRANT_annotations_by_mcp(df, mcp_terms=None, false_terms=None):
    """
    Filters a DataFrame based on annotations related to major coat proteins (MCP) across multiple columns.
    Only columns with 'name' in their titles are considered for filtering to find or exclude specific terms.

    Args:
        df (pandas.DataFrame): The DataFrame to filter.
        mcp_terms (list, optional): Terms related to MCP. Defaults to a predefined list if None.
        false_terms (list, optional): Terms to identify and exclude false positives. Defaults to a predefined list if None.

    Returns:
        pandas.DataFrame: A filtered DataFrame containing only rows likely related to MCPs, across specified columns.
    """
    if mcp_terms is None:
        mcp_terms = ["mcp", "major", "coat", "capsid"
        ]
    if false_terms is None:
        false_terms = [
            "minor", "fiber", "tropism", "non-structural", "envelope protein",
            "replicase", "polymerase", "regulatory protein", "accessory protein", "tail", "assembly",
            "protease", "encapsidation"
        ]
    
    # Identify columns with 'name' in their title to filter based on annotations within these columns
    name_columns = [col for col in df.columns if 'name' in col.lower()]
    
    # Filter rows by checking for MCP-related terms and excluding false positive terms across the identified columns
    
    def is_relevant_row(row):
        count_relevant = 0
        for col in name_columns:
            cell_content = str(row[col]).lower()
            if any(mcp_term in cell_content for mcp_term in mcp_terms) and not any(false_term in cell_content for false_term in false_terms):
                count_relevant += 1
                return True
        return False

    filtered_df = df[df.apply(is_relevant_row, axis=1)]

    def filter_df_by_strings(df, strings_to_exclude):
        # Create a boolean mask to filter out rows containing any of the strings
        mask = df.apply(lambda row: any(string in row.values for string in strings_to_exclude), axis=1)
        filtered_df = df[~mask]
        return filtered_df

    filtered_df = filter_df_by_strings(filtered_df, false_terms)
    
    return filtered_df



def filter_VIBRANT_genome_quality(VIBRANT_genome_quality, Quality):
    """
    Filters the VIBRANT genome quality based on the specified quality level.

    Args:
        VIBRANT_genome_quality (DataFrame): DataFrame containing VIBRANT genome quality data.
        Quality (str): The desired quality level to filter the data. Valid values are 'low', 'medium', 'high', and 'complete'.

    Returns:
        DataFrame: Filtered DataFrame containing VIBRANT genome quality data.

    """
    quality_args = ['low', 'medium', 'high', 'complete']
    quality_levels = ['low quality draft', 'medium quality draft', 'high quality draft', 'complete circular']
    filtered_genome_quality = VIBRANT_genome_quality[VIBRANT_genome_quality['Quality'].isin(quality_levels[quality_args.index(Quality):])]
    return filtered_genome_quality



def filter_scaffolds_by_contig_length(JGI_global_metadata, min_contig_length=None, max_contig_length=None):
    """
    Filter the JGI global metadata DataFrame based on the length of the contigs
    """
    if min_contig_length is not None and max_contig_length is not None:
        filtered_JGI_global_metadata = JGI_global_metadata[(JGI_global_metadata['contig_length'] >= min_contig_length) & (JGI_global_metadata['contig_length'] <= max_contig_length)]
    elif min_contig_length is not None:
        filtered_JGI_global_metadata = JGI_global_metadata[JGI_global_metadata['contig_length'] >= min_contig_length]
    elif max_contig_length is not None:
        filtered_JGI_global_metadata = JGI_global_metadata[JGI_global_metadata['contig_length'] <= max_contig_length]
    else:
        filtered_JGI_global_metadata = JGI_global_metadata
    
    return filtered_JGI_global_metadata


def write_VIBRANT_contigs_w_mcp_fasta(filtered_VIBRANT_contigs_w_mcp, output23_file_path, output_directory, VIBRANT_input_name,VIBRANT_phages_proteins):
    """
    Writes the contigs with MCP (Major Capsid Protein) fasta file from the filtered VIBRANT contigs.

    Args:
        filtered_VIBRANT_contigs_w_mcp (DataFrame): A DataFrame containing the filtered VIBRANT contigs with MCP.
        output23_file_path (str): The file path of the output23 fasta file.
        output_directory (str): The directory where the output file will be saved.
        VIBRANT_input_name (str): The name of the VIBRANT input.

    Returns:
        None

    Raises:
        FileNotFoundError: If the output23 fasta file does not exist.

    """
    count_genomes_with_MCP = 0
    count_genomes_without_MCP = 0
    contigs_with_mcp_IDs = filtered_VIBRANT_contigs_w_mcp['scaffold'].astype(str).tolist()
    contigs_with_mcp_output_path = output_directory + VIBRANT_input_name + '_contigs_with_mcp.fna'
    contigs_without_mcp_output_path = output_directory + VIBRANT_input_name + '_contigs_without_mcp.fna'
    with open(output23_file_path, "r") as output23_fasta_file:
        genomes_combined = SeqIO.parse(output23_fasta_file, "fasta")
        with open(contigs_with_mcp_output_path, "w") as output_file_with_mcp, open(contigs_without_mcp_output_path, "w") as output_file_without_mcp:
            for genome in genomes_combined:
                if genome.id in contigs_with_mcp_IDs:
                    count_genomes_with_MCP += 1
                    SeqIO.write(genome, output_file_with_mcp, "fasta")
                else:
                    count_genomes_without_MCP += 1
                    SeqIO.write(genome, output_file_without_mcp, "fasta")
    search_identifiers_in_VIBRANT_phages_proteins(contigs_with_mcp_output_path,VIBRANT_phages_proteins)
    print(str(count_genomes_with_MCP) + ' contigs with MCP fasta (.fna) written to:', contigs_with_mcp_output_path)
    search_identifiers_in_VIBRANT_phages_proteins(contigs_without_mcp_output_path,VIBRANT_phages_proteins)
    print(str(count_genomes_without_MCP) + ' contigs without MCP fasta (.fna) written to:', contigs_without_mcp_output_path)


def filter_and_write_genomes_by_quality(VIBRANT_genomes, VIBRANT_genome_quality, min_quality_by_usr, output_directory, VIBRANT_input_name,output23_file_path,VIBRANT_phages_proteins):
    """
    Filters and writes genomes based on their quality.

    Args:
        VIBRANT_genomes (list): List of VIBRANT genomes.
        VIBRANT_genome_quality (DataFrame): DataFrame containing genome quality information.
        min_quality_by_usr (str): Minimum quality threshold specified by the user.
        output_directory (str): Output directory path.
        VIBRANT_input_name (str): Name of the VIBRANT input.
        VIBRANT_phages_proteins (list): List of dictionaries representing the VIBRANT phage protein sequences.

    Returns:
        None
    """
    count_genomes_above_quality = 0
    count_genomes_bellow_quality = 0
    contigs_above_quality_IDs = VIBRANT_genome_quality['scaffold'].astype(str).tolist()
    contigs_bellow_quality_IDs = [genome.id for genome in VIBRANT_genomes if genome.id not in contigs_above_quality_IDs]
    contigs_above_quality_path = output_directory + VIBRANT_input_name + '_above_' + min_quality_by_usr + '_quality_genomes_' + '.fna'
    contigs_bellow_quality_path = output_directory + VIBRANT_input_name + '_bellow_' + min_quality_by_usr + '_quality_genomes_' + '.fna'
    with open(output23_file_path, "r") as output23_fasta_file:
        genomes_combined = SeqIO.parse(output23_fasta_file, "fasta")
        with open(contigs_above_quality_path, "w") as above_quality_output_file, open(contigs_bellow_quality_path, "w") as bellow_quality_output_file:
            for genome in genomes_combined:
                if genome.id in contigs_above_quality_IDs:
                    count_genomes_above_quality += 1
                    SeqIO.write(genome, above_quality_output_file, "fasta")
                else:
                    count_genomes_bellow_quality += 1
                    SeqIO.write(genome, bellow_quality_output_file, "fasta")
                    
        print( str(count_genomes_above_quality) + '(out of ' + str(count_genomes_above_quality + count_genomes_bellow_quality) + ') contigs of at least ' + min_quality_by_usr + ' quality (.fna) written to:', contigs_above_quality_path)
        print( str(count_genomes_bellow_quality) + '(out of ' + str(count_genomes_above_quality + count_genomes_bellow_quality) + ') contigs below ' + min_quality_by_usr + ' quality (.fna) written to:', contigs_bellow_quality_path)
    search_identifiers_in_VIBRANT_phages_proteins(contigs_above_quality_path,VIBRANT_phages_proteins)
    search_identifiers_in_VIBRANT_phages_proteins(contigs_bellow_quality_path,VIBRANT_phages_proteins)


def filter_and_write_MCP_sequences(filtered_VIBRANT_contigs_w_mcp, VIBRANT_phages_proteins, output_directory):
    """
    Filters the VIBRANT phage protein sequences based on a list of MCP IDs,
    writes the filtered sequences to a FASTA file, and returns the filtered sequences.

    Args:
        filtered_VIBRANT_contigs_w_mcp (DataFrame): A DataFrame containing the filtered VIBRANT contigs with MCPs.
        VIBRANT_phages_proteins (list): A list of dictionaries representing the VIBRANT phage protein sequences.
        output_directory (str): The directory where the output file will be saved.
        output_name (str): The name of the output file.

    Returns:
        list: A list of SeqRecord objects representing the filtered MCP sequences.

    """
    MCPs_IDS = filtered_VIBRANT_contigs_w_mcp['protein'].tolist()
    VIBRANT_phages_proteins_ids = [record['id'] for record in VIBRANT_phages_proteins]
    MCPs_sequences = [SeqRecord(Seq(record['sequence']), id=record['id'], description='') for record in VIBRANT_phages_proteins if record['id'] in MCPs_IDS]
    not_MCPs_sequences = [SeqRecord(Seq(record['sequence']), id=record['id'], description='') for record in VIBRANT_phages_proteins if record['id'] not in MCPs_IDS]
    mpcs_output_path = os.path.join(output_directory, "Proteins_identified_as_MCPs.faa")
    not_mpcs_output_path = os.path.join(output_directory, "Proteins_identified_as_not_MCPs.faa")
    SeqIO.write(not_MCPs_sequences, not_mpcs_output_path, "fasta")
    SeqIO.write(MCPs_sequences, mpcs_output_path, "fasta")
    num_mcp_sequences = len(MCPs_sequences)
    print(str(num_mcp_sequences) + " MCPs fasta file created at: " + mpcs_output_path)
    print(str(len(not_MCPs_sequences)) + " non-MCPs fasta file created at: " + not_mpcs_output_path)
    
    return MCPs_sequences


def write_contigs_by_contig_length(VIBRANT_genomes, JGI_global_metadata, filtered_JGI_global_metadata, output_directory, VIBRANT_input_name,output23_file_path,VIBRANT_phages_proteins):
    """
    Write contigs to separate files based on their contig length.

    Args:
        VIBRANT_genomes (list): List of VIBRANT genomes.
        JGI_global_metadata (pandas.DataFrame): DataFrame containing global metadata.
        filtered_JGI_global_metadata (pandas.DataFrame): DataFrame containing filtered global metadata.
        output_directory (str): Directory to save the output files.
        VIBRANT_input_name (str): Name of the VIBRANT input.
        output23_file_path (str): Path to the output23 fasta file.

    Returns:
        None
    """
    count_contigs_within_range = 0
    count_contigs_out_of_range = 0
    all_IDs = JGI_global_metadata['genome_id']
    
    contigs_within_range_IDs = filtered_JGI_global_metadata['genome_id'].astype(str).tolist()
    contigs_out_of_range_IDs = [genome.id for genome in VIBRANT_genomes if genome.id not in contigs_within_range_IDs]
    
    min_contig_length = filtered_JGI_global_metadata['contig_length'].min()
    max_contig_length = filtered_JGI_global_metadata['contig_length'].max()

    df_excluded_contigs = JGI_global_metadata[~JGI_global_metadata['genome_id'].isin(filtered_JGI_global_metadata['genome_id'])]
    min_contig_length_excluded = df_excluded_contigs['contig_length'].min()
    max_contig_length_excluded = df_excluded_contigs['contig_length'].max()
    
    contigs_within_range_path = output_directory + VIBRANT_input_name + '_contigs_within_' + str(min_contig_length) + 'kb_to_' + str(max_contig_length) + 'kb_' + '.fna'
    contigs_out_of_range_path = output_directory + VIBRANT_input_name + '_contigs_within_' + str(min_contig_length_excluded) + 'kb_to_' + str(max_contig_length_excluded) + 'kb_' + '.fna'
    
    with open(output23_file_path, "r") as output23_fasta_file:
        genomes_combined = SeqIO.parse(output23_fasta_file, "fasta")
        with open(contigs_within_range_path, "w") as contigs_within_range_file, open(contigs_out_of_range_path, "w") as contigs_out_of_range_file:
            for genome in genomes_combined:
                if genome.id in contigs_within_range_IDs:
                    count_contigs_within_range += 1
                    SeqIO.write(genome, contigs_within_range_file, "fasta")
                else:
                    count_contigs_out_of_range += 1
                    SeqIO.write(genome, contigs_out_of_range_file, "fasta")
                    
        print( str(count_contigs_within_range) + '(out of initial ' + str(len(all_IDs)) + ') contigs within ' + str(min_contig_length) + 'kb to ' + str(max_contig_length) + 'kb written to:', contigs_within_range_path)
        print( str(count_contigs_out_of_range) + '(out of initial ' + str(len(all_IDs)) + ') contigs NOT within ' + str(min_contig_length) + 'kb to ' + str(max_contig_length) + 'kb written to:', contigs_out_of_range_path)
        search_identifiers_in_VIBRANT_phages_proteins(contigs_within_range_path,VIBRANT_phages_proteins)
        search_identifiers_in_VIBRANT_phages_proteins(contigs_out_of_range_path,VIBRANT_phages_proteins)

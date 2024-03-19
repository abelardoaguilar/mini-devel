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
        for col in name_columns:
            cell_content = str(row[col]).lower()
            if any(mcp_term in cell_content for mcp_term in mcp_terms) and not any(false_term in cell_content for false_term in false_terms):
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

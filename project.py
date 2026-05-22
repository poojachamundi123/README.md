import nltk
nltk.download('punkt_tab')

# !pip install biopython

import pandas as pd
import requests
import gzip
import shutil
import os

def download_disgenet():
    try:
        df = pd.read_csv('/content/curated_gene_disease_associations.tsv', sep='\t')
        df.head()

        high_conf = df[df['score'] >= 0.7]

        # Ensure the 'data' directory exists before saving
        os.makedirs('./data', exist_ok=True)

        # Rename the 'geneSymbol' column to 'gene_symbol' before saving
        # Assuming the original file has a column named 'geneSymbol'
        high_conf = high_conf.rename(columns={'geneSymbol': 'gene_symbol'})
        high_conf = high_conf.rename(columns={'source': 'originalDB'})

        high_conf.to_csv('./data/disgenet_high_conf.csv', index=False)  # Save to the data directory

    except Exception as e:
        print(f"Error downloading or extracting file: {e}")
        return None

    return high_conf

# Example code for downloading GEO datasets
# !pip install GEOparse
import GEOparse
import pandas as pd
import numpy as np

def download_geo_dataset(accession, output_file):
    print(f"Downloading {accession}...")
    gse = GEOparse.get_GEO(geo=accession, destdir="./data")

    # Extract expression data
    expression_dict = {}
    for gsm_name, gsm in gse.gsms.items():
        expression_dict[gsm_name] = gsm.table.iloc[:, 1].values

    # Create expression matrix
    expression_df = pd.DataFrame(expression_dict)

    # Add gene identifiers
    expression_df.index = gse.gsms[list(gse.gsms.keys())[0]].table.iloc[:, 0].values

    # Save to file
    expression_df.to_csv(output_file)
    print(f"Saved {accession} to {output_file}")

    return expression_df

# Download all three required GEO datasets
expression_dfs = {}
expression_dfs['GSE24759'] = download_geo_dataset('GSE24759', './data/GSE24759_expression.csv')
# expression_dfs['GSE6042459'] = download_geo_dataset('GSE6042459', './data/GSE6042459_expression.csv')
expression_dfs['GSE100150'] = download_geo_dataset('GSE100150', './data/GSE100150_expression.csv')

def extract_m92_genes():
    """
    Extract the complete list of 30 genes in the M9.2 erythroid module.

    Returns:
        list: The 30 gene symbols in the M9.2 module
    """
    # Complete list of 30 genes in the M9.2 module
    m92_genes = [
        "ALAS2", "SLC4A1", "BCL2L1", "CA1", "FECH",
        "HBD", "HBM", "AHSP", "EPB42", "GYPA",
        "GYPB", "HEMGN", "SPTA1", "SPTB", "ANK1",
        "SLC25A37", "TMOD1", "SNCA", "TMCC2", "GLRX5",
        "EPB49", "FAM210B", "SELENBP1", "DCAF12", "RAB3IL1",
        "LARGE1", "FBXO7", "TRIB2", "TRIM58", "HEPACAM2"
    ]

    # Save the gene list
    os.makedirs('./data', exist_ok=True)
    pd.DataFrame({"gene_symbol": m92_genes}).to_csv("./data/m92_genes.csv", index=False)

    print(f"Extracted and saved {len(m92_genes)} M9.2 module genes")
    return m92_genes

import requests
import pandas as pd
import io  # Import the io module

def download_string_data(gene_list, species=9606):  # 9606 is human
    session = get_session()
    # Convert gene list to STRING IDs
    genes_str = "%0d".join(gene_list)
    url = f"https://string-db.org/api/tsv/get_string_ids?identifiers={genes_str}&species={species}"

    string_ids_response = session.get(url)

    # Use io.StringIO to create a file-like object from the response text
    string_ids_df = pd.read_csv(io.StringIO(string_ids_response.text), sep='\t')

    # Get interactions between these proteins
    string_ids = string_ids_df['stringId'].tolist()
    ids_str = "%0d".join(string_ids)

    interaction_url = f"https://string-db.org/api/tsv/network?identifiers={ids_str}&species={species}"
    interaction_response = session.get(interaction_url)

    # Use io.StringIO again for the interaction data
    interactions_df = pd.read_csv(io.StringIO(interaction_response.text), sep='\t')
    interactions_df.to_csv("./data/string_interactions.csv", index=False)

    return interactions_df

def get_pubmed_abstracts_for_genes(gene_list, max_per_gene=100):
    """
    Enhanced PubMed abstract retrieval with complete content and improved metadata.

    Args:
        gene_list: A list of gene symbols.
        max_per_gene: The maximum number of abstracts to retrieve per gene.

    Returns:
        A pandas DataFrame containing the complete abstracts with rich metadata.
    """
    from Bio import Entrez, Medline
    import time
    import random  # For sampling if too many results

    Entrez.email = "sidharthareddy114@gmail.com"  # Replace with your email

    all_abstracts = []
    successful_genes = 0

    for gene in gene_list:
        try:
            print(f"Retrieving abstracts for gene: {gene}")

            # Use more precise search terms for higher quality results
            search_term = f"{gene}[Gene] AND human[Organism] AND (\"journal article\"[Publication Type] OR review[Publication Type])"
            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_per_gene*2)  # Get more for sampling
            record = Entrez.read(handle)
            handle.close()

            # Get the list of PMIDs
            ids = record["IdList"]

            if not ids:
                print(f"No papers found for gene {gene}")
                continue

            # If we have too many results, sample diverse papers
            if len(ids) > max_per_gene:
                # Try to get a mix of recent and older papers
                ids = ids[:int(max_per_gene*0.7)] + random.sample(ids[int(max_per_gene*0.7):], int(max_per_gene*0.3))
                ids = ids[:max_per_gene]  # Ensure we don't exceed max_per_gene

            # Fetch the abstracts for these PMIDs in batches of 20
            gene_abstracts = []
            for i in range(0, len(ids), 20):
                batch_ids = ids[i:i+20]

                # Respect NCBI's rate limits
                time.sleep(1)

                handle = Entrez.efetch(db="pubmed", id=",".join(batch_ids),
                                      rettype="medline", retmode="text")
                records = list(Medline.parse(handle))
                handle.close()

                # Process each record with comprehensive metadata
                for record in records:
                    if "AB" in record:  # Only include if abstract is available
                        abstract = record["AB"]
                        title = record.get("TI", "No title")
                        journal = record.get("JT", "Unknown journal")
                        year = record.get("DP", "").split()[0] if "DP" in record else "Unknown year"
                        authors = record.get("AU", ["Unknown"])
                        first_author = authors[0] if authors else "Unknown"
                        keywords = "; ".join(record.get("MH", [])) if "MH" in record else ""

                        # Create a comprehensive text entry with full context
                        text_sections = [
                            f"Title: {title}",
                            f"Gene: {gene}",
                            f"Authors: {', '.join(authors[:3])}{'...' if len(authors) > 3 else ''}",
                            f"Journal: {journal} ({year})",
                            f"Keywords: {keywords[:100]}{'...' if len(keywords) > 100 else ''}",
                            f"Abstract: {abstract}"
                        ]

                        text = "\n".join(text_sections)

                        gene_abstracts.append({
                            "gene_symbol": gene,
                            "abstract": abstract,
                            "text": text,
                            "title": title,
                            "journal": journal,
                            "year": year,
                            "pmid": record.get("PMID", ""),
                            "authors": "; ".join(authors[:5]),
                            "keywords": keywords
                        })

            # Add all abstracts for this gene
            all_abstracts.extend(gene_abstracts)
            if gene_abstracts:
                successful_genes += 1

            print(f"Retrieved {len(gene_abstracts)} high-quality abstracts for {gene}")

        except Exception as e:
            print(f"Error retrieving abstracts for {gene}: {str(e)}")

    result_df = pd.DataFrame(all_abstracts)
    print(f"Total gene coverage: {successful_genes}/{len(gene_list)} genes")
    print(f"Total abstracts retrieved: {len(result_df)}")

    # Add source column for consistency
    result_df['source'] = "PubMed"
    return result_df

def prepare_disgenet_for_finetuning(df):
    """
    Enhanced preparation of DisGeNET data with improved formatting and data quality.

    Args:
        df: DataFrame containing DisGeNET data with gene-disease associations

    Returns:
        DataFrame with enriched, well-formatted entries for fine-tuning
    """
    if df is None or df.empty:
        print("Warning: DisGeNET data is empty or None. Returning empty DataFrame.")
        return pd.DataFrame(columns=["text", "gene_symbol", "abstract", "source"])

    # Verify and prepare required columns
    required_cols = ["gene_symbol", "geneId", "diseaseName", "diseaseId", "score"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Warning: DisGeNET data missing required columns: {missing_cols}")
        for col in missing_cols:
            df[col] = "unknown"

    # Create comprehensive, well-structured text entries
    texts = []
    for _, row in df.iterrows():
        # Create detailed entry with structured sections
        sections = []

        # Core information section
        sections.append(f"Gene Information: {row['gene_symbol']} (ID: {row['geneId']}) is associated with {row['diseaseName']} (Disease ID: {row['diseaseId']}).")

        # Evidence strength section
        sections.append(f"Association Evidence: This gene-disease association has a confidence score of {row['score']} on a scale of 0-1.")

        # Source information section
        if 'originalDB' in df.columns and not pd.isna(row.get('originalDB')):
            sections.append(f"Data Source: This association is documented in {row['originalDB']}.")

        # Publication evidence section
        if 'pmid' in df.columns and not pd.isna(row.get('pmid')):
            sections.append(f"Literature Evidence: This has been reported in scientific literature (PMID: {row['pmid']}).")

        # Create meaningful abstract from disease info (for consistency with PubMed entries)
        if 'diseaseType' in df.columns and 'diseaseName' in df.columns:
            abstract = f"Disease type: {row.get('diseaseType', 'Unknown')}. Disease name: {row['diseaseName']}. "
            if 'diseaseSemanticType' in df.columns and not pd.isna(row.get('diseaseSemanticType')):
                abstract += f"Semantic classification: {row['diseaseSemanticType']}. "
            if 'description' in df.columns and not pd.isna(row.get('description')):
                abstract += f"Description: {row['description']}"
        else:
            abstract = f"Information about {row['diseaseName']} associated with gene {row['gene_symbol']}."

        # Combine sections into well-structured text
        text = "\n".join(sections)

        texts.append({
            "text": text,
            "gene_symbol": row['gene_symbol'],
            "abstract": abstract,
            "source": "DisGeNET"
        })

    result_df = pd.DataFrame(texts)
    print(f"Prepared {len(result_df)} enriched DisGeNET entries with improved formatting")
    return result_df

def create_domain_adaptation_dataset():
    """
    Create a high-quality, balanced dataset for domain adaptation by combining data from multiple sources.
    Includes extensive data quality checks and balancing.

    Returns:
        DataFrame containing the enhanced, balanced dataset for domain adaptation
    """
    print("Creating enhanced domain adaptation dataset...")

    # Ensure data directory exists
    os.makedirs('./data', exist_ok=True)

    # 1. Load and process DisGeNET data
    try:
        disgenet_df = pd.read_csv('./data/disgenet_high_conf.csv')
        print(f"Loaded DisGeNET data: {len(disgenet_df)} entries")
    except FileNotFoundError:
        print("Warning: DisGeNET file not found. Check download_disgenet function.")
        disgenet_df = pd.DataFrame()

    disgenet_texts = prepare_disgenet_for_finetuning(disgenet_df)

    # 2. Load and process M9.2 genes
    try:
        m92_genes = pd.read_csv("./data/m92_genes.csv")['gene_symbol'].tolist()
        print(f"Loaded M9.2 genes: {len(m92_genes)} genes")
    except FileNotFoundError:
        print("Warning: M9.2 genes file not found. Extracting genes...")
        m92_genes = extract_m92_genes()

    # 3. Retrieve PubMed abstracts with improved quality
    print("Fetching high-quality PubMed abstracts for M9.2 genes...")
    pubmed_texts = get_pubmed_abstracts_for_genes(m92_genes, max_per_gene=25)

    # 4. Check coverage of genes
    covered_genes = pubmed_texts['gene_symbol'].unique()
    missing_genes = set(m92_genes) - set(covered_genes)
    if missing_genes:
        print(f"Warning: No PubMed abstracts found for {len(missing_genes)} genes: {missing_genes}")

    # 5. Ensure consistent columns for both data sources
    common_columns = ['text', 'gene_symbol', 'abstract', 'source']
    for col in common_columns:
        if col not in disgenet_texts.columns:
            disgenet_texts[col] = None
        if col not in pubmed_texts.columns:
            pubmed_texts[col] = None

    # 6. Apply data balancing
    disgenet_count = len(disgenet_texts)
    pubmed_count = len(pubmed_texts)

    # Balance the dataset to avoid source bias
    target_count = max(min(disgenet_count, pubmed_count) * 2, 100)  # Ensure minimum size
    print(f"Balancing dataset - DisGeNET: {disgenet_count}, PubMed: {pubmed_count}, Target per source: {target_count//2}")

    # Sample or duplicate entries to balance sources
    if disgenet_count > target_count // 2:
        print(f"Sampling DisGeNET entries from {disgenet_count} to {target_count//2}")
        disgenet_texts = disgenet_texts.sample(n=target_count//2, random_state=42)
    elif disgenet_count < target_count // 2:
        # Duplicate some entries to reach target count
        needed = (target_count//2) - disgenet_count
        print(f"Adding {needed} additional DisGeNET entries through selective duplication")
        additional = disgenet_texts.sample(n=needed, random_state=42, replace=True)
        disgenet_texts = pd.concat([disgenet_texts, additional])

    if pubmed_count > target_count // 2:
        print(f"Sampling PubMed entries from {pubmed_count} to {target_count//2}")
        # Ensure we maintain gene diversity by stratified sampling
        pubmed_texts = pubmed_texts.groupby('gene_symbol', group_keys=False).apply(
            lambda x: x.sample(min(len(x), max(1, int((target_count//2) * len(x) / pubmed_count))), random_state=42)
        ).reset_index(drop=True)
    elif pubmed_count < target_count // 2:
        # Duplicate some entries to reach target count
        needed = (target_count//2) - pubmed_count
        print(f"Adding {needed} additional PubMed entries through selective duplication")
        additional = pubmed_texts.sample(n=needed, random_state=42, replace=True)
        pubmed_texts = pd.concat([pubmed_texts, additional])

    # 7. Combine balanced datasets
    print("Combining balanced data sources...")
    all_texts = pd.concat([disgenet_texts[common_columns], pubmed_texts[common_columns]], ignore_index=True)

    # 8. Quality control on the final dataset
    # Replace empty strings with meaningful values
    all_texts['text'] = all_texts['text'].apply(lambda x: x if x and len(str(x).strip()) > 0 else "No text available")
    all_texts['abstract'] = all_texts['abstract'].apply(lambda x: x if x and len(str(x).strip()) > 0 else "No abstract available")
    all_texts['gene_symbol'] = all_texts['gene_symbol'].apply(lambda x: x if x and len(str(x).strip()) > 0 else "Unknown")

    # 9. Add a quality score column (simple version based on text length)
    all_texts['quality_score'] = all_texts['text'].str.len() / 100
    all_texts['quality_score'] = all_texts['quality_score'].clip(0, 10)  # Scale 0-10

    # 10. Add unique ID for tracking
    all_texts['entry_id'] = [f"{row.source[0]}{i:04d}" for i, row in all_texts.iterrows()]

    # 11. Shuffle the dataset for better training
    all_texts = all_texts.sample(frac=1, random_state=42).reset_index(drop=True)

    # 12. Print comprehensive statistics
    print("\n=== Final Dataset Statistics ===")
    print(f"Total entries: {len(all_texts)}")
    print(f"Source distribution:\n{all_texts['source'].value_counts()}")
    print(f"Gene coverage: {all_texts['gene_symbol'].nunique()} unique genes")
    print(f"Average text length: {all_texts['text'].str.len().mean():.1f} characters")
    print(f"Empty fields: text={sum(all_texts['text']=='No text available')}, "
          f"abstract={sum(all_texts['abstract']=='No abstract available')}, "
          f"gene_symbol={sum(all_texts['gene_symbol']=='Unknown')}")

    # 13. Save the enhanced dataset
    all_texts.to_csv('./data/domain_adaptation_data_enhanced.csv', index=False)
    print("Enhanced domain adaptation dataset saved to ./data/domain_adaptation_data_enhanced.csv")

    # 14. Save a simplified version with just the essential columns
    essential_cols = ['text', 'gene_symbol', 'abstract', 'source']
    all_texts[essential_cols].to_csv('./data/domain_adaptation_data.csv', index=False)
    print("Simplified version saved to ./data/domain_adaptation_data.csv")

    return all_texts

def validate_domain_dataset(filepath='./data/domain_adaptation_data.csv'):
    """
    Validate and analyze a domain adaptation dataset to identify quality issues.
    Provides detailed statistics and highlights potential problems.

    Args:
        filepath: Path to the CSV file containing the domain adaptation dataset

    Returns:
        DataFrame containing the validated dataset and print comprehensive statistics
    """
    try:
        df = pd.read_csv(filepath)
        print(f"Successfully loaded dataset from {filepath}")
    except Exception as e:
        print(f"Error loading dataset: {str(e)}")
        return None

    # Basic dataset properties
    print(f"\n=== Dataset Overview ===")
    print(f"Shape: {df.shape}")
    print(f"Columns: {', '.join(df.columns)}")

    # Check for missing values
    missing = df.isna().sum()
    print(f"\n=== Missing Values ===")
    for col in df.columns:
        print(f"{col}: {missing[col]} ({missing[col]/len(df)*100:.1f}%)")

    # Check for empty strings
    empty = {}
    for col in df.columns:
        if df[col].dtype == 'object':
            empty[col] = (df[col] == '').sum()
    print(f"\n=== Empty Strings ===")
    for col, count in empty.items():
        print(f"{col}: {count} ({count/len(df)*100:.1f}%)")

    # Source distribution
    if 'source' in df.columns:
        print(f"\n=== Source Distribution ===")
        source_counts = df['source'].value_counts()
        for source, count in source_counts.items():
            print(f"{source}: {count} ({count/len(df)*100:.1f}%)")

    # Gene coverage
    if 'gene_symbol' in df.columns:
        unique_genes = df['gene_symbol'].unique()
        print(f"\n=== Gene Coverage ===")
        print(f"Unique genes: {len(unique_genes)}")
        gene_counts = df['gene_symbol'].value_counts()
        print(f"Top 5 genes by frequency:")
        for gene, count in gene_counts.head(5).items():
            print(f"  {gene}: {count} entries")
        print(f"Bottom 5 genes by frequency:")
        for gene, count in gene_counts.tail(5).items():
            print(f"  {gene}: {count} entries")

    # Content length analysis
    if 'text' in df.columns:
        df['text_length'] = df['text'].astype(str).str.len()
        print(f"\n=== Text Length Statistics ===")
        print(f"Min: {df['text_length'].min()}")
        print(f"Max: {df['text_length'].max()}")
        print(f"Mean: {df['text_length'].mean():.1f}")
        print(f"Median: {df['text_length'].median():.1f}")

        # Flag potentially problematic content
        short_texts = df[df['text_length'] < 50]
        print(f"\n{len(short_texts)} entries with very short text (<50 chars)")

        # Sample of content
        print(f"\n=== Content Samples ===")
        print("Random sample of 2 entries:")
        sample = df.sample(min(2, len(df)))
        for i, row in sample.iterrows():
            print(f"\nEntry {i}:")
            print(f"Source: {row.get('source', 'N/A')}")
            print(f"Gene: {row.get('gene_symbol', 'N/A')}")
            print(f"Text snippet: {str(row.get('text', 'N/A'))[:100]}...")

    # Detect duplicates
    duplicates = df.duplicated().sum()
    print(f"\n=== Duplicates ===")
    print(f"Exact duplicates: {duplicates} ({duplicates/len(df)*100:.1f}%)")

    # Return the dataset for further processing
    return df

def gather_evidence_for_gene_criterion(gene, criterion, min_pmids=3, max_pmids=10):
    """
    Gather evidence from literature for a specific gene and evaluation criterion.

    Args:
        gene (str): The gene symbol (e.g., "ALAS2")
        criterion (str): The evaluation criterion (e.g., "erythroid_relevance")
        min_pmids (int): Minimum number of PMIDs to gather
        max_pmids (int): Maximum number of PMIDs to gather

    Returns:
        dict: Evidence including summary, PMIDs, and estimated score
    """
    import requests
    import xml.etree.ElementTree as ET
    import time
    import re
    from nltk.tokenize import sent_tokenize
    session = get_session()

    # Convert criterion to search terms
    criterion_search_terms = {
        "erythroid_relevance": ["erythroid", "erythrocyte", "red blood cell", "erythropoiesis", "hemoglobin"],
        "clinical_biomarker_status": ["biomarker", "clinical marker", "diagnostic", "prognostic", "clinical utility"],
        "potential_biomarker_value": ["potential biomarker", "biomarker candidate", "clinical potential", "diagnostic potential"],
        "leukocyte_biology_relevance": ["leukocyte", "white blood cell", "immune cell", "leukocyte function"],
        "drug_target_potential": ["drug target", "therapeutic target", "pharmacological target", "inhibitor", "drug development"],
        "immune_disease_relevance": ["immune disease", "autoimmune", "inflammation", "immune disorder", "immunodeficiency"]
    }

    search_terms = criterion_search_terms.get(criterion, [criterion.replace("_", " ")])

    # Build search query
    search_query = f"{gene}[gene] AND ({' OR '.join(search_terms)})"

    # Search PubMed
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={search_query}&retmax={max_pmids*2}&sort=relevance"

    try:
        search_response = session.get(search_url)
        search_root = ET.fromstring(search_response.content)

        # Get PMIDs
        pmids = [id_elem.text for id_elem in search_root.findall(".//Id")]

        if not pmids:
            # If no specific results, try broader search
            broader_search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={gene}[gene]&retmax={max_pmids}&sort=relevance"
            broader_response = session.get(broader_search_url)
            broader_root = ET.fromstring(broader_response.content)
            pmids = [id_elem.text for id_elem in broader_root.findall(".//Id")]

        # If still no results, return default values
        if not pmids:
            return {
                "summary": f"Limited evidence found for {gene} related to {criterion}.",
                "pmids": [],
                "abstracts": [],
                "estimated_score": 3  # Default low score due to lack of evidence
            }

        # Fetch abstracts (in batches to respect API limits)
        all_abstracts = []
        batch_size = 5

        for i in range(0, min(len(pmids), max_pmids), batch_size):
            batch_pmids = pmids[i:i+batch_size]
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={','.join(batch_pmids)}&retmode=xml"

            fetch_response = session.get(fetch_url)
            fetch_root = ET.fromstring(fetch_response.content)

            # Process each article
            for article in fetch_root.findall(".//PubmedArticle"):
                article_pmid = article.find(".//PMID").text

                # Get abstract text
                abstract_text = ""
                abstract_elems = article.findall(".//AbstractText")

                if abstract_elems:
                    for elem in abstract_elems:
                        if elem.text:
                            abstract_text += elem.text + " "

                # Get title
                title_elem = article.find(".//ArticleTitle")
                title = title_elem.text if title_elem is not None and title_elem.text else "No title"

                # Get year
                year_elem = article.find(".//PubDate/Year")
                year = year_elem.text if year_elem is not None else "Unknown year"

                if abstract_text:
                    all_abstracts.append({
                        "pmid": article_pmid,
                        "title": title,
                        "year": year,
                        "abstract": abstract_text
                    })

            # Respect NCBI API rate limits
            time.sleep(0.5)

        # Filter abstracts to find the most relevant ones
        relevant_abstracts = []
        for abstract_data in all_abstracts:
            abstract = abstract_data["abstract"].lower()
            title = abstract_data["title"].lower()

            # Check if both gene and criterion terms are mentioned
            gene_mentioned = gene.lower() in abstract or gene.lower() in title
            criterion_mentioned = any(term.lower() in abstract or term.lower() in title for term in search_terms)

            if gene_mentioned and criterion_mentioned:
                relevant_abstracts.append(abstract_data)

        # If no relevant abstracts, use all abstracts
        if not relevant_abstracts:
            relevant_abstracts = all_abstracts

        # Get used PMIDs
        used_pmids = [abstract["pmid"] for abstract in relevant_abstracts[:min(len(relevant_abstracts), max_pmids)]]

        # Generate evidence summary from abstracts
        summary = generate_evidence_summary(gene, criterion, relevant_abstracts)

        # Estimate score based on evidence strength
        estimated_score = estimate_score_from_evidence(gene, criterion, relevant_abstracts)

        return {
            "summary": summary,
            "pmids": used_pmids,
            "abstracts": relevant_abstracts[:max_pmids],
            "estimated_score": estimated_score
        }

    except Exception as e:
        print(f"Error gathering evidence for {gene}/{criterion}: {str(e)}")
        return {
            "summary": f"Error gathering evidence for {gene} related to {criterion}.",
            "pmids": [],
            "abstracts": [],
            "estimated_score": 5  # Default neutral score due to error
        }

def generate_evidence_summary(gene, criterion, abstracts, max_sentences=10):
    """Generate a concise summary from the abstracts for a gene and criterion."""
    if not abstracts:
        return f"No published evidence found for {gene} in relation to {criterion}."

    # Extract relevant sentences from abstracts
    all_sentences = []

    # Convert criterion to keywords for matching
    criterion_keywords = criterion.replace("_", " ").split()

    # Get all sentences that mention the gene or criterion keywords
    for abstract_data in abstracts:
        abstract = abstract_data["abstract"]
        sentences = sent_tokenize(abstract)

        for sentence in sentences:
            if gene in sentence or any(keyword in sentence.lower() for keyword in criterion_keywords):
                clean_sentence = sentence.strip()
                if clean_sentence and len(clean_sentence) > 20:  # Avoid fragments
                    all_sentences.append(clean_sentence)

    # Limit number of sentences
    selected_sentences = all_sentences[:max_sentences]

    # Create a coherent summary
    if selected_sentences:
        evidence_text = " ".join(selected_sentences)
        summary = f"Based on published literature, {gene} shows the following relevance to {criterion.replace('_', ' ')}: {evidence_text}"
    else:
        # Fallback if no good sentences were found
        summary = f"Published literature mentions {gene} in the context of {criterion.replace('_', ' ')}, but specific details are limited."

    return summary

def estimate_score_from_evidence(gene, criterion, abstracts):
    """Estimate a score (0-10) based on the strength of evidence in the abstracts."""
    if not abstracts:
        return 3  # Default low score for no evidence

    # Scoring factors
    factors = {
        "num_abstracts": min(len(abstracts), 10) / 10 * 3,  # 0-3 points based on number of papers
        "recency": 0,
        "directness": 0,
        "consistency": 0
    }

    # Check recency (more recent publications score higher)
    years = []
    for abstract in abstracts:
        try:
            year = int(abstract["year"])
            years.append(year)
        except:
            continue

    if years:
        current_year = 2024
        avg_year = sum(years) / len(years)
        recency_score = min(3, max(0, 3 - (current_year - avg_year) / 5))
        factors["recency"] = recency_score

    # Check directness (how explicitly the abstracts connect gene to criterion)
    criterion_terms = criterion.replace("_", " ").split()
    directness_count = 0

    for abstract in abstracts:
        abstract_text = abstract["abstract"].lower()
        gene_mentioned = gene.lower() in abstract_text
        criterion_mentioned = any(term.lower() in abstract_text for term in criterion_terms)

        # Check for co-occurrence in same sentence
        sentences = sent_tokenize(abstract_text)
        for sentence in sentences:
            if gene.lower() in sentence.lower() and any(term.lower() in sentence.lower() for term in criterion_terms):
                directness_count += 1
                break

    factors["directness"] = min(2, directness_count / 2)

    # Check consistency (do abstracts agree with each other)
    # This is simplified - in practice would need more sophisticated NLP
    factors["consistency"] = 2  # Default assumption of consistency

    # Calculate total score (0-10)
    total_score = sum(factors.values())

    # Round to nearest integer and ensure within range
    return max(0, min(10, round(total_score)))

def generate_gene_statements(gene, n=5, include_fabricated=True):
    """
    Generate statements about a gene with a controlled mixture of factual and non-factual content.

    Args:
        gene (str): The gene symbol (e.g., "ALAS2")
        n (int): Number of statements to generate per category (factual and non-factual)
        include_fabricated (bool): Whether to explicitly include fabricated statements

    Returns:
        dict: Dictionary with 'factual' and 'non_factual' statements
    """
    import random
    import requests
    import time

    # Function to fetch basic gene information from NCBI Gene
    def get_gene_info(gene_symbol):
        session = get_session()
        try:
            # Search for gene ID
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_symbol}[Symbol]+AND+human[Organism]&retmode=json"
            search_response = session.get(search_url)
            search_data = search_response.json()

            if int(search_data.get('esearchresult', {}).get('count', 0)) == 0:
                return None

            gene_id = search_data['esearchresult']['idlist'][0]

            # Get gene summary
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json"
            summary_response = session.get(summary_url)
            summary_data = summary_response.json()

            gene_info = summary_data['result'][gene_id]

            return {
                'gene_id': gene_id,
                'name': gene_info.get('name', ''),
                'description': gene_info.get('description', ''),
                'summary': gene_info.get('summary', ''),
                'chromosome': gene_info.get('chromosome', ''),
                'location': gene_info.get('maplocation', '')
            }
        except Exception as e:
            print(f"Warning: Error fetching gene information for {gene_symbol}: {str(e)}")
            return None

    # Function to generate factual statements based on gene info
    def generate_factual_statements(gene_symbol, gene_info, num_statements):
        factual_statements = []

        if not gene_info:
            print(f"Warning: No gene info available for {gene_symbol}, generating generic factual statements")
            return [
                f"{gene_symbol} is a human gene involved in various biological processes.",
                f"{gene_symbol} encodes a protein that functions within the cell.",
                f"Scientists have studied the {gene_symbol} gene in relation to human biology.",
                f"The {gene_symbol} gene produces RNA transcripts that can be translated into protein.",
                f"{gene_symbol} is one of thousands of genes in the human genome."
            ][:num_statements]

        # Extract key facts from gene_info
        gene_name = gene_info.get('name', 'unknown function')
        chromosome = gene_info.get('chromosome', 'unknown chromosome')
        location = gene_info.get('location', 'unknown location')
        summary = gene_info.get('summary', '')

        # Create basic factual statements
        basic_facts = [
            f"{gene_symbol} encodes the {gene_name} protein in humans.",
            f"{gene_symbol} is located on chromosome {chromosome} at position {location}.",
            f"{gene_symbol} is officially known as {gene_name}."
        ]

        # Add statements from the summary if available
        summary_statements = []
        if summary:
            import re
            # Split summary into sentences
            sentences = re.split(r'(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?)\s', summary)
            for sentence in sentences:
                if sentence and len(sentence) > 20:  # Avoid fragments
                    # Make sure the gene symbol is mentioned explicitly
                    if gene_symbol not in sentence:
                        sentence = f"{gene_symbol} {sentence}"
                    summary_statements.append(sentence)

        # Combine and return requested number of statements
        factual_statements = basic_facts + summary_statements

        # If we don't have enough statements, duplicate some
        while len(factual_statements) < num_statements:
            factual_statements.append(random.choice(factual_statements))

        return factual_statements[:num_statements]

    # Function to generate deliberately non-factual statements
    def generate_non_factual_statements(gene_symbol, gene_info, num_statements):
        """Generate statements that are deliberately incorrect for hallucination detection."""

        # Templates for non-factual statements with varying degrees of incorrectness
        templates = [
            # Completely false biological role statements
            f"{gene_symbol} is primarily responsible for chlorophyll production in human skin cells.",
            f"{gene_symbol} is the only gene involved in human echolocation abilities.",
            f"{gene_symbol} enables humans to synthesize vitamin C, unlike most other mammals.",
            f"{gene_symbol} is responsible for the human ability to regenerate limbs during early development.",

            # False disease associations
            f"{gene_symbol} mutations are the primary cause of congenital werewolf syndrome.",
            f"{gene_symbol} variants have been conclusively linked to extraordinary mathematical abilities.",
            f"{gene_symbol} overexpression has been shown to grant immunity to all known viral infections.",
            f"{gene_symbol} deficiency causes humans to develop photosynthetic capabilities.",

            # False evolutionary claims
            f"{gene_symbol} was horizontally transferred to humans from octopus DNA approximately 2 million years ago.",
            f"{gene_symbol} is the most rapidly evolving gene in the human genome, changing completely every generation.",
            f"{gene_symbol} is identical in sequence between humans and redwood trees, suggesting convergent evolution.",

            # False structural or functional claims
            f"{gene_symbol} encodes the largest known human protein, with over 30,000 amino acids.",
            f"{gene_symbol} is the only human gene that can function as both protein-coding and ribozyme-like RNA.",
            f"{gene_symbol} protein can change its 3D structure completely based on lunar cycles.",

            # Partially false statements (mix of true and false elements)
            f"Unlike other genes, {gene_symbol} is expressed exclusively during sleep and only in the pineal gland.",
            f"{gene_symbol} was first discovered in 1847 by Louis Pasteur during his studies on crystallography.",
            f"Recent CRISPR studies of {gene_symbol} revealed it contains coded messages in Morse code."
        ]

        # If we have gene info, create some plausible but false statements
        if gene_info:
            gene_name = gene_info.get('name', '')
            chromosome = gene_info.get('chromosome', '')

            # Create partially false statements based on true information
            partially_false = [
                f"{gene_symbol} encodes {gene_name}, which was discovered to be the primary enzyme in teleportation processes.",
                f"While {gene_symbol} is located on chromosome {chromosome}, it spontaneously relocates to mitochondrial DNA during cell division.",
                f"{gene_symbol} functions primarily in transmitting memories between generations via epigenetic mechanisms.",
                f"The protein encoded by {gene_symbol} has been shown to reverse the aging process in clinical trials.",
                f"{gene_symbol} expression increases 1000-fold when humans are exposed to zero gravity environments."
            ]

            templates.extend(partially_false)

        # Randomly select the needed number of templates
        selected_templates = random.sample(templates, min(num_statements, len(templates)))

        # If we need more statements than templates, generate some programmatically
        if len(selected_templates) < num_statements:
            # Create templates with false information
            extra_templates = [
                f"{gene_symbol} was named after the {random.choice(['Egyptian', 'Mayan', 'Atlantean', 'Martian'])} {random.choice(['god', 'goddess', 'deity', 'pharaoh'])} of {random.choice(['healing', 'wisdom', 'fertility', 'immortality'])}.",
                f"A rare variant of {gene_symbol} is found in {random.randint(5, 20)}% of {random.choice(['Olympic athletes', 'centenarians', 'opera singers', 'master chess players'])}, granting them superior {random.choice(['endurance', 'memory', 'reflex', 'sensory'])} abilities.",
                f"{gene_symbol} expression is {random.choice(['doubled', 'halved', 'completely silenced', 'maximized'])} during exposure to {random.choice(['full moon light', 'deep sea pressure', 'high altitude', 'specific musical frequencies'])}."
            ]

            # Add as many as needed
            while len(selected_templates) + len(extra_templates) < num_statements:
                extra_templates.append(f"{gene_symbol} has been shown to {random.choice(['communicate', 'synchronize', 'resonate', 'interface'])} with {random.choice(['quantum fields', 'electromagnetic waves', 'cosmic radiation', 'other genes across different organisms'])}.")

            selected_templates.extend(extra_templates[:num_statements - len(selected_templates)])

        return selected_templates[:num_statements]

    # Try to use the API first (if requested)
    api_statements = []
    gene_info = get_gene_info(gene)

    # Fallback to rule-based generation
    factual = generate_factual_statements(gene, gene_info, n)
    non_factual = generate_non_factual_statements(gene, gene_info, n)

    result = {
        'factual': factual,
        'non_factual': non_factual
    }

    return result

import json
from tqdm import tqdm

def create_instruction_response_pairs(gene_list, criteria):
    # This function requires manual curation or semi-automated approach
    pairs = []

    # Example of semi-automated approach for one criterion
    for gene in tqdm(gene_list):
        for criterion in criteria:
            # Query PubMed or other sources
            evidence = gather_evidence_for_gene_criterion(gene, criterion)

            # Generate instruction
            instruction = f"Evaluate the gene {gene} for its {criterion} on a scale of 0-10."

            # Create response (this requires domain expertise)
            # In practice, this would need to be prepared by experts
            response = f"Score: {evidence['estimated_score']}. {evidence['summary']} Evidence includes: {', '.join(evidence['pmids'])}."

            pairs.append({
                "instruction": instruction,
                "response": response,
                "gene": gene,
                "criterion": criterion,
                "score": evidence['estimated_score']
            })

    # Save the data
    with open('./data/instruction_response_pairs.json', 'w') as f:
        json.dump(pairs, f, indent=2)

    return pairs

def fact_check_statement(gene, statement, confidence_threshold=0.7):
    """
    Enhanced fact checking with clearer classification of factual vs non-factual statements.

    Args:
        gene (str): The gene symbol (e.g., "ALAS2")
        statement (str): The statement to fact-check
        confidence_threshold (float): Threshold for confidence in verification (0-1)

    Returns:
        tuple: (is_factual, correction, evidence, confidence_score)
    """
    import re
    from nltk.tokenize import sent_tokenize

    # Check if statement contains obvious non-factual markers
    def contains_non_factual_markers(text):
        non_factual_markers = [
            'teleport', 'immortal', 'psychic', 'supernatural', 'magic', 'divine',
            'werewolf', 'vampire', 'alien', 'photosynthesis in humans',
            'limb regeneration', 'lunar cycle', 'astrology', 'telepathy',
            'quantum consciousness', 'perpetual motion', 'time travel',
            'most rapidly evolving', 'horizontally transferred',
            'chlorophyll production in human', 'echolocation abilities',
            'synthesize vitamin C', 'regenerate limbs', 'immunity to all',
            'photosynthetic capabilities', 'octopus DNA', 'redwood trees',
            'largest known human protein', 'lunar cycles'
        ]

        return any(marker in text.lower() for marker in non_factual_markers)

    # Extract key concepts from the statement
    def extract_key_concepts(text):
        # Simple extraction of biological terms
        bio_terms = [
            'expressed', 'encodes', 'protein', 'enzyme', 'mutations', 'chromosome',
            'transcription', 'translation', 'domain', 'pathway', 'regulation',
            'disease', 'syndrome', 'disorder', 'function', 'structure',
            'binding', 'activation', 'inhibition', 'promoter', 'methylation',
            'phosphorylation', 'interaction', 'motif', 'homology'
        ]

        found_terms = [term for term in bio_terms if term in text.lower()]

        # Extract quoted or specific names
        specific_names = re.findall(r'"([^"]*)"', text)
        specific_names.extend(re.findall(r'([A-Z][A-Za-z0-9\-]+(?: [A-Z][A-Za-z0-9\-]+)*)', text))

        return {'bio_terms': found_terms, 'specific_names': specific_names}

    # Check if statement makes extreme claims
    def has_extreme_claims(text):
        extreme_markers = [
            'only gene', 'most important', 'always', 'never', 'all humans',
            'every cell', 'completely', 'exclusively', 'primary cause',
            'first discovered', 'conclusively', 'revolutionary', 'breakthrough',
            'miracle', 'cure', 'perfect', 'entirely', 'solely responsible',
            '100%', 'unique ability', 'only known', 'exponentially'
        ]

        return any(marker in text.lower() for marker in extreme_markers)

    # Perform the fact check
    try:
        # Quick check for obvious non-factual content
        if contains_non_factual_markers(statement):
            correction = f"The statement contains scientifically impossible claims about {gene}."
            return False, correction, ["Contains scientifically impossible claims"], 0.05

        # Extract key concepts
        concepts = extract_key_concepts(statement)

        # Look for extreme claims
        if has_extreme_claims(statement):
            confidence_score = 0.3  # Suspicious but not definitively false
            correction = f"The statement makes extreme or absolute claims about {gene} that should be treated with skepticism."
            return False, correction, ["Contains extreme claims"], confidence_score

        # For more nuanced checks, we'd normally query databases
        # But for the sake of creating hallucination data, we'll use a simple heuristic

        # Generate a detection confidence score (0-1)
        # Lower score = likely not factual
        confidence_score = 0.5  # Default middle value

        # Penalize for specific names that may be fabricated
        if concepts['specific_names']:
            confidence_score -= 0.1 * min(len(concepts['specific_names']), 3)

        # Reward for biological terms that suggest legitimate content
        if concepts['bio_terms']:
            confidence_score += 0.1 * min(len(concepts['bio_terms']), 3)

        # Penalize for long, complex statements (often hallucinated)
        if len(statement) > 200:
            confidence_score -= 0.1

        # Add some randomness to simulate real-world uncertainty
        import random
        confidence_score += random.uniform(-0.1, 0.1)

        # Clamp to valid range
        confidence_score = max(0.0, min(1.0, confidence_score))

        # Compare with threshold
        is_factual = confidence_score >= confidence_threshold

        if is_factual:
            return True, statement, ["Passed automated fact check"], confidence_score
        else:
            correction = f"Available evidence does not strongly support this statement about {gene}."
            return False, correction, ["Failed automated fact check"], confidence_score

    except Exception as e:
        print(f"Error in fact checking: {str(e)}")
        # Default to uncertain with low confidence when errors occur
        return False, f"Unable to verify this statement about {gene}.", ["Error in verification process"], 0.2

# Add at the top of your script
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Set up a session with retry capability
def get_session():
    session = requests.Session()
    retry = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504]
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session
session = get_session()

def create_hallucination_pairs(gene_list, min_pairs=50, max_per_gene=3):
    """
    Create a robust dataset of hallucination pairs with guaranteed non-factual content.

    Args:
        gene_list: List of gene symbols to process
        min_pairs: Minimum number of hallucination pairs to generate
        max_per_gene: Maximum number of hallucination pairs per gene

    Returns:
        List of hallucination pairs with source genes and evidence
    """
    import json
    import random
    from tqdm import tqdm

    hallucination_pairs = []
    processed_genes = 0

    print(f"Generating hallucination pairs for {len(gene_list)} genes...")

    # Process genes with progress bar
    for gene in tqdm(gene_list):
        # Generate statements with guaranteed non-factual content
        statements = generate_gene_statements(gene, n=max_per_gene+2, include_fabricated=True)

        # Track non-factual statements found for this gene
        gene_hallucinations = []

        # First, use the explicitly non-factual statements
        for statement in statements['non_factual']:
            # Fact check to get correction and confidence
            is_factual, correction, evidence, confidence = fact_check_statement(gene, statement, confidence_threshold=0.7)

            # Since these are explicitly created to be non-factual, if any are marked as factual,
            # we should lower our confidence threshold
            if is_factual:
                # Re-check with lower threshold
                is_factual, correction, evidence, confidence = fact_check_statement(gene, statement, confidence_threshold=0.3)

            # Add the statement if identified as non-factual
            if not is_factual:
                gene_hallucinations.append({
                    "hallucinated_statement": statement,
                    "corrected_statement": correction,
                    "gene": gene,
                    "evidence": evidence,
                    "confidence": confidence,
                    "source": "fabricated"
                })

        # Add factual statements that are misclassified (to add diversity)
        for statement in statements['factual']:
            # Check with strict threshold to find edge cases
            is_factual, correction, evidence, confidence = fact_check_statement(gene, statement, confidence_threshold=0.8)

            if not is_factual and len(gene_hallucinations) < max_per_gene:
                gene_hallucinations.append({
                    "hallucinated_statement": statement,
                    "corrected_statement": correction,
                    "gene": gene,
                    "evidence": evidence,
                    "confidence": confidence,
                    "source": "misclassified"
                })

        # Keep only up to max_per_gene hallucinations per gene
        gene_hallucinations = gene_hallucinations[:max_per_gene]

        # Add to overall collection
        if gene_hallucinations:
            hallucination_pairs.extend(gene_hallucinations)
            processed_genes += 1

        # Break early if we've collected enough examples
        if len(hallucination_pairs) >= min_pairs:
            print(f"Reached {len(hallucination_pairs)} hallucination pairs (minimum {min_pairs})")
            break

    # If we still don't have enough pairs, create some using pattern-based generation
    if len(hallucination_pairs) < min_pairs:
        additional_needed = min_pairs - len(hallucination_pairs)
        print(f"Only found {len(hallucination_pairs)} hallucination pairs, generating {additional_needed} more...")

        # Create additional hallucination pairs using pattern-based generation
        for i in range(additional_needed):
            # Select a random gene
            gene = random.choice(gene_list)

            # Create a clearly false statement
            templates = [
                f"{gene} has been proven to cause spontaneous teleportation in laboratory mice.",
                f"{gene} expression correlates with telepathic abilities in identical twins.",
                f"{gene} contains encrypted messages from ancient civilizations within its non-coding regions.",
                f"The protein encoded by {gene} can shift between solid and liquid states depending on ambient sound frequencies.",
                f"{gene} mutations have been linked to extraordinary abilities like night vision and water breathing."
            ]

            statement = random.choice(templates)
            correction = f"This statement about {gene} contains scientifically impossible claims that are not supported by evidence."

            hallucination_pairs.append({
                "hallucinated_statement": statement,
                "corrected_statement": correction,
                "gene": gene,
                "evidence": ["Generated as explicit hallucination"],
                "confidence": 0.05,
                "source": "guaranteed_false"
            })

    # Print statistics
    print(f"Created {len(hallucination_pairs)} hallucination pairs from {processed_genes} genes")
    print(f"Source distribution: {sum(1 for p in hallucination_pairs if p['source']=='fabricated')} fabricated, "
          f"{sum(1 for p in hallucination_pairs if p['source']=='misclassified')} misclassified, "
          f"{sum(1 for p in hallucination_pairs if p['source']=='guaranteed_false')} guaranteed false")

    # Save the data
    with open('./data/hallucination_pairs.json', 'w') as f:
        json.dump(hallucination_pairs, f, indent=2)

    # Save a simplified CSV version
    import pandas as pd
    df = pd.DataFrame(hallucination_pairs)
    if 'evidence' in df.columns:
        df['evidence'] = df['evidence'].apply(lambda x: ', '.join(x) if isinstance(x, list) else str(x))
    df.to_csv('./data/hallucination_pairs.csv', index=False)
    print(f"Saved hallucination pairs to JSON and CSV formats")

    return hallucination_pairs

!pip install biopython
from nltk.tokenize import sent_tokenize

# Phase 1: Setup and Data Collection
print("Phase 1: Setting up environment and collecting data...")

# Create directories
os.makedirs("./data", exist_ok=True)
os.makedirs("./models", exist_ok=True)
os.makedirs("./results", exist_ok=True)

# Download and process data
disgenet_data = download_disgenet() # Make sure download_disgenet is called and runs without errors
# Check if file was created
if not os.path.exists('./data/disgenet_high_conf.csv'):
    raise FileNotFoundError("disgenet_high_conf.csv not found. Check download_disgenet function.")

m92_genes = extract_m92_genes()


string_interactions = download_string_data(m92_genes)

# Phase 2: Creating Fine-Tuning Datasets
print("Phase 2: Creating fine-tuning datasets...")

domain_data = create_domain_adaptation_dataset()

criteria = [
    "erythroid_relevance",
    "clinical_biomarker_status",
    "potential_biomarker_value",
    "leukocyte_biology_relevance",
    "drug_target_potential",
    "immune_disease_relevance"
]
instruction_data = create_instruction_response_pairs(m92_genes, criteria)

print("Phase 3: Creating hallucination dataset...")
hallucination_data = create_hallucination_pairs(
    m92_genes,
    min_pairs=150,
    max_per_gene=3
)

validate_domain_dataset()

domain_data.head(50)

instruction_data

hallucination_data

"""# Phase 2: Model Fine-Tuning


* Step 1: Setup models for fine-tuning
		1) Deepseek
		2) Qwen
		3) Mistral
*	Step 2: Implement Stage 1 **domain adaptation fine-tuning**
*	Step 3: Implement Stage 2 **task-specific instruction tuning**
*	Step 4: Implement Stage 3 **factuality enhancement with DPO**
*	Step 5: Optimize models for deployment
*	Step 6: Plotting Graphs and diagrams related to Model Training and behavior



"""

# Install required packages
# !pip install -q transformers peft accelerate bitsandbytes wandb trl tiktoken
# # !pip install --upgrade fsspec==2025.3.2
# !pip install datasets==2.16.0
# !pip install fsspec==2025.3.2
# !pip install gcsfs==2025.3.2
import os
import json
import logging
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    Trainer,
    DataCollatorForLanguageModeling
)
from peft import (
    LoraConfig,
    get_peft_model,
    prepare_model_for_kbit_training,
    TaskType,
    PeftModel
)
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
from datasets import load_dataset
import bitsandbytes as bnb
from tqdm import tqdm
import wandb
from typing import Dict, List, Optional, Union, Any

!nvcc --version

from typing import List
import os
import torch
from transformers import TrainingArguments
from peft import LoraConfig, TaskType

class ModelConfig:
    """Optimized configuration for model fine-tuning on limited hardware (T100 GPU)"""

    def __init__(
        self,
        model_name: str,
        tokenizer_name: str = None,
        # Reduced LoRA parameters to save memory
        lora_r: int = 4,
        lora_alpha: int = 8,
        lora_dropout: float = 0.05,
        # Reduced sequence length
        max_seq_length: int = 512,
        target_modules: List[str] = None,
        # Quantization settings - keeping 4-bit for memory efficiency
        use_4bit: bool = True,
        use_nested_quant: bool = False,
        bnb_4bit_compute_dtype: torch.dtype = torch.bfloat16,
        bnb_4bit_quant_type: str = "nf4",
        # May need to disable flash attention if experiencing issues
        use_flash_attention: bool = False,
        gradient_checkpointing: bool = True,
        # Increased gradient accumulation to compensate for smaller batch size
        gradient_accumulation_steps: int = 8,
        output_dir: str = "./output",
        logging_dir: str = "./logs",
        # Reduced epochs to complete faster
        num_train_epochs: int = 1,
        # Smaller batch size
        per_device_train_batch_size: int = 1,
        per_device_eval_batch_size: int = 1,
        # Learning parameters
        learning_rate: float = 2e-4,
        weight_decay: float = 0.01,
        warmup_ratio: float = 0.03,
        optim: str = "paged_adamw_32bit",
        lr_scheduler_type: str = "cosine",
        # Reduced checkpoint saving to save disk space
        save_total_limit: int = 1,
        save_strategy: str = "steps",
        # Less frequent saving and evaluation
        save_steps: int = 1000,
        logging_steps: int = 200,
        evaluation_strategy: str = "steps",
        eval_steps: int = 1000,
        load_in_8bit: bool = False,
        # Precision settings - adjust based on your specific T100 capabilities
        fp16: bool = True,
        bf16: bool = False,
        seed: int = 42
    ):
        self.model_name = model_name
        self.tokenizer_name = tokenizer_name or model_name
        self.lora_r = lora_r
        self.lora_alpha = lora_alpha
        self.lora_dropout = lora_dropout
        self.max_seq_length = max_seq_length
        self.target_modules = target_modules
        self.use_4bit = use_4bit
        self.use_nested_quant = use_nested_quant
        self.bnb_4bit_compute_dtype = bnb_4bit_compute_dtype
        self.bnb_4bit_quant_type = bnb_4bit_quant_type
        self.use_flash_attention = use_flash_attention
        self.gradient_checkpointing = gradient_checkpointing
        self.gradient_accumulation_steps = gradient_accumulation_steps
        self.output_dir = output_dir
        self.logging_dir = logging_dir
        self.num_train_epochs = num_train_epochs
        self.per_device_train_batch_size = per_device_train_batch_size
        self.per_device_eval_batch_size = per_device_eval_batch_size
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay
        self.warmup_ratio = warmup_ratio
        self.optim = optim
        self.lr_scheduler_type = lr_scheduler_type
        self.save_total_limit = save_total_limit
        self.save_strategy = save_strategy
        self.save_steps = save_steps
        self.logging_steps = logging_steps
        self.evaluation_strategy = evaluation_strategy
        self.eval_steps = eval_steps
        self.load_in_8bit = load_in_8bit
        self.fp16 = fp16
        self.bf16 = bf16
        self.seed = seed

    def get_training_args(self, stage_name: str = None):
        """Get training arguments for Trainer"""
        output_dir = self.output_dir
        if stage_name:
            output_dir = os.path.join(output_dir, stage_name)

        return TrainingArguments(
            output_dir=output_dir,
            num_train_epochs=self.num_train_epochs,
            per_device_train_batch_size=self.per_device_train_batch_size,
            per_device_eval_batch_size=self.per_device_eval_batch_size,
            gradient_accumulation_steps=self.gradient_accumulation_steps,
            gradient_checkpointing=self.gradient_checkpointing,
            learning_rate=self.learning_rate,
            weight_decay=self.weight_decay,
            warmup_ratio=self.warmup_ratio,
            optim=self.optim,
            lr_scheduler_type=self.lr_scheduler_type,
            save_total_limit=self.save_total_limit,
            save_strategy=self.save_strategy,
            save_steps=self.save_steps,
            logging_steps=self.logging_steps,
            evaluation_strategy=self.evaluation_strategy,
            eval_steps=self.eval_steps,
            fp16=self.fp16,
            bf16=self.bf16,
            seed=self.seed,
            report_to="wandb"
        )

    def get_lora_config(self):
        """Get LoRA configuration for PEFT"""
        return LoraConfig(
            r=self.lora_r,
            lora_alpha=self.lora_alpha,
            lora_dropout=self.lora_dropout,
            target_modules=self.target_modules,
            bias="none",
            task_type=TaskType.CAUSAL_LM
        )

    def estimate_memory_usage(self, model_size_gb=None):
        """Estimate memory usage based on configured parameters"""
        # Basic estimation - adjust formula as needed based on empirical testing
        base_memory = 0.5  # Base overhead in GB

        # Model memory (with quantization applied)
        if model_size_gb:
            model_memory = model_size_gb * (0.25 if self.use_4bit else 0.5 if self.load_in_8bit else 1.0)
        else:
            model_memory = 0  # Unknown

        # Batch memory
        seq_factor = self.max_seq_length / 1024  # Normalize to 1k tokens
        batch_memory = seq_factor * self.per_device_train_batch_size * 0.4  # Estimated GB per normalized batch

        # Additional overhead for gradients etc.
        if not self.gradient_checkpointing:
            batch_memory *= 1.5

        total_estimate = base_memory + model_memory + batch_memory

        return {
            "estimated_total_gb": total_estimate,
            "model_memory_gb": model_memory,
            "batch_memory_gb": batch_memory,
            "base_overhead_gb": base_memory
        }

class DeepseekConfig(ModelConfig):
    """Configuration for Deepseek model fine-tuning"""

    def __init__(self, **kwargs):
        default_model = "deepseek-ai/deepseek-llm-7b-base"
        default_target_modules = ["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj", "up_proj", "down_proj"]

        kwargs.setdefault("model_name", default_model)
        kwargs.setdefault("target_modules", default_target_modules)
        super().__init__(**kwargs)

        # Deepseek-specific settings
        self.output_dir = os.path.join(self.output_dir, "deepseek")


class QwenConfig(ModelConfig):
    """Configuration for Qwen model fine-tuning"""

    def __init__(self, **kwargs):
        default_model = "Qwen/Qwen-7B"
        default_target_modules = ["c_attn", "c_proj", "w1", "w2"]

        kwargs.setdefault("model_name", default_model)
        kwargs.setdefault("target_modules", default_target_modules)
        super().__init__(**kwargs)

        # Qwen-specific settings
        self.output_dir = os.path.join(self.output_dir, "qwen")


class MistralConfig(ModelConfig):
    """Configuration for Mistral model fine-tuning"""

    def __init__(self, **kwargs):
        default_model = "mistralai/Mistral-7B-v0.1"
        default_target_modules = ["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj", "up_proj", "down_proj"]

        kwargs.setdefault("model_name", default_model)
        kwargs.setdefault("target_modules", default_target_modules)
        super().__init__(**kwargs)

        # Mistral-specific settings
        self.output_dir = os.path.join(self.output_dir, "mistral")

class ModelManager:
    """Manager for loading and preparing models for fine-tuning"""

    def __init__(self, config: ModelConfig):
        self.config = config
        self.model = None
        self.tokenizer = None
        self.peft_model = None

    def load_model_and_tokenizer(self):
        """Load base model and tokenizer"""
        print(f"Loading model: {self.config.model_name}")

        # Load tokenizer
        self.tokenizer = AutoTokenizer.from_pretrained(
            self.config.tokenizer_name,
            use_fast=True,
            trust_remote_code=True
        )

        # Set padding token if not set
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
            self.tokenizer.pad_token_id = self.tokenizer.eos_token_id

        # Ensure special tokens are set properly
        # For chat models, we might need specific tokens
        self.ensure_chat_tokens()

        # Load model with quantization if specified
        if self.config.use_4bit:
            self.model = AutoModelForCausalLM.from_pretrained(
                self.config.model_name,
                quantization_config=BitsAndBytesConfig(
                    load_in_4bit=True,
                    bnb_4bit_compute_dtype=self.config.bnb_4bit_compute_dtype,
                    bnb_4bit_quant_type=self.config.bnb_4bit_quant_type,
                    bnb_4bit_use_double_quant=self.config.use_nested_quant
                ),
                torch_dtype=self.config.bnb_4bit_compute_dtype,
                trust_remote_code=True,
                use_flash_attention_2=self.config.use_flash_attention
            )
        elif self.config.load_in_8bit:
            self.model = AutoModelForCausalLM.from_pretrained(
                self.config.model_name,
                load_in_8bit=True,
                torch_dtype=torch.float16,
                trust_remote_code=True,
                use_flash_attention_2=self.config.use_flash_attention
            )
        else:
            self.model = AutoModelForCausalLM.from_pretrained(
                self.config.model_name,
                torch_dtype=torch.bfloat16,
                trust_remote_code=True,
                use_flash_attention_2=self.config.use_flash_attention
            )

        # Prepare model for kbit training if needed
        if self.config.use_4bit or self.config.load_in_8bit:
            self.model = prepare_model_for_kbit_training(
                self.model,
                use_gradient_checkpointing=self.config.gradient_checkpointing
            )

        # Enable gradient checkpointing for memory efficiency
        if self.config.gradient_checkpointing:
            self.model.gradient_checkpointing_enable()

        return self.model, self.tokenizer

    def prepare_for_peft(self):
        """Prepare model for PEFT/LoRA fine-tuning"""
        if self.model is None:
            self.load_model_and_tokenizer()

        print("Preparing model for PEFT/LoRA fine-tuning")
        lora_config = self.config.get_lora_config()
        self.peft_model = get_peft_model(self.model, lora_config)
        self.peft_model.print_trainable_parameters()

        return self.peft_model

    def ensure_chat_tokens(self):
        """Ensure chat-specific tokens are properly set"""
        # Different models might have different chat formats
        model_type = self.identify_model_type()

        if model_type == "deepseek":
            self.tokenizer.chat_template = "{% for message in messages %}\n{% if message.role == 'user' %}{{message.content}}\n{% elif message.role == 'assistant' %}{{message.content}}\n{% endif %}{% endfor %}"
        elif model_type == "qwen":
            # Handle Qwen specific tokens
            if not hasattr(self.tokenizer, 'chat_template'):
                self.tokenizer.chat_template = "{% for message in messages %}\n{% if message['role'] == 'user' %}\n{{ '<|im_start|>user\n' + message['content'] + '<|im_end|>\n' }}\n{% elif message['role'] == 'assistant' %}\n{{ '<|im_start|>assistant\n' + message['content'] + '<|im_end|>\n' }}\n{% elif message['role'] == 'system' %}\n{{ '<|im_start|>system\n' + message['content'] + '<|im_end|>\n' }}\n{% endif %}\n{% endfor %}\n{% if add_generation_prompt %}\n{{ '<|im_start|>assistant\n' }}\n{% endif %}"
        elif model_type == "mistral":
            # Handle Mistral specific tokens
            if not hasattr(self.tokenizer, 'chat_template'):
                self.tokenizer.chat_template = "{% for message in messages %}\n{% if message['role'] == 'user' %}\n{{ '<|user|>\n' + message['content'] + '<|assistant|>\n' }}\n{% elif message['role'] == 'assistant' %}\n{{ message['content'] + '<|endoftext|>\n' }}\n{% elif message['role'] == 'system' %}\n{{ '<|system|>\n' + message['content'] + '<|user|>\n' }}\n{% endif %}\n{% endfor %}"

    def identify_model_type(self):
        """Identify the type of model based on model name"""
        model_name = self.config.model_name.lower()
        if "deepseek" in model_name:
            return "deepseek"
        elif "qwen" in model_name:
            return "qwen"
        elif "mistral" in model_name:
            return "mistral"
        else:
            return "unknown"

from transformers import BitsAndBytesConfig

class DatasetPreparation:
    """Utilities for preparing datasets for different training stages"""

    def __init__(self, tokenizer, max_seq_length=2048):
        self.tokenizer = tokenizer
        self.max_seq_length = max_seq_length

    def format_chat_message(self, example):
        """Format messages into chat template compatible with tokenizer"""
        messages = example["messages"]
        return {"text": self.tokenizer.apply_chat_template(
            messages,
            tokenize=False,
            add_generation_prompt=False
        )}

    def preprocess_function(self, examples):
        """Tokenize text examples"""
        batch_size = len(examples["text"])
        inputs = [f"{text}" for text in examples["text"]]
        model_inputs = self.tokenizer(
            inputs,
            max_length=self.max_seq_length,
            padding="max_length",
            truncation=True,
            return_tensors="pt"
        )

        # Create labels (needed for causal language modeling)
        labels = model_inputs["input_ids"].clone()

        # Mask padding tokens
        for i in range(batch_size):
            labels[i][model_inputs["attention_mask"][i] == 0] = -100

        model_inputs["labels"] = labels
        return model_inputs

    def load_domain_adaptation_dataset(self, file_path, val_split=0.1):
        """Load domain adaptation dataset"""
        print(f"Loading domain adaptation dataset from {file_path}")

        # Read the dataset
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Create formatted dataset
        formatted_data = {
            "text": [],
        }

        for item in data:
            # Format according to domain adaptation format
            formatted_item = self.format_chat_message(item)
            formatted_data["text"].append(formatted_item["text"])

        # Split into train and validation
        total_samples = len(formatted_data["text"])
        val_size = int(total_samples * val_split)

        train_data = {
            "text": formatted_data["text"][val_size:]
        }

        val_data = {
            "text": formatted_data["text"][:val_size]
        }

        # Convert to dataset objects
        train_dataset = Dataset.from_dict(train_data)
        val_dataset = Dataset.from_dict(val_data)

        # Tokenize
        train_dataset = train_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        val_dataset = val_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        return train_dataset, val_dataset

    def load_instruction_dataset(self, file_path, val_split=0.1):
        """Load instruction dataset"""
        print(f"Loading instruction dataset from {file_path}")

        # Read the dataset
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Create formatted dataset
        formatted_data = {
            "text": [],
        }

        for item in data:
            # Format according to instruction format
            formatted_item = self.format_chat_message(item)
            formatted_data["text"].append(formatted_item["text"])

        # Split into train and validation
        total_samples = len(formatted_data["text"])
        val_size = int(total_samples * val_split)

        train_data = {
            "text": formatted_data["text"][val_size:]
        }

        val_data = {
            "text": formatted_data["text"][:val_size]
        }

        # Convert to dataset objects
        train_dataset = Dataset.from_dict(train_data)
        val_dataset = Dataset.from_dict(val_data)

        # Tokenize
        train_dataset = train_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        val_dataset = val_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        return train_dataset, val_dataset

    def load_hallucination_dataset(self, file_path, val_split=0.1):
        """Load hallucination dataset"""
        print(f"Loading hallucination dataset from {file_path}")

        # Read the dataset
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Create formatted dataset
        formatted_data = {
            "text": [],
        }

        for item in data:
            # Format according to hallucination format
            formatted_item = self.format_chat_message(item)
            formatted_data["text"].append(formatted_item["text"])

        # Split into train and validation
        total_samples = len(formatted_data["text"])
        val_size = int(total_samples * val_split)

        train_data = {
            "text": formatted_data["text"][val_size:]
        }

        val_data = {
            "text": formatted_data["text"][:val_size]
        }

        # Convert to dataset objects
        train_dataset = Dataset.from_dict(train_data)
        val_dataset = Dataset.from_dict(val_data)

        # Tokenize
        train_dataset = train_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        val_dataset = val_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        return train_dataset, val_dataset

from datasets import Dataset

class TrainingManager:
    """Manager for handling different stages of model training"""

    def __init__(self, model_config, output_dir="./output"):
        self.model_config = model_config
        self.output_dir = os.path.join(output_dir, model_config.__class__.__name__.replace("Config", "").lower())
        os.makedirs(self.output_dir, exist_ok=True)

        self.model_manager = ModelManager(model_config)
        self.model, self.tokenizer = self.model_manager.load_model_and_tokenizer()
        self.peft_model = None

        self.dataset_prep = DatasetPreparation(
            self.tokenizer,
            max_seq_length=model_config.max_seq_length
        )

        # Initialize logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(f"{self.output_dir}/training.log"),
                logging.StreamHandler()
            ]
        )

        self.logger = logging.getLogger(__name__)

    def initialize_wandb(self, stage_name):
        """Initialize Weights & Biases logging"""
        model_name = self.model_config.model_name.split("/")[-1]
        wandb.init(
            project="genetic-llm-fine-tuning",
            name=f"{model_name}-{stage_name}",
            config={
                "model_name": self.model_config.model_name,
                "stage": stage_name,
                "lora_r": self.model_config.lora_r,
                "lora_alpha": self.model_config.lora_alpha,
                "lora_dropout": self.model_config.lora_dropout,
                "max_seq_length": self.model_config.max_seq_length,
                "learning_rate": self.model_config.learning_rate,
                "num_train_epochs": self.model_config.num_train_epochs,
                "per_device_train_batch_size": self.model_config.per_device_train_batch_size,
                "gradient_accumulation_steps": self.model_config.gradient_accumulation_steps
            }
        )

    def domain_adaptation_training(self, domain_data_path):
        """Stage 1: Domain Adaptation Fine-Tuning"""
        stage_name = "domain_adaptation"
        self.logger.info(f"Starting {stage_name} stage")

        # Initialize W&B
        self.initialize_wandb(stage_name)

        # Prepare model for PEFT if not already done
        if self.peft_model is None:
            self.peft_model = self.model_manager.prepare_for_peft()

        # Load datasets
        train_dataset, val_dataset = self.dataset_prep.load_domain_adaptation_dataset(
            domain_data_path
        )

        # Get training arguments
        training_args = self.model_config.get_training_args(stage_name)

        # Initialize data collator
        data_collator = DataCollatorForLanguageModeling(
            tokenizer=self.tokenizer,
            mlm=False
        )

        # Initialize trainer
        trainer = Trainer(
            model=self.peft_model,
            args=training_args,
            train_dataset=train_dataset,
            eval_dataset=val_dataset,
            data_collator=data_collator,
            tokenizer=self.tokenizer
        )

        # Start training
        self.logger.info("Starting training...")
        train_result = trainer.train()

        # Save model
        self.logger.info("Saving model...")
        trainer.save_model(f"{self.output_dir}/{stage_name}")

        # Save training metrics
        metrics = train_result.metrics
        trainer.log_metrics("train", metrics)
        trainer.save_metrics("train", metrics)

        # Evaluate
        self.logger.info("Evaluating model...")
        eval_results = trainer.evaluate()
        trainer.log_metrics("eval", eval_results)
        trainer.save_metrics("eval", eval_results)

        # Close W&B
        wandb.finish()

        return f"{self.output_dir}/{stage_name}"

    def instruction_tuning(self, instruction_data_path, base_model_path=None):
        """Stage 2: Task-Specific Instruction Tuning"""
        stage_name = "instruction_tuning"
        self.logger.info(f"Starting {stage_name} stage")

        # Initialize W&B
        self.initialize_wandb(stage_name)

        # Load base model from previous stage if provided
        if base_model_path:
            self.logger.info(f"Loading base model from {base_model_path}")
            self.peft_model = PeftModel.from_pretrained(
                self.model,
                base_model_path,
                is_trainable=True
            )
        else:
            # Prepare model for PEFT if not already done
            if self.peft_model is None:
                self.peft_model = self.model_manager.prepare_for_peft()

        # Load datasets
        train_dataset, val_dataset = self.dataset_prep.load_instruction_dataset(
            instruction_data_path
        )

        # Get training arguments
        training_args = self.model_config.get_training_args(stage_name)

        # Initialize data collator
        data_collator = DataCollatorForLanguageModeling(
            tokenizer=self.tokenizer,
            mlm=False
        )

        # Initialize trainer
        trainer = Trainer(
            model=self.peft_model,
            args=training_args,
            train_dataset=train_dataset,
            eval_dataset=val_dataset,
            data_collator=data_collator,
            tokenizer=self.tokenizer
        )

        # Start training
        self.logger.info("Starting training...")
        train_result = trainer.train()

        # Save model
        self.logger.info("Saving model...")
        trainer.save_model(f"{self.output_dir}/{stage_name}")

        # Save training metrics
        metrics = train_result.metrics
        trainer.log_metrics("train", metrics)
        trainer.save_metrics("train", metrics)

        # Evaluate
        self.logger.info("Evaluating model...")
        eval_results = trainer.evaluate()
        trainer.log_metrics("eval", eval_results)
        trainer.save_metrics("eval", eval_results)

        # Close W&B
        wandb.finish()

        return f"{self.output_dir}/{stage_name}"

    def factuality_enhancement_dpo(self, hallucination_data_path, base_model_path=None):
        """Stage 3: Factuality Enhancement with DPO"""
        stage_name = "factuality_dpo"
        self.logger.info(f"Starting {stage_name} stage")

        try:
            from trl import DPOTrainer
        except ImportError:
            print("Installing TRL for DPO training...")
            !pip install -q trl
            from trl import DPOTrainer

        # Initialize W&B
        self.initialize_wandb(stage_name)

        # Load base model from previous stage if provided
        if base_model_path:
            self.logger.info(f"Loading base model from {base_model_path}")
            self.peft_model = PeftModel.from_pretrained(
                self.model,
                base_model_path,
                is_trainable=True
            )
        else:
            # Prepare model for PEFT if not already done
            if self.peft_model is None:
                self.peft_model = self.model_manager.prepare_for_peft()

        # Load and prepare the hallucination dataset for DPO
        with open(hallucination_data_path, 'r') as f:
            hallucination_data = json.load(f)

        # Process the data into DPO format
        dpo_data = {
            "prompt": [],
            "chosen": [],
            "rejected": []
        }

        for item in hallucination_data:
            # For each hallucination item, we need:
            # - prompt: the user query
            # - chosen: the corrected response (factual)
            # - rejected: the hallucinated response (non-factual)
            user_message = next(msg for msg in item["messages"] if msg["role"] == "user")
            assistant_message = next(msg for msg in item["messages"] if msg["role"] == "assistant")

            hallucinated_statement = user_message["content"]

            corrected_statement = assistant_message["content"]

            # Assuming each item has hallucinated_statement, corrected_statement, gene and evidence
            prompt = f"Please review the following statement about the gene {item['gene']}: {hallucinated_statement}"

            # Create chosen and rejected responses
            chosen = f"The statement contains inaccurate information. {corrected_statement}"
            rejected = f"The statement appears to be accurate. {hallucinated_statement} is scientifically valid."

            dpo_data["prompt"].append(prompt)
            dpo_data["chosen"].append(chosen)
            dpo_data["rejected"].append(rejected)

        # Convert to dataset
        dpo_dataset = Dataset.from_dict(dpo_data)

        # Split into train and eval
        dpo_dataset = dpo_dataset.train_test_split(test_size=0.1)

        # Configure DPO training arguments
        dpo_training_args = TrainingArguments(
            output_dir=f"{self.output_dir}/{stage_name}",
            num_train_epochs=self.model_config.num_train_epochs,
            per_device_train_batch_size=self.model_config.per_device_train_batch_size//2, # Reduced batch size due to paired samples
            gradient_accumulation_steps=self.model_config.gradient_accumulation_steps*2, # Increased to compensate for smaller batch
            gradient_checkpointing=self.model_config.gradient_checkpointing,
            learning_rate=self.model_config.learning_rate/2, # Lower learning rate for DPO
            weight_decay=self.model_config.weight_decay,
            lr_scheduler_type=self.model_config.lr_scheduler_type,
            save_strategy=self.model_config.save_strategy,
            save_steps=self.model_config.save_steps,
            logging_steps=self.model_config.logging_steps,
            evaluation_strategy=self.model_config.evaluation_strategy,
            eval_steps=self.model_config.eval_steps,
            fp16=self.model_config.fp16,
            bf16=self.model_config.bf16,
            report_to="wandb"
        )

        # Initialize DPO trainer
        dpo_trainer = DPOTrainer(
            model=self.peft_model,
            args=dpo_training_args,
            train_dataset=dpo_dataset["train"],
            eval_dataset=dpo_dataset["test"],
            # tokenizer=self.tokenizer,
            # beta=0.1,  # DPO hyperparameter
            # max_length=self.model_config.max_seq_length,
            # max_prompt_length=self.model_config.max_seq_length // 2
        )

        # Start DPO training
        self.logger.info("Starting DPO training...")
        dpo_trainer.train()

        # Save model
        self.logger.info("Saving model...")
        dpo_trainer.save_model(f"{self.output_dir}/{stage_name}")

        # Close W&B
        wandb.finish()

        return f"{self.output_dir}/{stage_name}"

    def optimize_for_deployment(self, model_path):
        """Stage 5: Optimize models for deployment"""
        stage_name = "optimized"
        self.logger.info(f"Starting {stage_name} stage")

        # Load the fine-tuned PEFT model
        self.logger.info(f"Loading fine-tuned model from {model_path}")
        fine_tuned_model = PeftModel.from_pretrained(
            self.model,
            model_path
        )

        # Merge adapter weights with base model for efficiency
        self.logger.info("Merging adapter weights with base model...")
        merged_model = fine_tuned_model.merge_and_unload()

        # Save the merged model
        merged_model_path = f"{self.output_dir}/{stage_name}"
        self.logger.info(f"Saving merged model to {merged_model_path}")
        merged_model.save_pretrained(merged_model_path)
        self.tokenizer.save_pretrained(merged_model_path)

        # Optional: Further quantize the model for deployment
        self.logger.info("Quantizing model for deployment...")
        try:
            # You might want to use a different quantization approach
            # depending on your deployment environment
            from optimum.onnxruntime import ORTModelForCausalLM

            # Convert to ONNX format
            ort_model = ORTModelForCausalLM.from_pretrained(
                merged_model_path,
                export=True
            )

            # Save optimized model
            optimized_path = f"{self.output_dir}/{stage_name}_onnx"
            ort_model.save_pretrained(optimized_path)
            self.logger.info(f"Saved ONNX model to {optimized_path}")

            return optimized_path
        except:
            self.logger.warning("ONNX conversion failed, returning merged model instead.")
            return merged_model_path


class ModelEvaluator:
    """Utility for evaluating model performance across different stages"""

    def __init__(self, model_manager, output_dir="./evaluation"):
        self.model_manager = model_manager
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        self.logger = logging.getLogger(__name__)

    def evaluate_model(self, model_path, test_data_path, stage_name=None):
        """Evaluate model on test dataset"""
        self.logger.info(f"Evaluating model from {model_path} on {test_data_path}")

        # Load the model
        model, tokenizer = self.model_manager.model, self.model_manager.tokenizer
        peft_model = PeftModel.from_pretrained(
            model,
            model_path
        )

        # Load the test data
        with open(test_data_path, 'r') as f:
            test_data = json.load(f)

        # Generate responses and evaluate
        results = {
            "input": [],
            "reference": [],
            "prediction": [],
            "is_correct": []
        }

        for item in tqdm(test_data, desc="Evaluating"):
            # Extract input and reference
            input_text = item["messages"][0]["content"]
            reference = item["messages"][1]["content"]

            # Generate prediction
            inputs = tokenizer(
                input_text,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=self.model_manager.config.max_seq_length
            ).to(peft_model.device)

            # Generate output
            with torch.no_grad():
                outputs = peft_model.generate(
                    input_ids=inputs["input_ids"],
                    attention_mask=inputs["attention_mask"],
                    max_new_tokens=512,
                    do_sample=False,
                    num_beams=1
                )

            # Decode prediction
            prediction = tokenizer.decode(outputs[0], skip_special_tokens=True)

            # Simple exact match evaluation (this should be improved with proper metrics)
            is_correct = prediction.strip() == reference.strip()

            # Store results
            results["input"].append(input_text)
            results["reference"].append(reference)
            results["prediction"].append(prediction)
            results["is_correct"].append(is_correct)

        # Calculate overall accuracy
        accuracy = sum(results["is_correct"]) / len(results["is_correct"])

        # Save results
        output_file = f"{self.output_dir}/results_{stage_name or 'evaluation'}.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)

        self.logger.info(f"Evaluation complete. Accuracy: {accuracy:.4f}")
        self.logger.info(f"Results saved to {output_file}")

        return accuracy, results

    def plot_training_metrics(self, metrics_paths, stage_names=None):
        """Plot training metrics from different stages"""
        metrics_data = []

        # Load metrics from each stage
        for i, path in enumerate(metrics_paths):
            stage = stage_names[i] if stage_names and i < len(stage_names) else f"Stage {i+1}"

            try:
                # Load training metrics
                with open(f"{path}/trainer_state.json", 'r') as f:
                    train_data = json.load(f)

                # Extract metrics
                log_history = train_data.get("log_history", [])

                for entry in log_history:
                    if "loss" in entry:
                        metrics_data.append({
                            "stage": stage,
                            "step": entry.get("step", 0),
                            "loss": entry.get("loss", 0),
                            "learning_rate": entry.get("learning_rate", 0),
                            "epoch": entry.get("epoch", 0),
                            "type": "train"
                        })
                    elif "eval_loss" in entry:
                        metrics_data.append({
                            "stage": stage,
                            "step": entry.get("step", 0),
                            "loss": entry.get("eval_loss", 0),
                            "epoch": entry.get("epoch", 0),
                            "type": "eval"
                        })
            except Exception as e:
                self.logger.warning(f"Error loading metrics from {path}: {e}")

        # Convert to DataFrame
        df = pd.DataFrame(metrics_data)

        if df.empty:
            self.logger.warning("No metrics data to plot.")
            return

        # Create plots
        plt.figure(figsize=(12, 10))

        # Plot training loss
        plt.subplot(2, 1, 1)
        for stage in df["stage"].unique():
            stage_data = df[(df["stage"] == stage) & (df["type"] == "train")]
            if not stage_data.empty:
                plt.plot(stage_data["step"], stage_data["loss"], label=f"{stage} Training")

            stage_eval_data = df[(df["stage"] == stage) & (df["type"] == "eval")]
            if not stage_eval_data.empty:
                plt.plot(stage_eval_data["step"], stage_eval_data["loss"], label=f"{stage} Evaluation", linestyle="--")

        plt.xlabel("Training Steps")
        plt.ylabel("Loss")
        plt.title("Training and Evaluation Loss")
        plt.legend()
        plt.grid(True)

        # Plot learning rate
        plt.subplot(2, 1, 2)
        for stage in df["stage"].unique():
            stage_data = df[(df["stage"] == stage) & (df["type"] == "train") & ("learning_rate" in df.columns)]
            if not stage_data.empty:
                plt.plot(stage_data["step"], stage_data["learning_rate"], label=f"{stage}")

        plt.xlabel("Training Steps")
        plt.ylabel("Learning Rate")
        plt.title("Learning Rate Schedule")
        plt.legend()
        plt.grid(True)

        # Save the plot
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/training_metrics.png")
        plt.close()

        self.logger.info(f"Training metrics plot saved to {self.output_dir}/training_metrics.png")

class FineTuningPipeline:
    """Main interface for running the entire fine-tuning pipeline"""

    def __init__(
        self,
        model_type="deepseek",
        output_dir="./output",
        data_dir="./data",
        use_wandb=True,
        wandb_project="genetic-llm-fine-tuning",
        wandb_entity=None
    ):
        self.model_type = model_type.lower()
        self.output_dir = output_dir
        self.data_dir = data_dir
        self.use_wandb = use_wandb

        # Setup wandb
        if use_wandb:
            import wandb
            wandb.login()
            self.wandb_project = wandb_project
            self.wandb_entity = wandb_entity

        # Create model config based on model type
        if self.model_type == "deepseek":
            self.model_config = DeepseekConfig(output_dir=output_dir)
        elif self.model_type == "qwen":
            self.model_config = QwenConfig(output_dir=output_dir)
        elif self.model_type == "mistral":
            self.model_config = MistralConfig(output_dir=output_dir)
        else:
            raise ValueError(f"Unsupported model type: {model_type}. Choose from 'deepseek', 'qwen', or 'mistral'.")

        # Initialize training manager
        self.training_manager = TrainingManager(self.model_config)

        # Initialize evaluator
        self.evaluator = ModelEvaluator(self.training_manager.model_manager)

        # Setup paths
        self.domain_data_path = os.path.join(data_dir, "domain_adaptation_data.json")
        self.instruction_data_path = os.path.join(data_dir, "instruction_data.json")
        self.hallucination_data_path = os.path.join(data_dir, "hallucination_data.json")

    def prepare_datasets(self, domain_csv_path=None, instruction_json_path=None, hallucination_json_path=None):
        """Prepare datasets for training from raw data files"""

        # Process domain adaptation data
        if domain_csv_path and not os.path.exists(self.domain_data_path):
            print(f"Preparing domain adaptation data from {domain_csv_path}")
            df = pd.read_csv(domain_csv_path)

            domain_data = []
            for _, row in df.iterrows():
                # Create a chat message format
                domain_item = {
                    "messages": [
                        {"role": "user", "content": f"Gene: {row['gene_symbol']}\nInformation: {row['abstract']}"},
                        {"role": "assistant", "content": row['text']}
                    ]
                }
                domain_data.append(domain_item)

            # Save as JSON
            with open(self.domain_data_path, 'w') as f:
                json.dump(domain_data, f, indent=2)

            print(f"Domain adaptation data saved to {self.domain_data_path}")

        # Process instruction data
        if instruction_json_path and not os.path.exists(self.instruction_data_path):
            print(f"Preparing instruction data from {instruction_json_path}")

            with open(instruction_json_path, 'r') as f:
                raw_data = json.load(f)

            instruction_data = []
            for item in raw_data:
                # Create a chat message format
                instruction_item = {
                    "messages": [
                        {"role": "user", "content": item["instruction"]},
                        {"role": "assistant", "content": item["response"]}
                    ]
                }
                instruction_data.append(instruction_item)

            # Save as JSON
            with open(self.instruction_data_path, 'w') as f:
                json.dump(instruction_data, f, indent=2)

            print(f"Instruction data saved to {self.instruction_data_path}")

        # Process hallucination data
        if hallucination_json_path and not os.path.exists(self.hallucination_data_path):
            print(f"Preparing hallucination data from {hallucination_json_path}")

            with open(hallucination_json_path, 'r') as f:
                raw_data = json.load(f)

            hallucination_data = []
            for item in raw_data:
                # Create a chat message format
                hallucination_item = {
                    "messages": [
                        {"role": "user", "content": f"Review the following statement about the gene {item['gene']}: {item['hallucinated_statement']}"},
                        {"role": "assistant", "content": item["corrected_statement"]}
                    ]
                }
                hallucination_data.append(hallucination_item)

            # Save as JSON
            with open(self.hallucination_data_path, 'w') as f:
                json.dump(hallucination_data, f, indent=2)

            print(f"Hallucination data saved to {self.hallucination_data_path}")

    def run_pipeline(self, skip_stages=None):
        """Run the entire fine-tuning pipeline"""
        skip_stages = skip_stages or []

        # Check if datasets exist
        if not os.path.exists(self.domain_data_path):
            raise FileNotFoundError(f"Domain adaptation data not found at {self.domain_data_path}. Please run prepare_datasets first.")

        if not os.path.exists(self.instruction_data_path):
            raise FileNotFoundError(f"Instruction data not found at {self.instruction_data_path}. Please run prepare_datasets first.")

        if not os.path.exists(self.hallucination_data_path):
            raise FileNotFoundError(f"Hallucination data not found at {self.hallucination_data_path}. Please run prepare_datasets first.")

        # Stage 1: Domain Adaptation
        if "domain" not in skip_stages:
            print("Starting Stage 1: Domain Adaptation Fine-Tuning")
            domain_model_path = self.training_manager.domain_adaptation_training(self.domain_data_path)
        else:
            print("Skipping Stage 1: Domain Adaptation Fine-Tuning")
            domain_model_path = None

        # Stage 2: Instruction Tuning
        if "instruction" not in skip_stages:
            print("Starting Stage 2: Task-Specific Instruction Tuning")
            instruction_model_path = self.training_manager.instruction_tuning(
                self.instruction_data_path,
                base_model_path=domain_model_path
            )
        else:
            print("Skipping Stage 2: Task-Specific Instruction Tuning")
            instruction_model_path = domain_model_path

        # Stage 3: Factuality Enhancement with DPO
        if "dpo" not in skip_stages:
            print("Starting Stage 3: Factuality Enhancement with DPO")
            dpo_model_path = self.training_manager.factuality_enhancement_dpo(
                self.hallucination_data_path,
                base_model_path=instruction_model_path
            )
        else:
            print("Skipping Stage 3: Factuality Enhancement with DPO")
            dpo_model_path = instruction_model_path

        # Stage 5: Optimization for Deployment
        if "optimize" not in skip_stages:
            print("Starting Stage 5: Optimization for Deployment")
            final_model_path = self.training_manager.optimize_for_deployment(
                model_path=dpo_model_path or instruction_model_path or domain_model_path
            )
        else:
            print("Skipping Stage 5: Optimization for Deployment")
            final_model_path = dpo_model_path or instruction_model_path or domain_model_path

        # Stage 6: Evaluation and Plotting
        if "evaluate" not in skip_stages:
            print("Starting Stage 6: Evaluation and Plotting")

            # Evaluate on instruction dataset
            instruction_accuracy, _ = self.evaluator.evaluate_model(
                final_model_path,
                self.instruction_data_path,
                stage_name="instruction"
            )

            # Evaluate on hallucination dataset
            hallucination_accuracy, _ = self.evaluator.evaluate_model(
                final_model_path,
                self.hallucination_data_path,
                stage_name="hallucination"
            )

            # Plot training metrics
            metrics_paths = []
            stage_names = []

            if domain_model_path:
                metrics_paths.append(domain_model_path)
                stage_names.append("Domain Adaptation")

            if instruction_model_path and instruction_model_path != domain_model_path:
                metrics_paths.append(instruction_model_path)
                stage_names.append("Instruction Tuning")

            if dpo_model_path and dpo_model_path != instruction_model_path:
                metrics_paths.append(dpo_model_path)
                stage_names.append("DPO")

            self.evaluator.plot_training_metrics(metrics_paths, stage_names)

            print(f"Instruction Evaluation Accuracy: {instruction_accuracy:.4f}")
            print(f"Hallucination Evaluation Accuracy: {hallucination_accuracy:.4f}")
        else:
            print("Skipping Stage 6: Evaluation and Plotting")

        print(f"Fine-tuning pipeline completed. Final model saved at: {final_model_path}")
        return final_model_path

def convert_phase1_data_to_fine_tuning_format(output_dir="./data"):
    """
    Convert data from Phase 1 into the correct format for fine-tuning

    Args:
        output_dir (str): Directory to save the prepared data

    Returns:
        dict: Paths to the prepared data files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Paths
    domain_csv_path = "/content/data/domain_adaptation_data_enhanced.csv"
    instruction_json_path = "./data/instruction_response_pairs.json"
    hallucination_json_path = "./data/hallucination_pairs.json"

    # Output paths
    domain_output_path = os.path.join(output_dir, "domain_adaptation_data.json")
    instruction_output_path = os.path.join(output_dir, "instruction_data.json")
    hallucination_output_path = os.path.join(output_dir, "hallucination_data.json")

    # 1. Process domain adaptation data
    if os.path.exists(domain_csv_path):
        df = pd.read_csv(domain_csv_path)

        domain_data = []
        for _, row in df.iterrows():
            # Create a chat message format
            domain_item = {
                "messages": [
                    {"role": "user", "content": f"Gene: {row['gene_symbol']}\nInformation: {row['abstract']}"},
                    {"role": "assistant", "content": row['text']}
                ]
            }
            domain_data.append(domain_item)

        # Save as JSON
        with open(domain_output_path, 'w') as f:
            json.dump(domain_data, f, indent=2)

        print(f"Domain adaptation data saved to {domain_output_path}")
    else:
        print(f"Warning: Domain adaptation CSV not found at {domain_csv_path}")

    # 2. Process instruction data
    if os.path.exists(instruction_json_path):
        with open(instruction_json_path, 'r') as f:
            raw_data = json.load(f)

        instruction_data = []
        for item in raw_data:
            # Create a chat message format
            instruction_item = {
                "messages": [
                    {"role": "user", "content": item["instruction"]},
                    {"role": "assistant", "content": item["response"]}
                ]
            }
            instruction_data.append(instruction_item)

        # Save as JSON
        with open(instruction_output_path, 'w') as f:
            json.dump(instruction_data, f, indent=2)

        print(f"Instruction data saved to {instruction_output_path}")
    else:
        print(f"Warning: Instruction JSON not found at {instruction_json_path}")

    # 3. Process hallucination data
    if os.path.exists(hallucination_json_path):
        with open(hallucination_json_path, 'r') as f:
            raw_data = json.load(f)

        hallucination_data = []
        # Get gene symbols from domain adaptation data if available
        if os.path.exists(domain_csv_path):
            domain_df = pd.read_csv(domain_csv_path)
            gene_symbols = domain_df['gene_symbol'].unique().tolist()
        else:
            gene_symbols = [] # Handle case where domain data is missing

        for i, item in enumerate(raw_data):
            # Create a chat message format, adding 'gene' from gene_symbols
            hallucination_item = {
                "messages": [
                    {"role": "user", "content": f"Review the following statement about the gene {gene_symbols[i % len(gene_symbols)] if gene_symbols else 'unknown gene'}: {item['hallucinated_statement']}"},  # Use gene_symbols or 'unknown gene' if empty
                    {"role": "assistant", "content": item["corrected_statement"]}
                ],
                "gene": gene_symbols[i % len(gene_symbols)] if gene_symbols else 'unknown gene' # Add 'gene' to the item
            }
            hallucination_data.append(hallucination_item)

        # Save as JSON
        with open(hallucination_output_path, 'w') as f:
            json.dump(hallucination_data, f, indent=2)

        print(f"Hallucination data saved to {hallucination_output_path}")
    else:
        print(f"Warning: Hallucination JSON not found at {hallucination_json_path}")

    return {
        "domain_data": domain_output_path,
        "instruction_data": instruction_output_path,
        "hallucination_data": hallucination_output_path
    }

def main():
    """Main function to run the fine-tuning pipeline"""
    # Convert Phase 1 data to fine-tuning format
    data_paths = convert_phase1_data_to_fine_tuning_format()

    # Choose model type
    model_types = ["deepseek"] # , "qwen", "mistral"

    for model_type in model_types:
        print(f"\n{'=' * 50}")
        print(f"Starting fine-tuning pipeline for {model_type}")
        print(f"{'=' * 50}")

        # Initialize pipeline
        pipeline = FineTuningPipeline(
            model_type=model_type,
            output_dir=f"./output/{model_type}",
            data_dir="./data"
        )

        # Run pipeline
        pipeline.run_pipeline(skip_stages=["domain", "instruction"])

    print("\nFine-tuning completed for all models!")


if __name__ == "__main__":
    main()

# !pip install transformers_stream_generator

# !pip install --upgrade transformers
# !pip install flash-attn --force-reinstall

# !pip install transformers==4.31.0
from transformers import TrainingArguments

# 1. Clear CUDA cache
import torch
torch.cuda.empty_cache()

# 2. Deallocate specific models/tensors
import gc
# Delete any large variables you don't need right now
# del model  # If you have a model variable
# del dataset, dataloader  # If you have these
gc.collect()  # Force garbage collection
torch.cuda.empty_cache()  # Clear cache again after garbage collection

# !pip install trl transformers accelerate peft datasets bitsandbytes

























  

    def get_training_args(self, stage_name: str = None):
        """Get training arguments for Trainer"""
        output_dir = self.output_dir
        if stage_name:
            output_dir = os.path.join(output_dir, stage_name)

        return TrainingArguments(
            output_dir=output_dir,
            num_train_epochs=self.num_train_epochs,
            per_device_train_batch_size=self.per_device_train_batch_size,
            per_device_eval_batch_size=self.per_device_eval_batch_size,
            gradient_accumulation_steps=self.gradient_accumulation_steps,
            gradient_checkpointing=self.gradient_checkpointing,
            learning_rate=self.learning_rate,
            weight_decay=self.weight_decay,
            warmup_ratio=self.warmup_ratio,
            optim=self.optim,
            lr_scheduler_type=self.lr_scheduler_type,
            save_total_limit=self.save_total_limit,
            save_strategy=self.save_strategy,
            save_steps=self.save_steps,
            logging_steps=self.logging_steps,
            evaluation_strategy=self.evaluation_strategy,
            eval_steps=self.eval_steps,
            fp16=self.fp16,
            bf16=self.bf16,
            seed=self.seed,
            report_to="wandb"
        )

    def get_lora_config(self):
        """Get LoRA configuration for PEFT"""
        return LoraConfig(
            r=self.lora_r,
            lora_alpha=self.lora_alpha,
            lora_dropout=self.lora_dropout,
            target_modules=self.target_modules,
            bias="none",
            task_type=TaskType.CAUSAL_LM
        )

    def estimate_memory_usage(self, model_size_gb=None):
        """Estimate memory usage based on configured parameters"""
        # Basic estimation - adjust formula as needed based on empirical testing
        base_memory = 0.5  # Base overhead in GB

        # Model memory (with quantization applied)
        if model_size_gb:
            model_memory = model_size_gb * (0.25 if self.use_4bit else 0.5 if self.load_in_8bit else 1.0)
        else:
            model_memory = 0  # Unknown

        # Batch memory
        seq_factor = self.max_seq_length / 1024  # Normalize to 1k tokens
        batch_memory = seq_factor * self.per_device_train_batch_size * 0.4  # Estimated GB per normalized batch

        # Additional overhead for gradients etc.
        if not self.gradient_checkpointing:
            batch_memory *= 1.5

        total_estimate = base_memory + model_memory + batch_memory

        return {
            "estimated_total_gb": total_estimate,
            "model_memory_gb": model_memory,
            "batch_memory_gb": batch_memory,
            "base_overhead_gb": base_memory
        }

class DeepseekConfig(ModelConfig):
    """Configuration for Deepseek model fine-tuning"""

    def __init__(self, **kwargs):
        default_model = "deepseek-ai/deepseek-llm-7b-base"
        default_target_modules = ["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj", "up_proj", "down_proj"]

        kwargs.setdefault("model_name", default_model)
        kwargs.setdefault("target_modules", default_target_modules)
        super().__init__(**kwargs)

        # Deepseek-specific settings
        self.output_dir = os.path.join(self.output_dir, "deepseek")


class QwenConfig(ModelConfig):
    """Configuration for Qwen model fine-tuning"""

    def __init__(self, **kwargs):
        default_model = "Qwen/Qwen-7B"
        default_target_modules = ["c_attn", "c_proj", "w1", "w2"]

        kwargs.setdefault("model_name", default_model)
        kwargs.setdefault("target_modules", default_target_modules)
        super().__init__(**kwargs)

        # Qwen-specific settings
        self.output_dir = os.path.join(self.output_dir, "qwen")


class MistralConfig(ModelConfig):
    """Configuration for Mistral model fine-tuning"""

    def __init__(self, **kwargs):
        default_model = "mistralai/Mistral-7B-v0.1"
        default_target_modules = ["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj", "up_proj", "down_proj"]

        kwargs.setdefault("model_name", default_model)
        kwargs.setdefault("target_modules", default_target_modules)
        super().__init__(**kwargs)

        # Mistral-specific settings
        self.output_dir = os.path.join(self.output_dir, "mistral")

class ModelManager:
    """Manager for loading and preparing models for fine-tuning"""

    def __init__(self, config: ModelConfig):
        self.config = config
        self.model = None
        self.tokenizer = None
        self.peft_model = None

    def load_model_and_tokenizer(self):
        """Load base model and tokenizer"""
        print(f"Loading model: {self.config.model_name}")

        # Load tokenizer
        self.tokenizer = AutoTokenizer.from_pretrained(
            self.config.tokenizer_name,
            use_fast=True,
            trust_remote_code=True
        )

        # Set padding token if not set
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
            self.tokenizer.pad_token_id = self.tokenizer.eos_token_id

        # Ensure special tokens are set properly
        # For chat models, we might need specific tokens
        self.ensure_chat_tokens()

        # Load model with quantization if specified
        if self.config.use_4bit:
            self.model = AutoModelForCausalLM.from_pretrained(
                self.config.model_name,
                quantization_config=BitsAndBytesConfig(
                    load_in_4bit=True,
                    bnb_4bit_compute_dtype=self.config.bnb_4bit_compute_dtype,
                    bnb_4bit_quant_type=self.config.bnb_4bit_quant_type,
                    bnb_4bit_use_double_quant=self.config.use_nested_quant
                ),
                torch_dtype=self.config.bnb_4bit_compute_dtype,
                trust_remote_code=True,
                use_flash_attention_2=self.config.use_flash_attention
            )
        elif self.config.load_in_8bit:
            self.model = AutoModelForCausalLM.from_pretrained(
                self.config.model_name,
                load_in_8bit=True,
                torch_dtype=torch.float16,
                trust_remote_code=True,
                use_flash_attention_2=self.config.use_flash_attention
            )
        else:
            self.model = AutoModelForCausalLM.from_pretrained(
                self.config.model_name,
                torch_dtype=torch.bfloat16,
                trust_remote_code=True,
                use_flash_attention_2=self.config.use_flash_attention
            )

        # Prepare model for kbit training if needed
        if self.config.use_4bit or self.config.load_in_8bit:
            self.model = prepare_model_for_kbit_training(
                self.model,
                use_gradient_checkpointing=self.config.gradient_checkpointing
            )

        # Enable gradient checkpointing for memory efficiency
        if self.config.gradient_checkpointing:
            self.model.gradient_checkpointing_enable()

        return self.model, self.tokenizer

    def prepare_for_peft(self):
        """Prepare model for PEFT/LoRA fine-tuning"""
        if self.model is None:
            self.load_model_and_tokenizer()

        print("Preparing model for PEFT/LoRA fine-tuning")
        lora_config = self.config.get_lora_config()
        self.peft_model = get_peft_model(self.model, lora_config)
        self.peft_model.print_trainable_parameters()

        return self.peft_model

    def ensure_chat_tokens(self):
        """Ensure chat-specific tokens are properly set"""
        # Different models might have different chat formats
        model_type = self.identify_model_type()

        if model_type == "deepseek":
            self.tokenizer.chat_template = "{% for message in messages %}\n{% if message.role == 'user' %}{{message.content}}\n{% elif message.role == 'assistant' %}{{message.content}}\n{% endif %}{% endfor %}"
        elif model_type == "qwen":
            # Handle Qwen specific tokens
            if not hasattr(self.tokenizer, 'chat_template'):
                self.tokenizer.chat_template = "{% for message in messages %}\n{% if message['role'] == 'user' %}\n{{ '<|im_start|>user\n' + message['content'] + '<|im_end|>\n' }}\n{% elif message['role'] == 'assistant' %}\n{{ '<|im_start|>assistant\n' + message['content'] + '<|im_end|>\n' }}\n{% elif message['role'] == 'system' %}\n{{ '<|im_start|>system\n' + message['content'] + '<|im_end|>\n' }}\n{% endif %}\n{% endfor %}\n{% if add_generation_prompt %}\n{{ '<|im_start|>assistant\n' }}\n{% endif %}"
        elif model_type == "mistral":
            # Handle Mistral specific tokens
            if not hasattr(self.tokenizer, 'chat_template'):
                self.tokenizer.chat_template = "{% for message in messages %}\n{% if message['role'] == 'user' %}\n{{ '<|user|>\n' + message['content'] + '<|assistant|>\n' }}\n{% elif message['role'] == 'assistant' %}\n{{ message['content'] + '<|endoftext|>\n' }}\n{% elif message['role'] == 'system' %}\n{{ '<|system|>\n' + message['content'] + '<|user|>\n' }}\n{% endif %}\n{% endfor %}"

    def identify_model_type(self):
        """Identify the type of model based on model name"""
        model_name = self.config.model_name.lower()
        if "deepseek" in model_name:
            return "deepseek"
        elif "qwen" in model_name:
            return "qwen"
        elif "mistral" in model_name:
            return "mistral"
        else:
            return "unknown"

from transformers import BitsAndBytesConfig

class DatasetPreparation:
    """Utilities for preparing datasets for different training stages"""

    def __init__(self, tokenizer, max_seq_length=2048):
        self.tokenizer = tokenizer
        self.max_seq_length = max_seq_length

    def format_chat_message(self, example):
        """Format messages into chat template compatible with tokenizer"""
        messages = example["messages"]
        return {"text": self.tokenizer.apply_chat_template(
            messages,
            tokenize=False,
            add_generation_prompt=False
        )}

    def preprocess_function(self, examples):
        """Tokenize text examples"""
        batch_size = len(examples["text"])
        inputs = [f"{text}" for text in examples["text"]]
        model_inputs = self.tokenizer(
            inputs,
            max_length=self.max_seq_length,
            padding="max_length",
            truncation=True,
            return_tensors="pt"
        )

        # Create labels (needed for causal language modeling)
        labels = model_inputs["input_ids"].clone()

        # Mask padding tokens
        for i in range(batch_size):
            labels[i][model_inputs["attention_mask"][i] == 0] = -100

        model_inputs["labels"] = labels
        return model_inputs

    def load_domain_adaptation_dataset(self, file_path, val_split=0.1):
        """Load domain adaptation dataset"""
        print(f"Loading domain adaptation dataset from {file_path}")

        # Read the dataset
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Create formatted dataset
        formatted_data = {
            "text": [],
        }

        for item in data:
            # Format according to domain adaptation format
            formatted_item = self.format_chat_message(item)
            formatted_data["text"].append(formatted_item["text"])

        # Split into train and validation
        total_samples = len(formatted_data["text"])
        val_size = int(total_samples * val_split)

        train_data = {
            "text": formatted_data["text"][val_size:]
        }

        val_data = {
            "text": formatted_data["text"][:val_size]
        }

        # Convert to dataset objects
        train_dataset = Dataset.from_dict(train_data)
        val_dataset = Dataset.from_dict(val_data)

        # Tokenize
        train_dataset = train_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        val_dataset = val_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        return train_dataset, val_dataset

    def load_instruction_dataset(self, file_path, val_split=0.1):
        """Load instruction dataset"""
        print(f"Loading instruction dataset from {file_path}")

        # Read the dataset
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Create formatted dataset
        formatted_data = {
            "text": [],
        }

        for item in data:
            # Format according to instruction format
            formatted_item = self.format_chat_message(item)
            formatted_data["text"].append(formatted_item["text"])

        # Split into train and validation
        total_samples = len(formatted_data["text"])
        val_size = int(total_samples * val_split)

        train_data = {
            "text": formatted_data["text"][val_size:]
        }

        val_data = {
            "text": formatted_data["text"][:val_size]
        }

        # Convert to dataset objects
        train_dataset = Dataset.from_dict(train_data)
        val_dataset = Dataset.from_dict(val_data)

        # Tokenize
        train_dataset = train_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        val_dataset = val_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        return train_dataset, val_dataset

    def load_hallucination_dataset(self, file_path, val_split=0.1):
        """Load hallucination dataset"""
        print(f"Loading hallucination dataset from {file_path}")

        # Read the dataset
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Create formatted dataset
        formatted_data = {
            "text": [],
        }

        for item in data:
            # Format according to hallucination format
            formatted_item = self.format_chat_message(item)
            formatted_data["text"].append(formatted_item["text"])

        # Split into train and validation
        total_samples = len(formatted_data["text"])
        val_size = int(total_samples * val_split)

        train_data = {
            "text": formatted_data["text"][val_size:]
        }

        val_data = {
            "text": formatted_data["text"][:val_size]
        }

        # Convert to dataset objects
        train_dataset = Dataset.from_dict(train_data)
        val_dataset = Dataset.from_dict(val_data)

        # Tokenize
        train_dataset = train_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        val_dataset = val_dataset.map(
            self.preprocess_function,
            batched=True,
            num_proc=4,
            remove_columns=["text"]
        )

        return train_dataset, val_dataset

from datasets import Dataset

class TrainingManager:
    """Manager for handling different stages of model training"""

    def __init__(self, model_config, output_dir="./output"):
        self.model_config = model_config
        self.output_dir = os.path.join(output_dir, model_config.__class__.__name__.replace("Config", "").lower())
        os.makedirs(self.output_dir, exist_ok=True)

        self.model_manager = ModelManager(model_config)
        self.model, self.tokenizer = self.model_manager.load_model_and_tokenizer()
        self.peft_model = None

        self.dataset_prep = DatasetPreparation(
            self.tokenizer,
            max_seq_length=model_config.max_seq_length
        )

        # Initialize logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(f"{self.output_dir}/training.log"),
                logging.StreamHandler()
            ]
        )

        self.logger = logging.getLogger(__name__)

    def initialize_wandb(self, stage_name):
        """Initialize Weights & Biases logging"""
        model_name = self.model_config.model_name.split("/")[-1]
        wandb.init(
            project="genetic-llm-fine-tuning",
            name=f"{model_name}-{stage_name}",
            config={
                "model_name": self.model_config.model_name,
                "stage": stage_name,
                "lora_r": self.model_config.lora_r,
                "lora_alpha": self.model_config.lora_alpha,
                "lora_dropout": self.model_config.lora_dropout,
                "max_seq_length": self.model_config.max_seq_length,
                "learning_rate": self.model_config.learning_rate,
                "num_train_epochs": self.model_config.num_train_epochs,
                "per_device_train_batch_size": self.model_config.per_device_train_batch_size,
                "gradient_accumulation_steps": self.model_config.gradient_accumulation_steps
            }
        )

    def domain_adaptation_training(self, domain_data_path):
        """Stage 1: Domain Adaptation Fine-Tuning"""
        stage_name = "domain_adaptation"
        self.logger.info(f"Starting {stage_name} stage")

        # Initialize W&B
        self.initialize_wandb(stage_name)

        # Prepare model for PEFT if not already done
        if self.peft_model is None:
            self.peft_model = self.model_manager.prepare_for_peft()

        # Load datasets
        train_dataset, val_dataset = self.dataset_prep.load_domain_adaptation_dataset(
            domain_data_path
        )

        # Get training arguments
        training_args = self.model_config.get_training_args(stage_name)

        # Initialize data collator
        data_collator = DataCollatorForLanguageModeling(
            tokenizer=self.tokenizer,
            mlm=False
        )

        # Initialize trainer
        trainer = Trainer(
            model=self.peft_model,
            args=training_args,
            train_dataset=train_dataset,
            eval_dataset=val_dataset,
            data_collator=data_collator,
            tokenizer=self.tokenizer
        )

        # Start training
        self.logger.info("Starting training...")
        train_result = trainer.train()

        # Save model
        self.logger.info("Saving model...")
        trainer.save_model(f"{self.output_dir}/{stage_name}")

        # Save training metrics
        metrics = train_result.metrics
        trainer.log_metrics("train", metrics)
        trainer.save_metrics("train", metrics)

        # Evaluate
        self.logger.info("Evaluating model...")
        eval_results = trainer.evaluate()
        trainer.log_metrics("eval", eval_results)
        trainer.save_metrics("eval", eval_results)

        # Close W&B
        wandb.finish()

        return f"{self.output_dir}/{stage_name}"

    def instruction_tuning(self, instruction_data_path, base_model_path=None):
        """Stage 2: Task-Specific Instruction Tuning"""
        stage_name = "instruction_tuning"
        self.logger.info(f"Starting {stage_name} stage")

        # Initialize W&B
        self.initialize_wandb(stage_name)

        # Load base model from previous stage if provided
        if base_model_path:
            self.logger.info(f"Loading base model from {base_model_path}")
            self.peft_model = PeftModel.from_pretrained(
                self.model,
                base_model_path,
                is_trainable=True
            )
        else:
            # Prepare model for PEFT if not already done
            if self.peft_model is None:
                self.peft_model = self.model_manager.prepare_for_peft()

        # Load datasets
        train_dataset, val_dataset = self.dataset_prep.load_instruction_dataset(
            instruction_data_path
        )

        # Get training arguments
        training_args = self.model_config.get_training_args(stage_name)

        # Initialize data collator
        data_collator = DataCollatorForLanguageModeling(
            tokenizer=self.tokenizer,
            mlm=False
        )

        # Initialize trainer
        trainer = Trainer(
            model=self.peft_model,
            args=training_args,
            train_dataset=train_dataset,
            eval_dataset=val_dataset,
            data_collator=data_collator,
            tokenizer=self.tokenizer
        )

        # Start training
 

       self.logger.info("Starting training...")
        train_result = trainer.train()

        # Save model
        self.logger.info("Saving model...")
        trainer.save_model(f"{self.output_dir}/{stage_name}")

        # Save training metrics
        metrics = train_result.metrics
        trainer.log_metrics("train", metrics)
        trainer.save_metrics("train", metrics)

        # Evaluate
        self.logger.info("Evaluating model...")
        eval_results = trainer.evaluate()
        trainer.log_metrics("eval", eval_results)
        trainer.save_metrics("eval", eval_results)

        # Close W&B
        wandb.finish()

        return f"{self.output_dir}/{stage_name}"

    def factuality_enhancement_dpo(self, hallucination_data_path, base_model_path=None):
        """Stage 3: Factuality Enhancement with DPO"""
        stage_name = "factuality_dpo"
        self.logger.info(f"Starting {stage_name} stage")

        try:
            from trl import DPOTrainer
        except ImportError:
            print("Installing TRL for DPO training...")
            !pip install -q trl
            from trl import DPOTrainer

        # Initialize W&B
        self.initialize_wandb(stage_name)

        # Load base model from previous stage if provided
        if base_model_path:
            self.logger.info(f"Loading base model from {base_model_path}")
            self.peft_model = PeftModel.from_pretrained(
                self.model,
                base_model_path,
                is_trainable=True
            )
        else:
            # Prepare model for PEFT if not already done
            if self.peft_model is None:
                self.peft_model = self.model_manager.prepare_for_peft()

        # Load and prepare the hallucination dataset for DPO
        with open(hallucination_data_path, 'r') as f:
            hallucination_data = json.load(f)

        # Process the data into DPO format
        dpo_data = {
            "prompt": [],
            "chosen": [],
            "rejected": []
        }

        for item in hallucination_data:
            # For each hallucination item, we need:
            # - prompt: the user query
            # - chosen: the corrected response (factual)
            # - rejected: the hallucinated response (non-factual)
            user_message = next(msg for msg in item["messages"] if msg["role"] == "user")
            assistant_message = next(msg for msg in item["messages"] if msg["role"] == "assistant")

            hallucinated_statement = user_message["content"]

            corrected_statement = assistant_message["content"]

            # Assuming each item has hallucinated_statement, corrected_statement, gene and evidence
            prompt = f"Please review the following statement about the gene {item['gene']}: {hallucinated_statement}"

            # Create chosen and rejected responses
            chosen = f"The statement contains inaccurate information. {corrected_statement}"
            rejected = f"The statement appears to be accurate. {hallucinated_statement} is scientifically valid."

            dpo_data["prompt"].append(prompt)
            dpo_data["chosen"].append(chosen)
            dpo_data["rejected"].append(rejected)

        # Convert to dataset
        dpo_dataset = Dataset.from_dict(dpo_data)

        # Split into train and eval
        dpo_dataset = dpo_dataset.train_test_split(test_size=0.1)

        # Configure DPO training arguments
        dpo_training_args = TrainingArguments(
            output_dir=f"{self.output_dir}/{stage_name}",
            num_train_epochs=self.model_config.num_train_epochs,
            per_device_train_batch_size=self.model_config.per_device_train_batch_size//2, # Reduced batch size due to paired samples
            gradient_accumulation_steps=self.model_config.gradient_accumulation_steps*2, # Increased to compensate for smaller batch
            gradient_checkpointing=self.model_config.gradient_checkpointing,
            learning_rate=self.model_config.learning_rate/2, # Lower learning rate for DPO
            weight_decay=self.model_config.weight_decay,
            lr_scheduler_type=self.model_config.lr_scheduler_type,
            save_strategy=self.model_config.save_strategy,
            save_steps=self.model_config.save_steps,
            logging_steps=self.model_config.logging_steps,
            evaluation_strategy=self.model_config.evaluation_strategy,
            eval_steps=self.model_config.eval_steps,
            fp16=self.model_config.fp16,
            bf16=self.model_config.bf16,
            report_to="wandb"
        )

        # Initialize DPO trainer
        dpo_trainer = DPOTrainer(
            model=self.peft_model,
            args=dpo_training_args,
            train_dataset=dpo_dataset["train"],
            eval_dataset=dpo_dataset["test"],
            # tokenizer=self.tokenizer,
            # beta=0.1,  # DPO hyperparameter
            # max_length=self.model_config.max_seq_length,
            # max_prompt_length=self.model_config.max_seq_length // 2
        )

        # Start DPO training
        self.logger.info("Starting DPO training...")
        dpo_trainer.train()

        # Save model
        self.logger.info("Saving model...")
        dpo_trainer.save_model(f"{self.output_dir}/{stage_name}")

        # Close W&B
        wandb.finish()

        return f"{self.output_dir}/{stage_name}"

    def optimize_for_deployment(self, model_path):
        """Stage 5: Optimize models for deployment"""
        stage_name = "optimized"
        self.logger.info(f"Starting {stage_name} stage")

        # Load the fine-tuned PEFT model
        self.logger.info(f"Loading fine-tuned model from {model_path}")
        fine_tuned_model = PeftModel.from_pretrained(
            self.model,
            model_path
        )

        # Merge adapter weights with base model for efficiency
        self.logger.info("Merging adapter weights with base model...")
        merged_model = fine_tuned_model.merge_and_unload()

        # Save the merged model
        merged_model_path = f"{self.output_dir}/{stage_name}"
        self.logger.info(f"Saving merged model to {merged_model_path}")
        merged_model.save_pretrained(merged_model_path)
        self.tokenizer.save_pretrained(merged_model_path)

        # Optional: Further quantize the model for deployment
        self.logger.info("Quantizing model for deployment...")
        try:
            # You might want to use a different quantization approach
            # depending on your deployment environment
            from optimum.onnxruntime import ORTModelForCausalLM

            # Convert to ONNX format
            ort_model = ORTModelForCausalLM.from_pretrained(
                merged_model_path,
                export=True
            )

            # Save optimized model
            optimized_path = f"{self.output_dir}/{stage_name}_onnx"
            ort_model.save_pretrained(optimized_path)
            self.logger.info(f"Saved ONNX model to {optimized_path}")

            return optimized_path
        except:
            self.logger.warning("ONNX conversion failed, returning merged model instead.")
            return merged_model_path


class ModelEvaluator:
    """Utility for evaluating model performance across different stages"""

    def __init__(self, model_manager, output_dir="./evaluation"):
        self.model_manager = model_manager
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        self.logger = logging.getLogger(__name__)

    def evaluate_model(self, model_path, test_data_path, stage_name=None):
        """Evaluate model on test dataset"""
        self.logger.info(f"Evaluating model from {model_path} on {test_data_path}")

        # Load the model
        model, tokenizer = self.model_manager.model, self.model_manager.tokenizer
        peft_model = PeftModel.from_pretrained(
            model,
            model_path
        )

        # Load the test data
        with open(test_data_path, 'r') as f:
            test_data = json.load(f)

        # Generate responses and evaluate
        results = {
            "input": [],
            "reference": [],
            "prediction": [],
            "is_correct": []
        }

        for item in tqdm(test_data, desc="Evaluating"):
            # Extract input and reference
            input_text = item["messages"][0]["content"]
            reference = item["messages"][1]["content"]

            # Generate prediction
            inputs = tokenizer(
                input_text,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=self.model_manager.config.max_seq_length
            ).to(peft_model.device)

            # Generate output
            with torch.no_grad():
                outputs = peft_model.generate(
                    input_ids=inputs["input_ids"],
                    attention_mask=inputs["attention_mask"],
                    max_new_tokens=512,
                    do_sample=False,
                    num_beams=1
                )

            # Decode prediction
            prediction = tokenizer.decode(outputs[0], skip_special_tokens=True)

            # Simple exact match evaluation (this should be improved with proper metrics)
            is_correct = prediction.strip() == reference.strip()

            # Store results
            results["input"].append(input_text)
            results["reference"].append(reference)
            results["prediction"].append(prediction)
            results["is_correct"].append(is_correct)

        # Calculate overall accuracy
        accuracy = sum(results["is_correct"]) / len(results["is_correct"])

        # Save results
        output_file = f"{self.output_dir}/results_{stage_name or 'evaluation'}.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)

        self.logger.info(f"Evaluation complete. Accuracy: {accuracy:.4f}")
        self.logger.info(f"Results saved to {output_file}")

        return accuracy, results

    def plot_training_metrics(self, metrics_paths, stage_names=None):
        """Plot training metrics from different stages"""
        metrics_data = []

        # Load metrics from each stage
        for i, path in enumerate(metrics_paths):
            stage = stage_names[i] if stage_names and i < len(stage_names) else f"Stage {i+1}"

            try:
                # Load training metrics
                with open(f"{path}/trainer_state.json", 'r') as f:
                    train_data = json.load(f)

                # Extract metrics
                log_history = train_data.get("log_history", [])

                for entry in log_history:
                    if "loss" in entry:
                        metrics_data.append({
                            "stage": stage,
                            "step": entry.get("step", 0),
                            "loss": entry.get("loss", 0),
                            "learning_rate": entry.get("learning_rate", 0),
                            "epoch": entry.get("epoch", 0),
                            "type": "train"
                        })
                    elif "eval_loss" in entry:
                        metrics_data.append({
                            "stage": stage,
                            "step": entry.get("step", 0),
                            "loss": entry.get("eval_loss", 0),
                            "epoch": entry.get("epoch", 0),
                            "type": "eval"
                        })
            except Exception as e:
                self.logger.warning(f"Error loading metrics from {path}: {e}")

        # Convert to DataFrame
        df = pd.DataFrame(metrics_data)

        if df.empty:
            self.logger.warning("No metrics data to plot.")
            return

        # Create plots
        plt.figure(figsize=(12, 10))

        # Plot training loss
        plt.subplot(2, 1, 1)
        for stage in df["stage"].unique():
            stage_data = df[(df["stage"] == stage) & (df["type"] == "train")]
            if not stage_data.empty:
                plt.plot(stage_data["step"], stage_data["loss"], label=f"{stage} Training")

            stage_eval_data = df[(df["stage"] == stage) & (df["type"] == "eval")]
            if not stage_eval_data.empty:
                plt.plot(stage_eval_data["step"], stage_eval_data["loss"], label=f"{stage} Evaluation", linestyle="--")

        plt.xlabel("Training Steps")
        plt.ylabel("Loss")
        plt.title("Training and Evaluation Loss")
        plt.legend()
        plt.grid(True)

        # Plot learning rate
        plt.subplot(2, 1, 2)
        for stage in df["stage"].unique():
            stage_data = df[(df["stage"] == stage) & (df["type"] == "train") & ("learning_rate" in df.columns)]
            if not stage_data.empty:
                plt.plot(stage_data["step"], stage_data["learning_rate"], label=f"{stage}")

        plt.xlabel("Training Steps")
        plt.ylabel("Learning Rate")
        plt.title("Learning Rate Schedule")
        plt.legend()
        plt.grid(True)

        # Save the plot
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/training_metrics.png")
        plt.close()

        self.logger.info(f"Training metrics plot saved to {self.output_dir}/training_metrics.png")

class FineTuningPipeline:
    """Main interface for running the entire fine-tuning pipeline"""

    def __init__(
        self,
        model_type="deepseek",
        output_dir="./output",
        data_dir="./data",
        use_wandb=True,
        wandb_project="genetic-llm-fine-tuning",
        wandb_entity=None
    ):
        self.model_type = model_type.lower()
        self.output_dir = output_dir
        self.data_dir = data_dir
        self.use_wandb = use_wandb

        # Setup wandb
        if use_wandb:
            import wandb
            wandb.login()
            self.wandb_project = wandb_project
            self.wandb_entity = wandb_entity

        # Create model config based on model type
        if self.model_type == "deepseek":
            self.model_config = DeepseekConfig(output_dir=output_dir)
        elif self.model_type == "qwen":
            self.model_config = QwenConfig(output_dir=output_dir)
        elif self.model_type == "mistral":
            self.model_config = MistralConfig(output_dir=output_dir)
        else:
            raise ValueError(f"Unsupported model type: {model_type}. Choose from 'deepseek', 'qwen', or 'mistral'.")

        # Initialize training manager
        self.training_manager = TrainingManager(self.model_config)

        # Initialize evaluator
        self.evaluator = ModelEvaluator(self.training_manager.model_manager)

        # Setup paths
        self.domain_data_path = os.path.join(data_dir, "domain_adaptation_data.json")
        self.instruction_data_path = os.path.join(data_dir, "instruction_data.json")
        self.hallucination_data_path = os.path.join(data_dir, "hallucination_data.json")

    def prepare_datasets(self, domain_csv_path=None, instruction_json_path=None, hallucination_json_path=None):
        """Prepare datasets for training from raw data files"""

        # Process domain adaptation data
        if domain_csv_path and not os.path.exists(self.domain_data_path):
            print(f"Preparing domain adaptation data from {domain_csv_path}")
            df = pd.read_csv(domain_csv_path)

            domain_data = []
            for _, row in df.iterrows():
                # Create a chat message format
                domain_item = {
                    "messages": [
                        {"role": "user", "content": f"Gene: {row['gene_symbol']}\nInformation: {row['abstract']}"},
                        {"role": "assistant", "content": row['text']}
                    ]
                }
                domain_data.append(domain_item)

            # Save as JSON
            with open(self.domain_data_path, 'w') as f:
                json.dump(domain_data, f, indent=2)

            print(f"Domain adaptation data saved to {self.domain_data_path}")

        # Process instruction data
        if instruction_json_path and not os.path.exists(self.instruction_data_path):
            print(f"Preparing instruction data from {instruction_json_path}")

            with open(instruction_json_path, 'r') as f:
                raw_data = json.load(f)

            instruction_data = []
            for item in raw_data:
                # Create a chat message format
                instruction_item = {
                    "messages": [
                        {"role": "user", "content": item["instruction"]},
                        {"role": "assistant", "content": item["response"]}
                    ]
                }
                instruction_data.append(instruction_item)

            # Save as JSON
            with open(self.instruction_data_path, 'w') as f:
                json.dump(instruction_data, f, indent=2)

            print(f"Instruction data saved to {self.instruction_data_path}")

        # Process hallucination data
        if hallucination_json_path and not os.path.exists(self.hallucination_data_path):
            print(f"Preparing hallucination data from {hallucination_json_path}")

            with open(hallucination_json_path, 'r') as f:
                raw_data = json.load(f)

            hallucination_data = []
            for item in raw_data:
                # Create a chat message format
                hallucination_item = {
                    "messages": [
                        {"role": "user", "content": f"Review the following statement about the gene {item['gene']}: {item['hallucinated_statement']}"},
                        {"role": "assistant", "content": item["corrected_statement"]}
                    ]
                }
                hallucination_data.append(hallucination_item)

            # Save as JSON
            with open(self.hallucination_data_path, 'w') as f:
                json.dump(hallucination_data, f, indent=2)

            print(f"Hallucination data saved to {self.hallucination_data_path}")

    def run_pipeline(self, skip_stages=None):
        """Run the entire fine-tuning pipeline"""
        skip_stages = skip_stages or []

        # Check if datasets exist
        if not os.path.exists(self.domain_data_path):
            raise FileNotFoundError(f"Domain adaptation data not found at {self.domain_data_path}. Please run prepare_datasets first.")

        if not os.path.exists(self.instruction_data_path):
            raise FileNotFoundError(f"Instruction data not found at {self.instruction_data_path}. Please run prepare_datasets first.")

        if not os.path.exists(self.hallucination_data_path):
            raise FileNotFoundError(f"Hallucination data not found at {self.hallucination_data_path}. Please run prepare_datasets first.")

        # Stage 1: Domain Adaptation
        if "domain" not in skip_stages:
            print("Starting Stage 1: Domain Adaptation Fine-Tuning")
            domain_model_path = self.training_manager.domain_adaptation_training(self.domain_data_path)
        else:
            print("Skipping Stage 1: Domain Adaptation Fine-Tuning")
            domain_model_path = None

        # Stage 2: Instruction Tuning
        if "instruction" not in skip_stages:
            print("Starting Stage 2: Task-Specific Instruction Tuning")
            instruction_model_path = self.training_manager.instruction_tuning(
                self.instruction_data_path,
                base_model_path=domain_model_path
            )
        else:
            print("Skipping Stage 2: Task-Specific Instruction Tuning")
            instruction_model_path = domain_model_path

        # Stage 3: Factuality Enhancement with DPO
        if "dpo" not in skip_stages:
            print("Starting Stage 3: Factuality Enhancement with DPO")
            dpo_model_path = self.training_manager.factuality_enhancement_dpo(
                self.hallucination_data_path,
                base_model_path=instruction_model_path
            )
        else:
            print("Skipping Stage 3: Factuality Enhancement with DPO")
            dpo_model_path = instruction_model_path

        # Stage 5: Optimization for Deployment
        if "optimize" not in skip_stages:
            print("Starting Stage 5: Optimization for Deployment")
            final_model_path = self.training_manager.optimize_for_deployment(
                model_path=dpo_model_path or instruction_model_path or domain_model_path
            )
        else:
            print("Skipping Stage 5: Optimization for Deployment")
            final_model_path = dpo_model_path or instruction_model_path or domain_model_path

        # Stage 6: Evaluation and Plotting
        if "evaluate" not in skip_stages:
            print("Starting Stage 6: Evaluation and Plotting")

            # Evaluate on instruction dataset
            instruction_accuracy, _ = self.evaluator.evaluate_model(
                final_model_path,
                self.instruction_data_path,
                stage_name="instruction"
            )

            # Evaluate on hallucination dataset
            hallucination_accuracy, _ = self.evaluator.evaluate_model(
                final_model_path,
                self.hallucination_data_path,
                stage_name="hallucination"
            )

            # Plot training metrics
            metrics_paths = []
            stage_names = []

            if domain_model_path:
                metrics_paths.append(domain_model_path)
                stage_names.append("Domain Adaptation")

            if instruction_model_path and instruction_model_path != domain_model_path:
                metrics_paths.append(instruction_model_path)
                stage_names.append("Instruction Tuning")

            if dpo_model_path and dpo_model_path != instruction_model_path:
                metrics_paths.append(dpo_model_path)
                stage_names.append("DPO")

            self.evaluator.plot_training_metrics(metrics_paths, stage_names)

            print(f"Instruction Evaluation Accuracy: {instruction_accuracy:.4f}")
            print(f"Hallucination Evaluation Accuracy: {hallucination_accuracy:.4f}")
        else:
            print("Skipping Stage 6: Evaluation and Plotting")

        print(f"Fine-tuning pipeline completed. Final model saved at: {final_model_path}")
        return final_model_path

def convert_phase1_data_to_fine_tuning_format(output_dir="./data"):
    """
    Convert data from Phase 1 into the correct format for fine-tuning

    Args:
        output_dir (str): Directory to save the prepared data

    Returns:
        dict: Paths to the prepared data files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Paths
    domain_csv_path = "/content/data/domain_adaptation_data_enhanced.csv"
    instruction_json_path = "./data/instruction_response_pairs.json"
    hallucination_json_path = "./data/hallucination_pairs.json"

    # Output paths
    domain_output_path = os.path.join(output_dir, "domain_adaptation_data.json")
    instruction_output_path = os.path.join(output_dir, "instruction_data.json")
    hallucination_output_path = os.path.join(output_dir, "hallucination_data.json")

    # 1. Process domain adaptation data
    if os.path.exists(domain_csv_path):
        df = pd.read_csv(domain_csv_path)

        domain_data = []
        for _, row in df.iterrows():
            # Create a chat message format
            domain_item = {
                "messages": [
                    {"role": "user", "content": f"Gene: {row['gene_symbol']}\nInformation: {row['abstract']}"},
                    {"role": "assistant", "content": row['text']}
                ]
            }
            domain_data.append(domain_item)

        # Save as JSON
        with open(domain_output_path, 'w') as f:
            json.dump(domain_data, f, indent=2)

        print(f"Domain adaptation data saved to {domain_output_path}")
    else:
        print(f"Warning: Domain adaptation CSV not found at {domain_csv_path}")

    # 2. Process instruction data
    if os.path.exists(instruction_json_path):
        with open(instruction_json_path, 'r') as f:
            raw_data = json.load(f)

        instruction_data = []
        for item in raw_data:
            # Create a chat message format
            instruction_item = {
                "messages": [
                    {"role": "user", "content": item["instruction"]},
                    {"role": "assistant", "content": item["response"]}
                ]
            }
            instruction_data.append(instruction_item)

        # Save as JSON
        with open(instruction_output_path, 'w') as f:
            json.dump(instruction_data, f, indent=2)

        print(f"Instruction data saved to {instruction_output_path}")
    else:
        print(f"Warning: Instruction JSON not found at {instruction_json_path}")

    # 3. Process hallucination data
    if os.path.exists(hallucination_json_path):
        with open(hallucination_json_path, 'r') as f:
            raw_data = json.load(f)

        hallucination_data = []
        # Get gene symbols from domain adaptation data if available
        if os.path.exists(domain_csv_path):
            domain_df = pd.read_csv(domain_csv_path)
            gene_symbols = domain_df['gene_symbol'].unique().tolist()
        else:
            gene_symbols = [] # Handle case where domain data is missing

        for i, item in enumerate(raw_data):
            # Create a chat message format, adding 'gene' from gene_symbols
            hallucination_item = {
                "messages": [
                    {"role": "user", "content": f"Review the following statement about the gene {gene_symbols[i % len(gene_symbols)] if gene_symbols else 'unknown gene'}: {item['hallucinated_statement']}"},  # Use gene_symbols or 'unknown gene' if empty
                    {"role": "assistant", "content": item["corrected_statement"]}
                ],
                "gene": gene_symbols[i % len(gene_symbols)] if gene_symbols else 'unknown gene' # Add 'gene' to the item
            }
            hallucination_data.append(hallucination_item)

        # Save as JSON
        with open(hallucination_output_path, 'w') as f:
            json.dump(hallucination_data, f, indent=2)

        print(f"Hallucination data saved to {hallucination_output_path}")
    else:
        print(f"Warning: Hallucination JSON not found at {hallucination_json_path}")

    return {
        "domain_data": domain_output_path,
        "instruction_data": instruction_output_path,
        "hallucination_data": hallucination_output_path
    }

def main():
    """Main function to run the fine-tuning pipeline"""
    # Convert Phase 1 data to fine-tuning format
    data_paths = convert_phase1_data_to_fine_tuning_format()

    # Choose model type
    model_types = ["deepseek"] # , "qwen", "mistral"

    for model_type in model_types:
        print(f"\n{'=' * 50}")
        print(f"Starting fine-tuning pipeline for {model_type}")
        print(f"{'=' * 50}")

        # Initialize pipeline
        pipeline = FineTuningPipeline(
            model_type=model_type,
            output_dir=f"./output/{model_type}",
            data_dir="./data"
        )

        # Run pipeline
        pipeline.run_pipeline(skip_stages=["domain", "instruction"])

    print("\nFine-tuning completed for all models!")


if __name__ == "__main__":
    main()

# !pip install transformers_stream_generator

# !pip install --upgrade transformers
# !pip install flash-attn --force-reinstall

# !pip install transformers==4.31.0
from transformers import TrainingArguments

# 1. Clear CUDA cache
import torch
torch.cuda.empty_cache()

# 2. Deallocate specific models/tensors
import gc
# Delete any large variables you don't need right now
# del model  # If you have a model variable
# del dataset, dataloader  # If you have these
gc.collect()  # Force garbage collection
torch.cuda.empty_cache()  # Clear cache again after garbage collection

# !pip install trl transformers accelerate peft datasets bitsandbytes

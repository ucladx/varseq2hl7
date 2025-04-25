# creates unique variant IDs that Beaker likes in the following pattern:
# 2a, 2b, ..., 2z, 2A, 2B, ..., 2Z, 2aa, 2ab, ..., 2aZ, 2ba, 2bb
alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
def _get_variant_id(idx):
    if idx < 0:
        return ""
    quotient, remainder = divmod(idx, len(alphabet))
    return _get_variant_id(quotient - 1) + alphabet[remainder]

def get_variant_id(idx):
    return f"2{_get_variant_id(idx)}"

LOINCS = {
    "48018-6": ["Gene Studied Name", "CWE"],
    "81252-9": ["Discrete Genetic Variant", "CWE"],
    "48005-3": ["Amino Acid Change p.HGVS", "ST"],
    "48004-6": ["DNA Change c.HGVS", "CWE"],
    "48006-1": ["Molecular Consequence", "CWE"],
    "62374-4": ["Human Reference Sequence Assembly Version", "CWE"],
    "53037-8": ["Genetic Sequence Variation Clinical Significance", "CWE"],
    "48000-4": ["Chromosome", "ST"],
    "81258-6": ["Allelic Frequency", "NM"],
    "51969-4": ["Genetic Analysis Report", "ST"],
    "82121-5": ["Allelic Read Depth", "NM"],
    "94076-7": ["Tumor Mutational Burden", "NM"],
    "47998-0": ["Variant Display Name", "ST"],
    "47999-8": ["DNA Region", "ST"],
    "83005-9": ["Variant Category", "ST"],
    "81695-9": ["Microsatellite Instability", "CWE"],
    "69548-6": ["Genetic Variant Assessment", "ST"],
    "93364-8": ["Genetic Variant Diagnostic Significance", "TX"],
    "48002-0": ["Genomic Source Class", "CWE"],
    "48019-4": ["DNA Change [Type]", "CWE"],
    "51958-7": ["Transcript Reference Sequence", "CWE"],
    "81290-9": ["Genomic DNA Change g.HGVS", "ST"],
    "81254-5": ["Genomic Allele Start-End", "NR"],
    "53034-5": ["Allelic State", "CWE"],
    "69547-8": ["Genomic Reference Allele", "ST"],
    "69551-0": ["Genomic Alternate Allele", "ST"],
    "81300-6": ["Structural Variant Length", "NM"],
    "7400052": ["Protein Reference Sequence", "ST"],
    "7102415": ["Tumor Type", "CWE"],
    "7102423": ["20x Depth Percentage", "ST"],
    "7102424": ["200x Depth Percentage", "ST"],
    "7102425": ["500x Depth Percentage", "ST"],
    "7102426": ["Mean Depth", "ST"],
}

# Maps Sequence Ontology Terms to EPIC molecular consequence and DNA change type
SEQ_ONTOLOGY_MAP = {
    "transcript_ablation": "Transcript Ablation",
    "splice_acceptor_variant": "Splice Acceptor Variant",
    "splice_donor_variant": "Splice Donor Variant",
    "stop_gained": "Nonsense",
    "frameshift_variant": "Frameshift Variant",
    "stop_lost": "Stop Lost",
    "start_lost": "Start Lost",
    "transcript_amplification": "Transcript Amplification",
    "feature_elongation": "Feature Elongation",
    "feature_truncation": "Feature Truncation",
    "inframe_insertion": "Inframe Insertion",
    "disruptive_inframe_insertion": "Inframe Insertion",
    "inframe_deletion": "Inframe Deletion",
    "disruptive_inframe_deletion": "Inframe Deletion",
    "missense_variant": "Missense Variant",
    "protein_altering_variant": "Protein Altering Variant",
    "splice_donor_5th_base_variant": "Splice Donor Variant",
    "splice_region_variant": "Splice Region Variant",
    "splice_donor_region_variant": "Splice Region Variant",
    "splice_polypyrimidine_tract_variant": "Splice Region Variant",
    "incomplete_terminal_codon_variant": "Incomplete Terminal Codon Variant",
    "start_retained_variant": "Start Retained Variant",
    "stop_retained_variant": "Stop Retained Variant",
    "synonymous_variant": "Synonymous Variant",
    "coding_sequence_variant": "Coding Sequence Variant",
    "initiator_codon_variant": "Coding Sequence Variant",
    "mature_miRNA_variant": "Mature miRNA Variant",
    "5_prime_UTR_variant": "5 Prime UTR Variant",
    "3_prime_UTR_variant": "3 Prime UTR Variant",
    "non_coding_transcript_exon_variant": "Non Coding Transcript Exon Variant",
    "intron_variant": "Intron Variant",
    "NMD_transcript_variant": "NMD Transcript Variant",
    "non_coding_transcript_variant": "Non Coding Transcript Variant",
    "coding_transcript_variant": "Coding Sequence Variant",
    "upstream_gene_variant": "Regulatory Region Variant", # special case for TERT promoter --- upstream/downstream variants worth reporting are almost always in regulatory regions
    "downstream_gene_variant": "Regulatory Region Variant", # see above
    "TFBS_ablation": "TFBS Ablation",
    "TFBS_amplification": "TFBS Amplification",
    "TF_binding_site_variant": "TF Binding Site Variant",
    "regulatory_region_ablation": "Regulatory Region Ablation",
    "regulatory_region_amplification": "Regulatory Region Amplification",
    "regulatory_region_variant": "Regulatory Region Variant",
    "intergenic_variant": "Intergenic Variant",
    "sequence_variant": "Coding Sequence Variant",
}

LAB_CODES = {
    "UCLA Heme v2": "LAB9121R^Hematologic Malignancy Sequencing Panel^BKREAP^^^^^^HEMATOLOGIC MALIGNANCY SEQUENCING PANEL",
    "UCLA Pan-Cancer All v1": "LAB9055^Pan-cancer Solid Tumor Panel^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL",
}

import hl7.client
from flask import Flask, request, jsonify
from datetime import date
from textwrap import wrap

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
    "83005-9": ["Variant Category", "ST"],
    "7102415": ["Tumor Type", "CWE"],
    "81695-9": ["Microsatellite Instability", "CWE"],
    "69548-6": ["Genetic Variant Assessment", "ST"],
    "93364-8": ["Genetic Variant Diagnostic Significance", "ST"],
    "48002-0": ["Genomic Source Class", "CWE"],
    "48019-4": ["DNA Change [Type]", "CWE"],
    "51958-7": ["Transcript Reference Sequence", "CWE"],
    "7400052": ["Protein Reference Sequence", "ST"],
    "81290-9": ["Genomic DNA Change g.HGVS", "ST"],
    "81254-5": ["Genomic Allele Start-End", "NR"],
    "53034-5": ["Allelic State", "CWE"],
    "69547-8": ["Genomic Reference Allele", "ST"],
    "69551-0": ["Genomic Alternate Allele", "ST"],
    "7102423": ["20x Depth Percentage", "ST"],
    "7102424": ["200x Depth Percentage", "ST"],
    "7102425": ["500x Depth Percentage", "ST"],
    "7102426": ["Mean Depth", "ST"],
}

# Sequence Ontology -> EPIC molecular consequence and DNA change type
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
    "inframe_deletion": "Inframe Deletion",
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
    "mature_miRNA_variant": "Mature miRNA Variant",
    "5_prime_UTR_variant": "5 Prime UTR Variant",
    "3_prime_UTR_variant": "3 Prime UTR Variant",
    "non_coding_transcript_exon_variant": "Non Coding Transcript Exon Variant",
    "intron_variant": "Intron Variant",
    "NMD_transcript_variant": "NMD Transcript Variant",
    "non_coding_transcript_variant": "Non Coding Transcript Variant",
    "coding_transcript_variant": "Coding Sequence Variant",
    "upstream_gene_variant": "Upstream Gene Variant",
    "downstream_gene_variant": "Downstream Gene Variant",
    "TFBS_ablation": "TFBS Ablation",
    "TFBS_amplification": "TFBS Amplification",
    "TF_binding_site_variant": "TF Binding Site Variant",
    "regulatory_region_ablation": "Regulatory Region Ablation",
    "regulatory_region_amplification": "Regulatory Region Amplification",
    "regulatory_region_variant": "Regulatory Region Variant",
    "intergenic_variant": "Intergenic Variant",
    "sequence_variant": "Coding Sequence Variant",
}

def get_loinc_info(code):
    loinc_info = LOINCS.get(code)
    if loinc_info:
        loinc_name = loinc_info[0]
        loinc_type = loinc_info[1]
        return f"{loinc_type}|{code}^{loinc_name}^LOINC"
    else:
        raise RuntimeError("Invalid LOINC code")

class VarSeqInfo():
    def __init__(self, varseq_json, pt_info=None):
        self.obx_idx = 0
        self.varseq_json = varseq_json
        self.sample_info = self.varseq_json["sampleState"]
        self.sample_id = self.sample_info["sampleName"]
        self.mrn = self.sample_info["mrn"]
        self.pt_ln, self.pt_fn = self.get_pt_name()
        self.bday = self.get_date("dob")
        self.sex = self.sample_info["sex"][0]
        self.order_num = self.get_order_num()
        self.prov_id = self.get_prov_id()
        self.prov_ln, self.prov_fn = self.get_prov_name()
        self.date_ordered = self.get_date("dateOrdered")
        self.date_received = self.get_date("dateReceived")
        self.date_collected = self.get_date("dateCollected")
        self.date_sent = date.today().strftime("%Y%m%d")
        self.norm_sample_id = self.get_custom_field("N_SID")
        self.norm_order_num = self.get_custom_field("N_OrderID")
        self.norm_date_ordered = self.format_msh_date(self.get_custom_field("N_DateOrdered").split(" ")[0])
        self.norm_date_received = self.format_msh_date(self.get_custom_field("N_DateReceived").split(" ")[0])

    def get_obx_idx(self):
        self.obx_idx += 1
        return self.obx_idx

    def reset_obx_idx(self):
        self.obx_idx = 0
        return self.obx_idx

    def get_biomarkers(self):
        biomarkers = filter(lambda x: x["type"] == "VARIANT", self.varseq_json["biomarkers"])
        return list(biomarkers)

    def get_germline_variants(self):
        germline_biomarkers = filter(lambda x: x["tags"][0]["initials"] == "GL", self.get_all_variants())
        return list(germline_biomarkers) + self.varseq_json["germlineVariants"]

    def get_all_variants(self):
        all_variants = self.get_biomarkers() + self.varseq_json["germlineVariants"] + self.varseq_json["uncertainVariants"]
        variant_vafs = [(variant, variant["vaf"]) for variant in all_variants]
        variant_vafs.sort(key=lambda x: x[1], reverse=True)
        return [variant_vaf[0] for variant_vaf in variant_vafs]
        # return self.get_biomarkers() + self.varseq_json["germlineVariants"] + self.varseq_json["uncertainVariants"]

    def get_sig(self, sig_name):
        for biomarker in self.varseq_json["biomarkers"]:
            if biomarker["type"] == "SIGNATURE" and biomarker["geneName"] == sig_name:
                return biomarker["quantitativeValue"]
        return ""

    def get_msi(self): 
        return self.get_sig("MSI")

    def get_tmb(self):
        return self.get_sig("TMB")

    def get_custom_field(self, field_name):
        fields = self.varseq_json["customFields"]
        return fields[field_name] if field_name in fields else ""

    def get_prov_id(self):
        return self.get_custom_field("ProviderID")

    def get_order_num(self):
        return self.get_custom_field("OrderID")

    def get_pt_name(self):
        return self.varseq_json["sampleState"]["patientName"].split(" (")[0].split(",")

    def get_prov_name(self):
        return [fn_or_ln.strip() for fn_or_ln in self.varseq_json["sampleState"]["orderingPhysician"].split(",")][:2]

    def format_msh_date(self, date):
        [m, d, y] = date.split('/')
        return y + m + d

    def get_date(self, date_type):
        return self.format_msh_date(self.varseq_json["sampleState"][date_type])

    def get_tumor_type(self):
        tt = self.varseq_json["tumorType"]
        tt_abbrev = self.varseq_json["tumorTypeAbbrev"]
        return f"{tt} ({tt_abbrev})"

    def get_clin_sig(self, variant):
        if variant in self.varseq_json["biomarkers"]:
            return "Tier 1: Strong significance"
        else:
            return "Tier 3: Unknown clinical significance"
        
    def get_variant_type(self, variant):
        variant_type = variant["tags"][0]["initials"]
        if variant_type == "SM":
            return "Somatic"
        elif variant_type == "GL":
            return "Germline"
        else:
            return "Unknown"

    def insert_newlines(self, s, every=110):
        return '\.br\\'.join(wrap(s, every))

    def get_interp(self, variant):
        summary = variant["biomarkerSummary"]
        if summary:
            interp = summary["interpretation"]
            # TODO if we end up seeing more HTML elements, actually do this properly with regex (re.sub)
            interp = interp.replace("<p>", "").replace("</p>", "").replace("&nbsp;", " ")
            return self.insert_newlines(interp)
        else:
            return ""

    def get_chrom(self, variant):
        mutation = variant["mutation"]
        if mutation:
            return mutation["chr"]
        else:
            return ""

    def get_naf(self, variant):
        naf = variant["naf"]
        if naf == "null":
            return "<LoD"
        rounded_naf = round(naf, 2)
        if rounded_naf == 0:
            return "<LoD"
        else:
            return rounded_naf

    def get_consequence(self, variant):
        cons =  SEQ_ONTOLOGY_MAP.get(variant["sequenceOntology"])
        if cons:
            return cons
        else:
            raise Exception(f"SequenceOntology term {variant['sequenceOntology']} not found")

    def is_subst(self, variant):
        ref, alt = map(lambda b: b.replace('-', ''), variant["refAlt"].split("/"))
        return len(alt) == len(ref)

    def get_dna_change(self, variant):
        return "Substitution" if self.is_subst(variant) else "Insertion/Deletion"

    def get_variant_obx(self, variant_id, loinc_code, value):
        return f"""OBX|{self.get_obx_idx()}|{get_loinc_info(loinc_code)}|{variant_id}|{value}"""

    def get_variant_obx_function(self, variant_id):
        def wrapper_func(loinc_code, value):
            return self.get_variant_obx(variant_id, loinc_code, value)
        return wrapper_func

    def get_variant_obxs(self, variant, idx):
        display_name = variant["geneName"]
        variant_id = get_variant_id(idx)
        get_variant_obx = self.get_variant_obx_function(variant_id)
        pDot = variant["pDotThreeLetter"] if variant["pDotThreeLetter"] else "p.?"
        consequence = self.get_consequence(variant)
        dna_change = self.get_dna_change(variant)
        ref, alt = map(lambda b: b.replace('-', ''), variant["refAlt"].split("/"))
        obxs =  [
            get_variant_obx("47998-0", display_name), # Variant Display Name
            get_variant_obx("83005-9", "Simple"), # EPIC Variant Category (Simple, Complex, Fusion, etc.))
            get_variant_obx("81252-9", f"v1^{variant['hgvsWithGene']}^ClinVar-V"), # Discrete Genetic Variant
            get_variant_obx("48018-6", f"^{display_name}^"), # Variant Name (This is also a discrete field in EPIC)
            get_variant_obx("48005-3", f"{pDot}"), # Amino Acid Change p.HGVS
            get_variant_obx("48004-6", f"^{variant['cDot']}"), # DNA Change c.HGVS
            get_variant_obx("62374-4", "^GRCh38"), # Human Reference Sequence Assembly Version
            get_variant_obx("53037-8", f"^{self.get_clin_sig(variant)}"), # Genetic Sequence Variation Clinical Significance
            get_variant_obx("48000-4", self.get_chrom(variant)), # Chromosome
            get_variant_obx("81258-6", round(variant["vaf"], 2)), # Allelic Frequency
            get_variant_obx("82121-5", variant["altReadCount"]), # Allelic Read Depth
            get_variant_obx("69548-6", "Detected"), # Genetic Variant Assessment
            get_variant_obx("93364-8", self.get_interp(variant)), # Genetic Variant Diagnostic Significance
            get_variant_obx("48002-0", f"^{self.get_variant_type(variant)}"), # Genomic Source Class
            get_variant_obx("51958-7", f"{variant['transcriptName']}^{variant['transcriptName']}^RefSeq-T"), # Transcript Reference Sequence
            get_variant_obx("7400052", variant["proteinId"]),
            get_variant_obx("81290-9", variant["gDot"]), # Genomic DNA Change g.HGVS
            get_variant_obx("81254-5", f"{variant['start']}^{variant['stop']}"), # Genomic Allele Start-End
            get_variant_obx("69547-8", ref), # Genomic Reference Allele
            get_variant_obx("69551-0", alt), # Genomic Alternate Allele
            get_variant_obx("48006-1", f"^{consequence}"),
            get_variant_obx("48019-4", f"^{dna_change}"),
        ]
        return "\r\n".join(obxs)

    # While we send these variants in the tumor's HL7 message, we send them again here for the normal sample so we can get the normal allele frequency
    def get_normal_variant_obxs(self, variant, idx):
        display_name = variant["geneName"]
        variant_id = get_variant_id(idx)
        get_variant_obx = self.get_variant_obx_function(variant_id)
        obxs =  [
            get_variant_obx("47998-0", display_name), # Variant Display Name
            get_variant_obx("83005-9", "Simple"), # EPIC Variant Category (Simple, Complex, Fusion, etc.))
            get_variant_obx("81252-9", f"v1^{variant['hgvsWithGene']}^ClinVar-V"), # Discrete Genetic Variant
            get_variant_obx("48018-6", f"^{display_name}^"), # Variant Name (This is also a discrete field in EPIC)
            get_variant_obx("62374-4", "^GRCh38"), # Human Reference Sequence Assembly Version
            get_variant_obx("81258-6", self.get_naf(variant)), # Allelic Frequency
            get_variant_obx("69548-6", "Detected"), # Genetic Variant Assessment
            get_variant_obx("48002-0", f"^{self.get_variant_type(variant)}"), # Genomic Source Class
            get_variant_obx("53034-5", variant["zygosity"]),
        ]
        return "\r\n".join(obxs)

    def get_hl7_header(self):
        self.reset_obx_idx()
        bases_20x, bases_200x, bases_500x, covg_mean = self.get_covg_metrics()
        return f"""MSH|^~\&|RRH||Beaker||{self.date_sent}||ORU^R01|1|P|2.3||||||\r
PID|1||{self.mrn}^^^MRN^MRN||{self.pt_ln}^{self.pt_fn}^||{self.bday}|{self.sex}\r
ORC|RE\r
OBR|1|{self.order_num}|{self.sample_id}^Beaker|LAB9055^Pan-cancer Solid Tumor Panel^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL|||{self.date_ordered}|||||||||{self.prov_id}^{self.prov_ln}^{self.prov_fn}^^^^^^EPIC^^^^PROVID||||||{self.date_received}|||F\r
{self.get_variant_obx("2a", "7102415", f"^{self.get_tumor_type()}")}\r
{self.get_variant_obx("2a", "81695-9", f"^{self.get_msi()}")}\r
{self.get_variant_obx("2a", "94076-7", f"{self.get_tmb()}")}\r
{self.get_variant_obx("2a", "7102423", bases_20x)}\r
{self.get_variant_obx("2a", "7102424", bases_200x)}\r
{self.get_variant_obx("2a", "7102425", bases_500x)}\r
{self.get_variant_obx("2a", "7102426", covg_mean)}\r"""

    def get_normal_hl7_header(self):
        self.reset_obx_idx()
        return f"""MSH|^~\&|RRH||Beaker||{self.date_sent}||ORU^R01|1|P|2.3||||||\r
PID|1||{self.mrn}^^^MRN^MRN||{self.pt_ln}^{self.pt_fn}^||{self.bday}|{self.sex}\r
ORC|RE\r
OBR|1|{self.norm_order_num}|{self.norm_sample_id}^Beaker|LAB9056^Pan-cancer Panel, Comparator^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL|||{self.norm_date_ordered}|||||||||{self.prov_id}^{self.prov_ln}^{self.prov_fn}^^^^^^EPIC^^^^PROVID||||||{self.norm_date_received}|||F\r"""

    def get_hl7_variants(self):
        variants = self.get_all_variants()
        return "\r\n".join([self.get_variant_obxs(variant, idx) for idx, variant in enumerate(variants)])
    
    def get_normal_hl7_variants(self):
        variants = self.get_germline_variants()
        return "\r\n".join([self.get_normal_variant_obxs(variant, idx) for idx, variant in enumerate(variants)])

def get_hl7_msg(vsinfo):
    return vsinfo.get_hl7_header() + vsinfo.get_hl7_variants()

def get_normal_hl7_msg(vsinfo):
    header = vsinfo.get_normal_hl7_header()
    normal_variants = vsinfo.get_normal_hl7_variants()
    if normal_variants:
        return header + normal_variants
    else:
        return ""

def create_hl7_msgs(vs_json):
    vsinfo = VarSeqInfo(vs_json)
    tumor_msg, normal_msg = (get_hl7_msg(vsinfo), get_normal_hl7_msg(vsinfo))
    return tumor_msg, normal_msg

def send_hl7_msgs(vs_json):
    tumor_msg, normal_msg = create_hl7_msgs(vs_json)
    hostname = "interface-test.mednet.ucla.edu"
    port = 7199

    with open("./hl7_msg_T.txt", "w") as f:
        f.write(tumor_msg)
    # sends messages and prints the corresponding acks from Beaker
    print(hl7.client.MLLPClient(hostname, port).send_message(tumor_msg))
    if normal_msg:
        with open("./hl7_msg_N.txt", "w") as f:
            f.write(normal_msg)
        print(hl7.client.MLLPClient(hostname, port).send_message(normal_msg))
    return vs_json

app = Flask(__name__)

@app.route('/receivejson', methods=['POST'])
def receive_json():
    data = request.get_json()
    result = send_hl7_msgs(data)
    return jsonify(result), 200

# launches Flask server on localhost:5000
if __name__ == '__main__':
    app.run(debug=True)

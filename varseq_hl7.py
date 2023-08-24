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
    "7102416": ["Sample Library Average Size", "CWE"],
    "7102417": ["Run Loading Concentration", "CWE"],
    "7102418": ["Dna Nanodrop Concentration", "CWE"],
    "7102419": ["Dna Qubit Concentration", "CWE"],
    "7102420": ["Pool Library Concentration", "CWE"],
    "7102421": ["Pool Library Average Size", "CWE"],
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
    def __init__(self, varseq_json):
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
        variant_vaf_dict = {variant: variant["vaf"] for variant in all_variants}
        return sorted(all_variants, key=lambda variant: variant_vaf_dict[variant], reverse=True)
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

    # TODO test if EPIC can do this automatically with HTML tags
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

    # TODO replace with map from sequenceOntology -> EPIC term
    def get_effect(self, variant):
        mutation = variant["mutation"]
        if mutation:
            effect = mutation["effect"].replace("In-frame", "inframe")
            return effect
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
        obxs =  [
            get_variant_obx("47998-0", display_name), # Variant Display Name
            get_variant_obx("83005-9", "Simple"), # EPIC Variant Category (Simple, Complex, Fusion, etc.))
            get_variant_obx("81252-9", f"v1^{variant['hgvsWithGene']}^ClinVar-V"), # Discrete Genetic Variant
            get_variant_obx("48018-6", f"^{display_name}^"), # Variant Name (This is also a discrete field in EPIC)
            get_variant_obx("48005-3", f"{pDot}"), # Amino Acid Change p.HGVS
            get_variant_obx("48004-6", f"^{variant['cDot']}"), # DNA Change c.HGVS
            get_variant_obx("48006-1", f"^{variant['effect']}"), # Molecular Consequence (Missense, Nonsense, etc.)
            get_variant_obx("62374-4", "^GRCh38"), # Human Reference Sequence Assembly Version
            get_variant_obx("53037-8", f"^{self.get_clin_sig(variant)}"), # Genetic Sequence Variation Clinical Significance
            get_variant_obx("48000-4", self.get_chrom(variant)), # Chromosome
            get_variant_obx("81258-6", round(variant["vaf"], 2)), # Allelic Frequency
            get_variant_obx("82121-5", variant["altReadCount"]), # Allelic Read Depth
            get_variant_obx("69548-6", "Detected"), # Genetic Variant Assessment
            get_variant_obx("93364-8", self.get_interp(variant)), # Genetic Variant Diagnostic Significance
            get_variant_obx("48002-0", f"^{self.get_variant_type(variant)}"), # Genomic Source Class
        ]
        return "\r\n".join(obxs) + "\r\n"

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
        ]
        return "\r\n".join(obxs) + "\r\n"

    def get_hl7_header(self):
        self.reset_obx_idx()
        return f"""MSH|^~\&|RRH||Beaker||{self.date_sent}||ORU^R01|1|P|2.3||||||\r
PID|1||{self.mrn}^^^MRN^MRN||{self.pt_ln}^{self.pt_fn}^||{self.bday}|{self.sex}\r
ORC|RE\r
OBR|1|{self.order_num}|{self.sample_id}^Beaker|LAB9055^Pan-cancer Solid Tumor Panel^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL|||{self.date_ordered}|||||||||{self.prov_id}^{self.prov_ln}^{self.prov_fn}^^^^^^EPIC^^^^PROVID||||||{self.date_received}|||F\r
{self.get_variant_obx("2a", "7102415", f"^{self.get_tumor_type()}")}\r
{self.get_variant_obx("2a", "81695-9", f"^{self.get_msi()}")}\r
{self.get_variant_obx("2a", "94076-7", f"{self.get_tmb()}")}\r"""

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
    return vsinfo.get_normal_hl7_header() + vsinfo.get_normal_hl7_variants()

def create_hl7_msgs(vs_json):
    vsinfo = VarSeqInfo(vs_json)
    tumor_msg, normal_msg = (get_hl7_msg(vsinfo), get_normal_hl7_msg(vsinfo))
    return tumor_msg, normal_msg

def send_hl7_msgs(vs_json):
    tumor_msg, normal_msg = create_hl7_msgs(vs_json)
    hostname = "interface-test.mednet.ucla.edu"
    port = 17199

    # with open("./hl7_msg_N.txt", "w") as f:
    #     f.write(normal_msg)
    # with open("./hl7_msg_T.txt", "w") as f:
    #     f.write(tumor_msg)
    # sends messages and prints the corresponding acks from Beaker
    print(hl7.client.MLLPClient(hostname, port).send_message(tumor_msg))
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

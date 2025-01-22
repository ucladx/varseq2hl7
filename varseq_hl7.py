import hl7.client
from flask import Flask, request, jsonify
from datetime import date
from textwrap import wrap
from mappings import LOINCS, SEQ_ONTOLOGY_MAP, get_variant_id
import argparse

parser = argparse.ArgumentParser(description='HL7 Server')
parser.add_argument('--hostname', type=str, help='Hostname of the server')
parser.add_argument('--port', type=int, help='Port number of the server')
args = parser.parse_args()

HOSTNAME = args.hostname
PORT = args.port

class VarSeqInfo():
    def __init__(self, varseq_json):
        self.obx_idx = 0
        self.varseq_json = varseq_json
        self.sample_state = varseq_json["sampleState"]
        self.sample_id = self.sample_state["sampleName"].rstrip('R')
        self.mrn = self.get_mrn()
        self.pt_ln, self.pt_fn = self.get_pt_name()
        self.bday = self.get_date("dob")
        self.sex = self.get_sex(self.sample_state["sex"])
        self.prov_ln, self.prov_fn = self.get_prov_name()
        self.order_num = self.get_custom_field("OrderID")
        self.prov_id = self.get_prov_id()
        self.norm_sample_id = self.get_custom_field("N_SID").rstrip('R')
        self.norm_order_num = self.get_custom_field("N_OrderID")
        self.date_ordered = self.get_date("DateOrdered")
        self.date_received = self.get_date("DateReceived")
        self.norm_date_ordered = self.get_date("N_DateOrdered")
        self.norm_date_received = self.get_date("N_DateReceived")
        self.date_sent = self.get_date_sent()
        self.variants = self.get_all_variants()

    def next_obx_idx(self):
        self.obx_idx += 1
        return self.obx_idx

    def reset_obx_idx(self):
        self.obx_idx = 0
        return self.obx_idx

    def get_biomarker_variants(self):
        biomarkers = filter(lambda x: x["type"] == "VARIANT", self.varseq_json["biomarkers"])
        return list(biomarkers)

    def get_all_variants(self):
        biomarkers = self.get_biomarker_variants()
        biomarkers.sort(key=lambda x: x["vaf"], reverse=True)
        all_variants = self.varseq_json["germlineVariants"] + self.varseq_json["uncertainVariants"]
        all_variants.sort(key=lambda x: x["vaf"], reverse=True) # sort variants by descending VAF
        return biomarkers + all_variants

    def get_sig(self, sig_name):
        for biomarker in self.varseq_json["biomarkers"]:
            if biomarker["type"] == "SIGNATURE" and biomarker["geneName"] == sig_name:
                return biomarker["quantitativeValue"]
        return ""

    def get_sex(self, sex):
        if sex == "Female" or sex == "Male":
            return sex[0]
        else:
            return "O"

    def get_msi(self): 
        return self.get_sig("MSI")

    def get_tmb(self):
        return self.get_sig("TMB")

    def get_covg_metrics(self):
        covg_summary = self.varseq_json["coverageSummary"]
        bases_20x = covg_summary["basesAt20x"]
        bases_200x = covg_summary["basesAt200x"]
        bases_500x = covg_summary["basesAt500x"]
        covg_mean = round(covg_summary["meanDepth"])
        return [round(x, 2) for x in [bases_20x, bases_200x, bases_500x, covg_mean]]

    def get_custom_field(self, field_name):
        custom_fields = self.varseq_json["customFields"]
        if field_name in custom_fields:
            return custom_fields[field_name]
        else:
            return ""

    def get_pt_name(self):
        # patient names are in the format "Last, First (MRN)"
        name = self.varseq_json["sampleState"]["patientName"].split(" (")[0]
        fn, ln = name.split(",")[:2]
        return fn, ln

    def get_mrn(self):
        return str(self.sample_state["mrn"]).zfill(7)

    def get_prov_name(self):
        name = self.varseq_json["sampleState"]["orderingPhysician"]
        fn, ln = name.split(",")[:2]
        return fn.strip(), ln.strip()

    def get_prov_id(self):
        id =  str(self.get_custom_field("ProviderID"))
        return id.zfill(6)

    def format_header_date(self, date):
        split_date = date.split("/")
        if len(split_date) == 3:
            m, d, y = split_date
            if len(m) == 1: m = "0" + m
            if len(d) == 1: d = "0" + d
            return y + m + d
        else:
            raise Exception(f"Invalid date format: {date} found in JSON for {self.sample_id}, patient name: {self.pt_fn} {self.pt_ln}")

    def get_date(self, date_type):
        date = self.get_custom_field(date_type).split(" ")[0]
        # if date_type.startswith("N_") or date_type.startswith("Date"):
        #     date = self.get_custom_field(date_type).split(" ")[0]
        # else:
        #     date = self.varseq_json["sampleState"][date_type]
        return self.format_header_date(date)

    def get_date_sent(self):
        return date.today().strftime("%Y%m%d")

    def get_tumor_type(self):
        tt = self.varseq_json["tumorType"]
        tt_abbrev = self.varseq_json["tumorTypeAbbrev"]
        return f"{tt} ({tt_abbrev})"

    def get_clin_sig(self, variant):
        if variant in self.varseq_json["biomarkers"]:
            return "Tier 1: Strong significance"
        else:
            return "Tier 2: Potential significance"

    def get_variant_type(self, variant):
        tags = variant["tags"]
        if tags and len(tags) == 1:
            type_initials = tags[0]["initials"]
            if type_initials == "SM":
                return "Somatic"
            elif type_initials == "GL":
                return "Germline"
        return "Unknown"
    
    def insert_newlines(self, s, every=110):
        return '\.br\\'.join(wrap(s, every)).strip()

    def format_interp(self, interp):
        interp = interp.replace("<p>", "").replace("</p>", "").replace("<em>", "").replace("</em>", "").replace("&nbsp;", " ").replace("–", "-")
        return self.insert_newlines(interp)

    def get_interp(self, variant):
        summary = variant["biomarkerSummary"]
        if summary:
            return self.format_interp(summary["interpretation"])
        else:
            return ""

    def get_chrom(self, variant):
        mutation = variant["mutation"]
        return mutation["chr"] if mutation else ""

    def get_naf(self, variant):
        naf = variant["naf"]
        return round(naf, 2) if naf else 0

    def get_consequence(self, variant):
        cons =  SEQ_ONTOLOGY_MAP.get(variant["sequenceOntology"])
        if cons:
            return cons
        else:
            raise RuntimeError(f"SequenceOntology term {variant['sequenceOntology']} not found")

    def get_dna_change(self, variant):
        # count number of bases in ref and alt (excluding dashes) to determine if substitution or indel
        ref, alt = map(lambda b: b.replace('-', ''), variant["refAlt"].split("/"))
        return "Substitution" if len(alt) == len(ref) else "Insertion/Deletion"

    def get_coords(self, variant):
        # VarSeq usually uses 0 based coordinates for the starting position,
        # so we increment "start" to get 1-based coordinates
        # This is not true for insertions, so in that case we increment "stop" for congruency
        start, stop = int(variant['start']), int(variant['stop'])
        if start != stop:
            start += 1
        else: # insertion
            stop += 1
        return str(start), str(stop)

    def get_loinc_info(self, code):
        loinc_info = LOINCS.get(code)
        if loinc_info:
            loinc_name, loinc_type = loinc_info
            return f"{loinc_type}|{code}^{loinc_name}^LOINC"
        else:
            raise RuntimeError("Invalid LOINC code: " + code)

    def create_obx_segment(self, variant_id, loinc_code, value):
        return f"""OBX|{self.next_obx_idx()}|{self.get_loinc_info(loinc_code)}|{variant_id}|{value}"""

    def get_variant_obx_function(self, variant_idx):
        variant_id = get_variant_id(variant_idx)
        def wrapper_func(loinc_code, value):
            return self.create_obx_segment(variant_id, loinc_code, value)
        return wrapper_func

    def get_pdot(self, variant):
        return variant["pDot"] if variant["pDot"] else "p.?"

    def get_variant_display_name(self, variant):
        pdot = self.get_pdot(variant)
        pdot = "" if pdot == "p.?" else f"{pdot} - "
        molecular_cons = self.get_consequence(variant)
        if molecular_cons == "Nonsense":
            molecular_cons = "Nonsense Variant"
        return f"{variant['geneName']} {pdot}{variant['cDot']} {molecular_cons}"

    def get_tumor_variant_obxs(self, variant, variant_idx):
        ref, alt = variant["refAlt"].split("/")
        start, stop = self.get_coords(variant)
        create_obx_segment = self.get_variant_obx_function(variant_idx)
        obxs =  [
            create_obx_segment("48018-6", '^' + variant["geneName"] + '^'), # Variant Name (This is also a discrete field in EPIC)
            create_obx_segment("47998-0", self.get_variant_display_name(variant)), # Variant Display Name
            create_obx_segment("81252-9", 'v1^' + variant['hgvsWithGene'] + '^ClinVar-V'), # Discrete Genetic Variant
            create_obx_segment("48002-0", '^' + self.get_variant_type(variant)), # Genomic Source Class
            create_obx_segment("53037-8", '^' + self.get_clin_sig(variant)), # Genetic Sequence Variation Clinical Significance
            create_obx_segment("81258-6", round(variant["vaf"], 2)), # Allelic Frequency
            create_obx_segment("82121-5", variant["altReadCount"]), # Allelic Read Depth
            create_obx_segment("48005-3", self.get_pdot(variant)), # Amino Acid Change p.HGVS
            create_obx_segment("48019-4", '^' + self.get_dna_change(variant)),
            create_obx_segment("48004-6", '^' + variant['cDot']), # DNA Change c.HGVS
            create_obx_segment("81290-9", variant["gDot"]), # Genomic DNA Change g.HGVS
            create_obx_segment("7400052", variant["proteinId"]),
            create_obx_segment("51958-7", f"{variant['transcriptName']}^{variant['transcriptName']}^RefSeq-T"), # Transcript Reference Sequence
            create_obx_segment("48000-4", self.get_chrom(variant)), # Chromosome
            create_obx_segment("48006-1", '^' + self.get_consequence(variant)),
            create_obx_segment("69547-8", ref), # Genomic Reference Allele
            create_obx_segment("69551-0", alt), # Genomic Alternate Allele
            create_obx_segment("81254-5", f"{start}^{stop}"), # Genomic Allele Start-End
            create_obx_segment("62374-4", '^' + variant['assembly']), # Human Reference Sequence Assembly Version
            create_obx_segment("83005-9", "Simple"), # EPIC Variant Category (Simple, Complex, Fusion, etc.))
            create_obx_segment("69548-6", "Detected"), # Genetic Variant Assessment
        ]
        interp = self.get_interp(variant)
        if interp: obxs.append(create_obx_segment("93364-8", interp))
        zyg = variant.get("zygosity")
        if zyg: obxs.append(create_obx_segment("53034-5", '^' + zyg))
        return "\r\n".join(obxs)

    # While we send these variants in the tumor's HL7 message, we send them again here for the normal sample so we can get the normal allele frequency
    def get_normal_variant_obxs(self, variant, variant_idx):
        create_obx_segment = self.get_variant_obx_function(variant_idx)
        obxs =  [
            create_obx_segment("48018-6", '^' + variant["geneName"] + '^'), # Variant Name (This is also a discrete field in EPIC)
            create_obx_segment("47998-0", self.get_variant_display_name(variant)), # Variant Display Name
            create_obx_segment("81252-9", f"v1^{variant['hgvsWithGene']}^ClinVar-V"), # Discrete Genetic Variant
            create_obx_segment("48002-0", f"^{self.get_variant_type(variant)}"), # Genomic Source Class
            create_obx_segment("53037-8", '^' + self.get_clin_sig(variant)), # Genetic Sequence Variation Clinical Significance
            create_obx_segment("81258-6", self.get_naf(variant)), # Allelic Frequency
            create_obx_segment("48005-3", self.get_pdot(variant)), # Amino Acid Change p.HGVS
            create_obx_segment("48004-6", f"^{variant['cDot']}"), # DNA Change c.HGVS
            create_obx_segment("48000-4", self.get_chrom(variant)), # Chromosome
            create_obx_segment("48006-1", '^' + self.get_consequence(variant)),
            create_obx_segment("62374-4", '^' + variant['assembly']), # Human Reference Sequence Assembly Version
            create_obx_segment("83005-9", "Simple"), # EPIC Variant Category (Simple, Complex, Fusion, etc.))
            create_obx_segment("69548-6", "Detected"), # Genetic Variant Assessment
        ]
        zyg = variant.get("zygosity")
        if zyg: obxs.append(create_obx_segment("53034-5", '^' + zyg))
        return "\r\n".join(obxs)

    def get_tumor_msg_header(self):
        self.reset_obx_idx()
        bases_20x, bases_200x, bases_500x, covg_mean = self.get_covg_metrics()
        return f"""MSH|^~\&|RRH||Beaker||{self.date_sent}||ORU^R01|1|P|2.3||||||\r
PID|1||{self.mrn}^^^MRN^MRN||{self.pt_ln}^{self.pt_fn}^||{self.bday}|{self.sex}\r
ORC|RE\r
OBR|1|{self.order_num}|{self.sample_id}^Beaker|LAB9055^Pan-cancer Solid Tumor Panel^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL|||{self.date_ordered}|||||||||{self.prov_id}^{self.prov_ln}^{self.prov_fn}^^^^^^EPIC^^^^PROVID||||||{self.date_received}|||F\r
{self.create_obx_segment("2a", "7102415", f"^{self.get_tumor_type()}")}\r
{self.create_obx_segment("2a", "81695-9", f"^{self.get_msi()}")}\r
{self.create_obx_segment("2a", "94076-7", f"{self.get_tmb()}")}\r
{self.create_obx_segment("2a", "7102423", bases_20x)}\r
{self.create_obx_segment("2a", "7102424", bases_200x)}\r
{self.create_obx_segment("2a", "7102425", bases_500x)}\r
{self.create_obx_segment("2a", "7102426", covg_mean)}\r"""

    def get_normal_msg_header(self):
        self.reset_obx_idx()
        return f"""MSH|^~\&|RRH||Beaker||{self.date_sent}||ORU^R01|1|P|2.3||||||\r
PID|1||{self.mrn}^^^MRN^MRN||{self.pt_ln}^{self.pt_fn}^||{self.bday}|{self.sex}\r
ORC|RE\r
OBR|1|{self.norm_order_num}|{self.norm_sample_id}^Beaker|LAB9056^Pan-cancer Panel, Comparator^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL|||{self.norm_date_ordered}|||||||||{self.prov_id}^{self.prov_ln}^{self.prov_fn}^^^^^^EPIC^^^^PROVID||||||{self.norm_date_received}|||F\r"""

    def get_tumor_obxs(self):
        return "\r\n".join([self.get_tumor_variant_obxs(variant, idx) for idx, variant in enumerate(self.variants)])

    def get_normal_obxs(self):
        germline_variants = filter(lambda x: self.get_variant_type(x) == "Germline", self.variants)
        return "\r\n".join([self.get_normal_variant_obxs(variant, idx) for idx, variant in enumerate(germline_variants)])

    def get_tumor_msg(self):
        return self.get_tumor_msg_header() + self.get_tumor_obxs()

    def get_normal_msg(self):
        header = self.get_normal_msg_header()
        normal_obxs = self.get_normal_obxs()
        if normal_obxs:
            return header + normal_obxs
        else:
            return ""

def create_hl7_msgs(vs_json):
    vs_info = VarSeqInfo(vs_json)
    return vs_info.get_tumor_msg(), vs_info.get_normal_msg()

def send_hl7_msgs(vs_json):
    with hl7.client.MLLPClient(HOSTNAME, PORT) as client:
        tumor_msg, normal_msg = create_hl7_msgs(vs_json)
        # with open("tumor_msg.txt", "w") as f:
        #     f.write(tumor_msg)
        print(client.send_message(tumor_msg))
        print(f"Sent tumor message: {tumor_msg.splitlines()[3]}")
        if normal_msg:
            # with open("normal_msg.txt", "w") as f:
            #     f.write(normal_msg)
            print(client.send_message(normal_msg))
            print(f"Sent normal message: {normal_msg.splitlines()[3]}")
    return vs_json

app = Flask(__name__)
@app.route('/receivejson', methods=['POST'])
def receive_json():
    data = request.get_json()
    if data:
        result = send_hl7_msgs(data)
        return jsonify(result), 200
    else:
        return jsonify({"error": "No JSON data received"}), 400

# launches Flask server on localhost:5000
if __name__ == '__main__':
    app.run(debug=True)

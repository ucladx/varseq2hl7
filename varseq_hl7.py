import hl7.client
import json
import sys

with open(sys.argv[1], 'r') as f:
    varseq_json = json.load(f)

alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

def _get_variant_id(idx):
    if idx < 0:
        return ""
    quotient, remainder = divmod(idx, len(alphabet))
    # Recursive call for the quotient - 1 to prepend the previous letters if needed
    return _get_variant_id(quotient - 1) + alphabet[remainder]

def get_variant_id(idx):
    return f"2{_get_variant_id(idx)}"

obx_idx = 0
def get_obx_idx():
    global obx_idx
    obx_idx += 1
    return obx_idx

def get_biomarkers():
    biomarkers = filter(lambda x: x["type"] == "VARIANT", varseq_json["biomarkers"])
    return list(biomarkers)

def get_all_variants():
    return get_biomarkers() + varseq_json["germlineVariants"] + varseq_json["uncertainVariants"]

def get_sig(sig_name):
    for biomarker in varseq_json["biomarkers"]:
        if biomarker["type"] == "SIGNATURE" and biomarker["geneName"] == sig_name:
            return biomarker["quantitativeValue"]

def get_msi(): 
    return get_sig("MSI")

def get_tmb():
    return get_sig("TMB")

def get_custom_field(field_name):
    fields = varseq_json["customFields"]
    return fields[field_name] if field_name in fields else ""

def get_prov_id():
    return get_custom_field("providerID")

def get_order_num():
    return get_custom_field("orderNumber")

def get_pt_name():
    return varseq_json["sampleState"]["patientName"].split(" (")[0].split(",")

def get_prov_name():
    return [fn_or_ln.strip() for fn_or_ln in varseq_json["sampleState"]["orderingPhysician"].split(",")][:2]

def format_msh_date(date):
    [m, d, y] = date.split('/')
    return y + m + d

def get_tumor_type():
    tt = varseq_json["tumorType"]
    tt_abbrev = varseq_json["tumorTypeAbbrev"]
    return f"{tt} ({tt_abbrev})"

def get_clin_sig(variant):
    if variant in varseq_json["biomarkers"]:
        return "Tier 1: Strong significance"
    else:
        return "Tier 3: Unknown clinical significance"

def get_variant_type(variant):
    variant_type = variant["tags"][0]["initials"]
    if variant_type == "SM":
        return "Somatic"
    elif variant_type == "GL":
        return "Germline"
    else:
        return "Unknown"

def get_interp(variant):
    summary = variant["biomarkerSummary"]
    if summary:
        interp = summary["interpretation"]
        # TODO if we end up seeing more HTML elements, actually do this properly with re.sub
        interp = interp.replace("<p>", "").replace("</p>", "")
        return interp
    else:
        return ""

def get_chrom(variant):
    mutation = variant["mutation"]
    if mutation:
        return mutation["chr"]
    else:
        return ""

LOINCS = { "48018-6": ["Gene Studied Name", "CWE"],
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
}

def get_loinc_info(code):
    loinc_info = LOINCS.get(code)
    if loinc_info:
        loinc_name = loinc_info[0]
        loinc_type = loinc_info[1]
        return f"{loinc_type}|{code}^{loinc_name}^LOINC"
    else:
        raise RuntimeError("Invalid LOINC code")

def get_variant_obx(variant_id, loinc_code, value):
    return f"""OBX|{get_obx_idx()}|{get_loinc_info(loinc_code)}|{variant_id}|{value}"""

def create_header():
    sample_info = varseq_json["sampleState"]
    sample_id = sample_info["sampleName"]
    mrn = sample_info["mrn"]
    pt_ln, pt_fn = get_pt_name()
    bday = format_msh_date(sample_info["dob"])
    sex = sample_info["sex"]
    order_num = get_order_num()
    prov_id = get_prov_id()
    prov_ln, prov_fn = get_prov_name()
    date_ordered = format_msh_date(sample_info["dateOrdered"])
    date_received = format_msh_date(sample_info["dateReceived"])
    date_collected = format_msh_date(sample_info["dateCollected"])
    date_sent = format_msh_date(sample_info["lastModifiedDate"])
    return f"""MSH|^~\&|RRH||Beaker||{date_sent}||ORU^R01|1|P|2.3||||||\r
PID|1||{mrn}^^^MRN^MRN||{pt_ln}^{pt_fn}^||{bday}|{sex}\r
ORC|RE\r
OBR|1|{order_num}|{sample_id}^Beaker|LAB9055^Pan-cancer Solid Tumor Panel^BKREAP^^^^^^SOLID TUMOR PAN-CANCER PANEL|||{date_ordered}|||||||||{prov_id}^{prov_ln}^{prov_fn}^^^^^^EPIC^^^^PROVID||||||{date_received}|||F\r
{get_variant_obx("2a", "7102415", f"^{get_tumor_type()}")}\r
{get_variant_obx("2a", "81695-9", f"^{get_msi()}")}\r
{get_variant_obx("2a", "94076-7", f"{get_tmb()}")}\r"""

def get_variant_obx_function(variant_id):
    def wrapper(loinc_code, value):
        return get_variant_obx(variant_id, loinc_code, value)
    return wrapper

def get_variant_obxs(variant, idx):
    display_name = variant["geneName"]
    variant_id = get_variant_id(idx)
    get_variant_obx = get_variant_obx_function(variant_id)
    obxs =  [
        get_variant_obx("47998-0", display_name), # Variant Display Name
        get_variant_obx("83005-9", "Simple"), # EPIC Variant Category (Simple, Complex, Fusion, etc.))
        get_variant_obx("81252-9", f"v1^{variant['hgvsWithGene']}^ClinVar-V"), # Discrete Genetic Variant
        get_variant_obx("48018-6", f"^{display_name}^"), # Variant Name (This is also a discrete field in EPIC)
        get_variant_obx("48005-3", variant["pDotThreeLetter"]), # Amino Acid Change p.HGVS
        get_variant_obx("48004-6", f"^{variant['cDot']}"), # DNA Change c.HGVS
        get_variant_obx("48006-1", f"^{variant['effect']}"), # Molecular Consequence (Missense, Nonsense, etc.)
        get_variant_obx("62374-4", "^GRCh38"), # Human Reference Sequence Assembly Version
        get_variant_obx("53037-8", f"^{get_clin_sig(variant)}"), # Genetic Sequence Variation Clinical Significance
        get_variant_obx("48000-4", get_chrom(variant)), # Chromosome
        get_variant_obx("81258-6", round(variant["vaf"], 2)), # Allelic Frequency
        get_variant_obx("82121-5", variant["altReadCount"]), # Allelic Read Depth
        get_variant_obx("69548-6", "Detected"), # Genetic Variant Assessment
        get_variant_obx("93364-8", get_interp(variant)), # Genetic Variant Diagnostic Significance
        get_variant_obx("48002-0", f"^{get_variant_type(variant)}"), # Genomic Source Class
    ]
    return "\r\n".join(obxs) + "\r\n"

def create_hl7_msg():
    msg = create_header()
    for idx, variant in enumerate(get_all_variants()):
        msg += get_variant_obxs(variant, idx)
    return msg

def send_hl7_msg():
    msg = create_hl7_msg()
    # prints the response from the server
    print(hl7.client.MLLPClient("interface-test.mednet.ucla.edu", 7199).send_message(msg))
    with open("./hl7_msg.txt", "w") as f:
        f.write(msg)
    # print(msg)

send_hl7_msg()

# if __name__ == "__main__":
#     import sys
#     if len(sys.argv) != 2:
#         print("Usage: python3 varseq_hl7.py <path to varseq json>")
#         sys.exit(1)
#     send_hl7_msg(sys.argv[1])

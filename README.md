### Goal

Convert discrete genomic results in JSONs from [VarSeq VSClinical](https://www.goldenhelix.com/blog/versatile-ngs-exports-with-vsclinicals-custom-scripts/) to HL7v2 messages compatible with EPIC's [DiscreteGenomicResults](https://open.epic.com/Interface/HL7v2#DiscreteGenomicResults) interface.

### Scope

- Host a REST API to receive VSClinical JSONs containing discrete data on variants.
- Maintain a custom VSClinical report script (`report_cancer.js`) that enables VarSeq to send this JSON to the API.
- Map discrete variant data from the VSClinical JSON to LOINC codes supported by EPIC.
- Use the JSON and this mapping to construct and send an HL7v2 message to a pre-configured [EPIC interface](https://open.epic.com/Interface/HL7v2).

### Usage
Launch the Flask server which hosts the REST API:

`python3 varseq_hl7.py`

When this server receives a POST request with a VSClinical JSON, it will send an HL7v2 message representing the JSON's data to the HOSTNAME and PORT variables set in `varseq_hl7.py`

### Custom Report Script
To automate sending a JSON to this server from VSClinical, we use a customized report script.

This [GoldenHelix blog post](https://www.goldenhelix.com/blog/customizing-vsclinical-reports-with-new-output-types/) contains a tutorial on creating a custom report script.

`report_cancer.js` is an example custom report script. When a report is rendered in VarSeq using this script, it sends its JSON to the Flask server (by default assumed to be hosted at http://localhost:5000/)
### Goal

Convert discrete genomic results in JSONs from [VarSeq VSClinical](https://www.goldenhelix.com/blog/versatile-ngs-exports-with-vsclinicals-custom-scripts/) to [HL7v2 messages](https://open.epic.com/Interface/HL7v2#DiscreteGenomicResults) compatible with EPIC.

### Scope

- Host a REST API to receive VSClinical JSONs containing discrete data on variants.
- Maintain a custom script that enables VarSeq instances to send this JSON.
- Map discrete variant data from VSClinical specs to specs for EPIC Genomic Indicators.
- Package data into an HL7 message and send to a pre-configured EPIC interface.

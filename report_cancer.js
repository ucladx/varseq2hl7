// import { exec } from '/tmp/child_process.js';

function create_hl7_msg(reportData) {
	const output = exec("python3 ~iatol/hl7/create_hl7_msg.py " + JSON.stringify(reportData));
	console.log("output: " + output);
};

export default {
    description: 'Add NAF and send an HL7 message to Beaker',

    inputFileNamesTemplate: ['{sampleName}_enhanced.json'],
    outputFileNamesTemplate: ['{sampleName}_hl7.json'],

    async render({
        templatedInputs,
        api
    }) {
        const [sampleReportDataString] = templatedInputs;
        const reportData = JSON.parse(sampleReportDataString);
        console.log("Rendering report in Beaker...");
        // Add NAF and NDP to variants
        for(const variant of reportData.biomarkers.concat(reportData.germlineVariants.concat(reportData.uncertainVariants))) {
            variant.naf = null;
            variant.ndp = null;
            const { projectTableUuid, projectTableRecordId} = variant;
            if (!projectTableUuid) {
                continue;
            }

            const sources = await api.projectTableSources({uuid: projectTableUuid, algKey: 'variants'});
            if (sources.length<1){
                continue;
            }
            const variantInfoSource = sources[0];
            const recordNAF = await api.projectTableRecords({
                sourceUrl: variantInfoSource.url,
                recordId: projectTableRecordId,
                sampleId: reportData.sampleId,
                fieldSymbols: ['NAF']
            });
            const recordNDP = await api.projectTableRecords({
                sourceUrl: variantInfoSource.url,
                recordId: projectTableRecordId,
                sampleId: reportData.sampleId,
                fieldSymbols: ['NDP']
            });
            const NAF = recordNAF.find(r => r != null && r.length > 0 && r[0] != null);
            const NDP = recordNDP.find(r => r != null && r.length > 0 && r[0] != null);
            if(NAF != null) {
                variant.naf = NAF[0][0];
            }
            if(NDP != null) {
                variant.ndp = NDP[0];
            }
        }

        // OS call a Python script that creates and sends the HL7 message
        // create_hl7_msg(reportData);
        return [
            JSON.stringify(reportData)
        ];
    }
};
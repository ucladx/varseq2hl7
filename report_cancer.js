async function sendJSON(jsonData, url) {
    const response = await fetch(url, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(jsonData)
    });

    if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
    }

    return await response.json();
}


export default {
    description: 'Add SV VAF (if present) to variants, then send an HL7 message to Beaker',
    inputFileNamesTemplate: ['{sampleName}.json'],
    outputFileNamesTemplate: ['{sampleName}_enhanced.json'],

    async render({
        templatedInputs,
        api
    }) {
        const [sampleReportDataString] = templatedInputs;
        const reportData = JSON.parse(sampleReportDataString);
        console.log("Rendering report in Beaker...");
        // Add SV VAF to variants
        for(const variant of reportData.biomarkers.concat(reportData.germlineVariants.concat(reportData.uncertainVariants))) {
            variant.sv_vaf = null;
            variant.sv_reads = null;
            const { projectTableUuid, projectTableRecordId} = variant;
            if (!projectTableUuid) {
                continue;
            }

            const sources = await api.projectTableSources({uuid: projectTableUuid, algKey: 'variants'});
            if (sources.length<1){
                continue;
            }
            const variantInfoSource = sources[0];
            const recordSR = await api.projectTableRecords({
                sourceUrl: variantInfoSource.url,
                recordId: projectTableRecordId,
                sampleId: reportData.sampleId,
                fieldSymbols: ['SR']
            });
            const SR = recordSR.find(r => r != null && r.length > 0 && r[0] != null);
            if(SR != null) {
                const split_reads = String(SR[0]).split(',');
                const wt_reads = parseInt(split_reads[0]);
                const sv_reads = parseInt(split_reads[1]);
                variant.sv_reads = sv_reads;
                variant.sv_vaf = sv_reads / (wt_reads + sv_reads);
            }
        }

        // Send the JSON to a Flask server that will send a corresponding HL7 message to Beaker
        const url = 'http://localhost:5000/receivejson';
        sendJSON(reportData, url)
            .then(data => console.log(data))
            .catch(error => console.error(`Error: ${error}`));
        return [
            JSON.stringify(reportData)
        ];
    }
};
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
    description: 'Add NAF and send an HL7 message to Beaker',
    inputFileNamesTemplate: ['{sampleName}.json'],
    outputFileNamesTemplate: ['{sampleName}_enhanced.json'],

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
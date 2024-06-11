// @ts-check

export default {
    description: 'Import variants tagged SM/GL and add TMB/MSI',
    icon: 'import',
    fileFilter: undefined,

    async run({
        selected, inputs, api: { table, project, cancer, core }, args,
    }) {
        try {
            // args as string[] input list, extended by vspipeline run command context
            const [evalId] = args;
            const evaluationId = parseInt(evalId);
			const sample = await project.currentSampleState();
			await addTMB(sample, cancer, evaluationId);
			await addMSI(sample, cancer, evaluationId);

            console.log("Adding variants tagged SM and GL...");
            const variantSomaticIds = await cancer.addProjectVariants({
                origin: "SOMATIC",
                recordSetInitials: "SM",
                maxRecordCount: 500,
                evaluationId: evaluationId
            });

            const variantGermlineIds = await cancer.addProjectVariants({
                origin: "GERMLINE_CONFIRMED",
                recordSetInitials: "GL",
                maxRecordCount: 500,
                evaluationId: evaluationId
            });

            console.log(`Added ${variantSomaticIds.length} somatic and ${variantGermlineIds.length} germline variants.`);

//            const cnvBiomarkerIds = await cancer.addProjectCnvs({
//                recordSetInitials: "PF",
//                maxRecordCount: 500,
//                evaluationId: evaluationId
//            });
//            console.log(`Added ${cnvBiomarkerIds.length} CNVs.`);

        } catch (error) {
            console.log(error);
            throw error;
        }
    }
};

async function addTMB(sample, cancer, evaluationId) {
    try {
		const [tmbScoreField = {}] = sample.extraField.filter(f => f.fieldName == "TMB");
		let { fieldValue: tmbScore } = tmbScoreField;

        if (tmbScore === undefined) {
            console.log(`Failed to parse TMB: couldn't find values for ${tmbScore}`);
            return;
        }

        let status = 'Low';
		if (parseFloat(tmbScore) >= 20.0) {
			status = 'High';
		}

		const units = 'mut/Mb'
        console.log(`Adding TMB with status ${status} for value ${tmbScore}`);
        await cancer.addSignature({
            signature: 'TMB',
            status,
            quantitativeValue: tmbScore,
            units,
            evaluationId
        });
    } catch (error) {
        console.error(`Failed to parse TMB: ${error.toString()}`);
    }
}

async function addMSI(sample, cancer, evaluationId) {
    try {
		const [msiScoreField = {}] = sample.extraField.filter(f => f.fieldName == "MSI");
		let { fieldValue: msiScore } = msiScoreField;

        if (msiScore === undefined) {
            console.log(`Failed to parse MSI: couldn't find values for ${msiScore}`);
            return;
        }

        let status = 'Stable';
		const msiScoreFloat = parseFloat(msiScore);
		if (msiScoreFloat <= 6.5) {
			status = 'Stable';
		}
		else if (msiScoreFloat >= 8.5) {
			status = 'High';
		}
		else {
			status = 'Indeterminate';
		}

		const units = 'sites'
        console.log(`Adding MSI with status ${status} for value ${msiScore}`);
        await cancer.addSignature({
            signature: 'MSI',
            status,
            quantitativeValue: msiScore + '%',
            units,
            evaluationId
        });
    } catch (error) {
        console.error(`Failed to parse TMB: ${error.toString()}`);
    }
}
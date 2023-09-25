// Test Import All Variants
// Creates a new evaluation and adds all variants from the project to it, then creates a JSON
let project_path = process.argv[2];
let output_path = project_path + "/report";
const JSON_SCRIPT_NAME = 'PCP Report Script v1';

const { evaluationId } = await cancer.createEvaluation({});
const tumorType = await cancer.getTumorType({ evaluationId });
console.log(`Created evaluation with ID ${evaluationId} and tumor type ${tumorType.tumorType}.`);

console.log("Adding project variants...");
const variantBiomarkerIds = await cancer.addProjectVariants({
	evaluationId: evaluationId
});
console.log(`Added ${variantBiomarkerIds.length} variants.`);

console.log("Creating JSON...");

await cancer.runReportScript({ evaluationId, templateName: JSON_SCRIPT_NAME, outputPath: output_path });

console.log("Done creating JSON.");

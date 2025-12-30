export function getAggrCsvConfigSpec() {
    return {
        format: "Cell Ranger aggr config CSV for 10x Cloud",
        description: "Aggregation CSV for combining multiple Cell Ranger count or multi analyses into a single feature-barcode matrix and analysis.",
        structure: "Simple CSV with headers",
        columns: {
            sample_id: {
                description: "A unique identifier for this sample.",
                required: true,
                type: "string",
                example: "Sample1",
                notes: [
                    "For count analysis, there is only one sample and its ID is the analysis name.",
                    "For multi analyses, there can be multiple samples and the ID can be inferred from the analysis files. e.g. For this file: per_sample_outs\\Sample2\\count\\sample_cloupe.cloupe the sample_id would be Sample2. Note you must use 'Sample2' without modifications such as changing case or adding prefixes. Sample names do not need to be unique across analyses.",
                ],
            },
            molecule_h5: {
                description: "Analysis ID in txg://analyses/<analysisID> format. This references the molecule_info.h5 file from a completed Cell Ranger count or multi analysis.",
                required: true,
                type: "string",
                example: "txg://analyses/abc123xyz",
                notes: [
                    "Must reference a completed analysis",
                    "Use 'list_analyses' tool to find analysis IDs",
                    "All analyses must use the same references",
                ],
            },
        },
        rules: [
            "First row must contain headers: sample_id,molecule_h5",
            "Each subsequent row represents one sample to aggregate",
            "Minimum 2 samples required for aggregation",
        ],
        examples: {
            basic_aggr: {
                description: "Basic aggregation of three analyses",
                csv: `sample_id,molecule_h5
sample1,txg://analyses/analysis_id_1
sample2,txg://analyses/analysis_id_1
sample3,txg://analyses/analysis_id_2
sample4,txg://analyses/analysis_id_3`,
            },
        },
    };
}

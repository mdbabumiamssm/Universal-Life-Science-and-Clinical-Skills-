# 10x Genomics Cloud MCP Server

## Overview

This extension enables Claude to interact with 10x Genomics Cloud platform, allowing users to launch, monitor, and manage single-cell analysis pipelines directly through conversation. Users can execute Cell Ranger analysis, upload data, track pipeline progress, retrieve results, and perform quality control checks.

## Prerequisites

### Software requirements

- 10x Cloud account and 10x Cloud access token

### Data requirements

These data inputs are needed for analyzing data on the 10x Cloud:

- **Data**: FASTQ files from 10x Chromium single cell data from these products:
  - Universal 3' Gene Expression
  - Universal 5' Gene Expression
  - Flex Gene Expression
- **Optional**: custom reference (alternatively, use pre-built human and mouse references)
- **Optional**: multi config CSV for multi pipeline, Feature Barcode reference for count pipeline, aggregation CSV for the aggr pipeline (alternatively, the LLM can help you create the CSV files, but may need refinement to ensure the parameters and formatting are correct)

## MCP server set up instructions

### 1. Set up a 10x Cloud account

First-time users of the 10x Cloud platform must create a new account. If you already have a 10x Cloud account, there is no need to create a new one.

### 2. Download 10x Genomics Cloud extension

There are two ways to access the server:

**Option 1**: In Claude Desktop, go to Settings > Extensions > Browse extensions > Look for the 10x Genomics Cloud. Click "Install" and follow the prompts.

**Option 2**: Download the latest 10x Genomics Cloud Analysis MCP server for Claude Desktop from the [Releases page](https://github.com/10XGenomics/txg-mcp/releases). Double-click the .mcpb file downloaded from GitHub. Alternatively, go to Settings > Extensions > Drag and drop the .mcpb file to install. An installation pop-up should open in Claude Desktop. Then click "Install" and follow the prompts.

### 3. Enter your 10x Cloud Analysis Access Token

Your access token can be found in the Security section of your Cloud Analysis Account Settings.

### 4. Enable the extension

Switch the toggle from "Disabled" to "Enabled".

You should now be able to use the 10x Cloud Analysis MCP server in Claude chats.

### Optional: Claude Desktop Extensions

For operations that involve uploading files to the Cloud (e.g., CSV files), you may have an improved experience with the MCP server workflow by connecting additional Desktop extensions, such as the **Filesystem** extension. These are available in Settings > Extensions > Desktop extensions. Configure on install, or go to Settings > Extensions > Configure to toggle "Enabled" and set the allowed directories that the filesystem server can access.

## Example data

To test the MCP server, you can select from a number of single cell datasets on the [10x Genomics public datasets website](https://www.10xgenomics.com/datasets) and compare the MCP server-created analysis results with the posted outputs on the datasets website (e.g., web_summary.html).

## LLM prompt tips and best practices

Before we get started with use cases, here are some tips for constructing chat prompts. It is important to be specific in your requests and queries and encourage the LLM to explain the decision-making process. The LLM can help you get started, but as with a normal conversation, it will likely require iteration to successfully create, run, and analyze a dataset.

### Include specific details in your prompts

- **Dataset details**: Include information about the dataset species (e.g., human, mouse), sample type (e.g., cells or nuclei, disease state), and experimental setup (e.g., sample preparation, replicates, experiment conditions). These details are useful when determining how to analyze and interpret the data.

- **Product and modality information**: Include the data product family (e.g., Universal 3' Gene Expression, Universal 5' Gene Expression, Flex Gene Expression), modality (e.g., gene expression, cell multiplexing, antibody, antigen, CRISPR), and sample plexy (single vs. multiple samples). It is also helpful to indicate the library type (e.g., "GEM-X 3' Gene Expression v4") if known.

- **Pipeline and parameters**: Include what pipeline to run (e.g., count, multi, aggr), as well as specific parameters, if any (e.g., create BAM files for 3' and 5' Gene Expression, but not for Flex Gene Expression analyses). Or if unsure, indicate in the prompt that you would like assistance in selecting the optimal pipeline and parameters. The include_introns parameter is set to true and chemistry detection is set to auto by default.

## Use cases: Analysis with the MCP server

Below, we provide example prompts and follow-up questions that you can use as templates to get started with common 10x Cloud and data analysis use cases.

### General MCP server questions and cloud project management

The MCP server has many functions, including the following:

- Create a project
- Upload input files (e.g., FASTQ, CSV)
- Create an analysis and run it on 10x Cloud
- Check analysis run status
- Download analysis output files to your computer
- List project, analysis, files, pre-built references, annotation models
- Check multi config CSV specifications

Go to your Cloud Account to further manage projects and analyses. These functions are not yet supported with the MCP server:

- Cancel an in-progress analysis
- Delete an analysis

This is a good prompt to begin with following the initial extension setup:

```
What 10x cloud tools do you have access to?
```

Example follow-up prompts:

```
What project settings can I update via this MCP server?
```

```
List pre-built transcriptome references.
```

```
List the annotation models available for cell annotation
```

```
List the current analyses in all my 10x cloud projects.
```

```
List the project files of the current project.
```

```
List the analysis outputs from the "3prime-count" analysis.
```

### Set up a Cell Ranger count analysis

```
Create a new project called: "first-mcp-project". The analysis I want to run should be called: "3p-GEX-count". I have 1 sample of single cell Universal 3' Gene Expression human data with 2 lanes of sequencing. The library and chemistry type is NextGEM 3' Gene Expression v3.

The FASTQ files to upload to the Cloud are in this folder: /Users/<user.name>/Desktop/3pGEX-count. I want to run the Cell Ranger count pipeline.

Use the support documentation from these websites to help me set up an analysis: https://www.10xgenomics.com/support/software/cloud-analysis/latest and https://www.10xgenomics.com/support/software/cell-ranger/latest. Are there any parameters for the count pipeline that I should consider for my analysis?
```

Example follow-up questions:

```
Could you check the analysis run status?
```

```
What is the total file size of the outputs from this analysis?
```

```
Could you download all the outputs to my computer? Save them to my Desktop/.
```

### Set up a Cell Ranger multi analysis

In addition to creating the multi analysis, you can also learn about the configuration options directly in your chat:

```
How do I create a multi config csv for a cellranger multi analysis?
```

An example with a detailed prompt. The level of detail you provide is subjective, but more information may improve or speed up the analysis creation process. This particular example includes information provided with this 10x public dataset.

```
The analysis I want to run should be called: "Flex-GEX-multi", and I want it added to an existing project called "first-claude-mcp-CR-9-flex". I have 1 sample of single nuclei Flex Gene Expression human data with 2 lanes of sequencing. The library and chemistry type is GEM-X Flex Gene Expression.

Fresh frozen human kidney tissue was obtained from Avaden BioSciences by 10x Genomics. The tissue was sectioned into two 44 mg samples and nuclei were isolated using the Chromium Nuclei Isolation Kit (CG000505). A total of 4,180,000 nuclei were obtained, which were aliquoted into a sample containing approximately 2 million nuclei. The sample was fixed for one hour at room temperature using the Fixation of Cells & Nuclei for GEM-X Flex Gene Expression protocol (CG000782). Following fixation, probe hybridization was performed according to the user guide. 300,000 nuclei aliquots were used for hybridization. After hybridization, the samples were washed and filtered according to the protocol. The sample was then resuspended, counted, and loaded into a Chromium GEM-X Chip. 4,000 nuclei were targeted. Gene Expression libraries for GEM-X Flex were generated using the GEM-X Flex Gene Expression Reagent Kit for Singleplex samples (CG000786). Libraries were sequenced on an Illumina NovaSeq 6000, with a mean read depth of approximately 20,000 reads per cell using a paired-end, dual indexing sequencing scheme.

The FASTQ files are already uploaded to this project: "first-claude-mcp-CR-9-flex". I want to run the Cell Ranger multi pipeline. Use the support documentation from these websites to help me set up an analysis: https://www.10xgenomics.com/support/software/cloud-analysis/latest and https://www.10xgenomics.com/support/software/cell-ranger/latest. Are there any parameters for the multi pipeline that I should consider for my analysis?
```

Example of setting up a more complicated Cell Ranger multi analysis, in addition to iterating on the analysis with the same inputs but changing some parameters.

```
I want to set up another cellranger multi analysis in the mcp-project-30sept25 project. The analysis should be called "3pgex-cmo-crispr-multi". This dataset is Universal 3' Gene Expression v3 data with a Cellplex (Cell Multiplexing) library and a CRISPR library. There is 1 sample, and 3 feature types: Gene Expression, Multiplexing Capture, and CRISPR Guide Capture. The CellPlex experiment uses 4 CMO tags (CMO308, CMO309, CMO310, and CMO311). This is a human sample, so I want to use a pre-built human transcriptome reference.

The FASTQ files for all these libraries are located in: /Users/<user.name>/Desktop/3pgex-crispr-cmo. The feature reference is in the same path, called "feature_ref.csv". What analysis considerations or parameters should I think about before running this analysis?
```

Note that analysis iteration within Claude should work on analyses created by the MCP, but not analyses created by other methods (e.g., directly in the Cloud web app or TXG CLI).

```
Can you start a new run with the same 3 libraries (GEX, CMO, and CRISPR), but this time treat it as 1 sample pooled from 4 tubes. In this case, the samples part of the config CSV should look like this:

[samples]
sample_id,cmo_ids,description
gex_crispr_cmo_4tags,CMO308|CMO309|CMO310|CMO311,4tags
```

### Explain inputs, parameters, and pipelines

> **Note**: Think of these conversations as starting points. However, it is important to check the original resources to verify analysis recommendations and data interpretation.

Example prompts:

```
Can you explain what the R1, R2, I1, and I2 FASTQ files are and provide the resources used to compile this information?
```

```
What is the difference between the Cell Ranger count and multi pipelines? Which is best for Flex Gene Expression data analysis? Provide the resources used to compile this information.
```

```
What are the main pipeline analysis steps used to analyze my data with the count pipeline? Provide the resources used to compile this information.
```

### Create a batch of analyses for multiple samples

> **Tip**: Set up the Filesystem Desktop Extension for this use case, as it makes it easier for Claude to find the files.

```
I have many samples of data on my computer at /Users/<path>/many_samples. Can you upload all of these into a new project and start cellranger count analyses for all of them at the same time?
```

### QC and interpreting analysis results

> **Note**: Think of these conversations as starting points. However, it is important to check the original resources to verify analysis recommendations and data interpretation.

> **Note**: Claude Desktop cannot read downloaded HTML files > 1MB. However, you can ask the LLM to read other output CSV files or type metric results into the chat and ask for interpretation guidance.

```
Read the metrics_summary.csv file and explain the results.
```

```
Read the cell_types.csv file. Which are the most common cell types in my dataset?
```

## Feedback and troubleshooting

If you have feedback about this feature or need help troubleshooting, reach out to support@10xgenomics.com.

## Privacy Policy

https://www.10xgenomics.com/legal/privacy-policy

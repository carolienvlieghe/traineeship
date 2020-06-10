# Electronic traineeship notebook 
## General info repository
This repository contains different scripts used for automated analysis of flow cytometry data.
## Pipeline
### [csv2fcs](/csv2fcs.R)
This script converts csv files to fcs files (flow cytometry standard).

The input data for the experiment is a csv file created in Kaluza: after pre-analysis (compensation, excluding debris, etc.) use the 'Export compensated event data' option in Kaluza to save the (compensated!) events of interest in a csv file. For this experiment we will export the events of a purified CD19 positive cell population (B-lineage).

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window will open where you can select multiple csv files to convert.
- A second dialog window opens where you can select or create the folder to save the converted fcs files.

The newly created fcs files are to be found in your folder of choice.

### [normalization](/normalization.R)
This script normalizes the data from the normal bone marrow samples (NBM) before we merge them into one reference NBM sample. Note that this script was written specifically for the B-ALL tube. To apply this to other tubes, the reference values of the negative population for each marker (reference.neg.pop.mode) should be adjusted.

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window will open where you can select multiple fcs files to normalize.
- A second dialog window opens where you can select or create the folder
to save the normalized fcs files.
- A pop up will appear, telling you what do to do in the next step. 

One fcs file after another is automatically processed. In a first step, the linear data is transformed into logicle with the logicleTransform() function. A plot window opens, containing a density plot for each marker to check the logicle transformation. If the transformation looks incorrect, it might help to adjust the parameter 'w' within the function.
After the logicle transformation, per marker a plot window opens where you can select the peak of the negative population with the cursor. The blue line is the value of the reference negative population. Note that for CD81, CD38 and CD19 there is no true negative population. Here the peak of positive population is selected. For CD45 the mid peak needs to be selected.

After going through these steps, the normalized fcs files are to be found in your folder of choice.

### [FCSmerge](/FCSmerge.R)
This script merges fcs files into one fcs file. It can be used e.g. to merge the normalized normal bone marrow samples into one reference NBM file.

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window will open where you can select multiple fcs files to merge.
- A pop up will appear, asking you if you want to select another file (this way you are able to select files from multiple folders). If you choose yes, the previous dialog window will open and you can select more fcs files. If you choose no, the script proceeds to the next step.
- A dialog window opens to save your merged fcs file: select a folder, name the file and click the 'Save' button.

### [FCS_tag_merge](/FCS_tag_merge.R)
This script was created to be able to give fcs files a tag and then merge them together into one fcs file. The tag allows you to identify each dataset in post analysis. The tags in this script are numeric, the value corresponds with a certain sample type:
- Diagnosis (Dg) = 250
- Follow up (FU) = 500
- Normal bone marrow (reference) = 750

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window opens, asking you to select the NBM file: select the file or press 'Cancel' if no NBM file needs to be merged.
- A dialog window opens, asking you to select the FU file: select the file or press 'Cancel' if no FU file needs to be merged.
- A dialog window opens, asking you to select the Dg file: select the file or press 'Cancel' if no Dg file needs to be merged.
- A dialog window opens to save your tagged & merged fcs file: select a folder, name the file and click the 'Save' button.

The easiest way to identify the different datasets in Kaluza is to add a histogram plot on the parameter 'TAG', there you can auto-gate the 3 different datasets and name them accordingly.

### [frozen_flowSOM_NBM](/frozen_flowSOM_NBM.R)
This script creates 24 different flowSOM representations of a reference NBM. These representation are revised and the best representation is 'frozen' for analysis of diagnostic and follow up samples.

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window opens, asking you to select the reference NBM file.
- A second dialog window opens where you can select or create the folder to save the data from each of the 24 flowSOM runs. One run contains an fcs file with the flowSOM coordinates included, a pdf file with the flowSOM MST & grid (metaclusters included) and an Rdata file containing the flowSOM result. This last one is used to plot the new data from the diagnosis and follow up samples on the same MST (FROZEN!).
- A plot window opens, containing a density plot for each marker to check the logicle transformation. As described above, the 'w' parameter can be adjusted when the transformation looks incorrect.

Afterwards, the different representations can be revised in Kaluza. The fcs files contain 3 more parameters:
- flowSOM.X: containing the X-coordinates of the MST for every event.
- flowSOM.Y: containing the Y-coordinates of the MST for every event.
- Metaclustering Consensus: containing the metacluster number for every event.

A good starting point to identify populations, is to use the Metaclustering Concensus in a histogram. Every metacluster can be gated and named accordingly after analysis of the metaclusters. The flowSOM MST is represented by plotting flowSOM.X & flowSOM.Y.

A few important variables in the script used as arguments for flowSOM parameters that can be adjusted are:
- parameters.to.use: a vector with the channels that are used in the flowSOM() function. The vector contains the column numbers corresponding with the parameters that need to be used. Select the channels that may contribute to identify cell subsets. 
- x.dim & y.dim: these values refer to the grid size, e.g. setting both x.dim and y.dim is set to 11, results in a grid containing 121 nodes (11 rows and 11 columns). Decreasing these values will result in a smaller grid, hence a more simplified MST but could cluster rare population together with other another population that is actually different. Where increasing these values will purify the population more but could separate cells that are similar into different clusters. Note that it is desired to overcluster a little bit, this can be compensated by the metaclustering value.
- nb.metacluster: a number representing the expected number of metaclusters often based on prior knowledge.

### [frozen_dg_fu_nbm](/frozen_dg_fu_nbm.R)
This script plots new data onto the frozen flowSOM. Note that it is advisable to run a 'free' flowSOM as well, since adding new data to an existing flowSOM object will map all events to an existing cluster. Previous to running this script, the 'tag & merge fcs files' is used to merge datasets from diagnostic and follow up samples.

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window opens, asking you to select the merged fcs file.
- A second dialog window opens, asking you to upload the Rdata file. This is the file containing the information to recreate the frozen flowSOM analysis.
- A plot window opens, containing a density plot for each marker to check the logicle transformation. As described above, the 'w' parameter can be adjusted when the transformation looks incorrect.
- A dialog window opens to save your flowSOM analysis: select a folder, name the analysis and click the 'Save' button.

Two files are created:
- A pdf file containing the MST of the diagnosis sample, the follow up sample and merged NBM sample.
- A fcs file containing 3 additional parameters: flowSOM.X, flowSOM.Y and metaclustering consensus (see above for more explanation)

For the post analysis in Kaluza, it is easiest to make a histogram with the tag parameter to distinguish between the different datasets (see above) and a second histogram with the metaclustering consensus to identify the cell subsets.

### [free_flowSOM](/free_flowSOM.R)
This script runs a 'free' flowSOM on the data and is used in addition to the frozen flowSOM script.

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window opens, asking you to select the fcs file you want to process.
- A plot window opens, containing a density plot for each marker to check the logicle transformation. As described above, the 'w' parameter can be adjusted when the transformation looks incorrect.
- A dialog window opens to save your flowSOM analysis: select a folder, name the analysis and click the 'Save' button.

Two files are created:
- A pdf file containing the MST of the diagnosis sample, the follow up sample and merged NBM sample.
- A fcs file containing 3 additional parameters: flowSOM.X, flowSOM.Y and metaclustering consensus (see above for more explanation).

### [tSNE_flowSOM_umap](/tSNE_flowSOM_umap.R) & [tSNE_flowSOM](/tSNE_flowSOM.R)
These script are based on the original template provided by Beckman Coulter. The first script runs tSNE, flowSOM & umap simultaneously on an fcs file, the second one only runs tSNE & flowSOM. Note that the tSNE and umap analysis have a remakably longer run time.

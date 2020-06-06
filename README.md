# Electronic traineeship notebook 
## General info repository
This repository contains different scripts used for the analysis of flow cytometry data.
The .Rtemplate files were provided by Beckman Coulter and can be used in the Rplugin in Kaluza. They can be viewed in a text editor like visual studio code, notepad++, ... 
The templates are used as a starting point to write scripts for a comparable analysis in Rstudio (independent of the plugin and Kaluza).
## Pipeline
### csv2fcs
This script converts csv files to fcs files (flow cytometry standard).

The input data for the experiment is a csv file created in Kaluza: after pre-analysis (compensation, excluding debri, etc.) use the 'Export compensated event data' option in Kaluza to save the (compensated!) events of interest in a csv file. For this experiment we will export the events of a purified CD19 positive cell population (B-lineage).

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window will open where you can select multiple csv files to convert.
- A second dialog window opens where you can select or create the folder to save the converted fcs files.

The newly created fcs files are to be found in your folder of choice.

### preprocessing_normalization
This script normalizes the data from the normal bone marrow samples (NBM) before we merge them into one reference NBM sample. Note that this script was written specifically for the B-ALL tube. To apply this to other tubes, the reference values of the negative population for each marker (reference.neg.pop.mode) should be adjusted.

In Rstudio press the 'Source' button or press 'Ctrl+Shift+Enter' to run the entire script:
- A dialog window will open where you can select multiple fcs files to normalize.
- A second dialog window opens where you can select or create the folder to save the normalized fcs files.
- A pop up will appear, telling you what do to do in the next step. 

One fcs file after another is automatically processed. In a first step, the linear data is transformed into logicle with the logicleTransform() function. A plot window opens, containing a density plot for each marker to check the logicle transformation. If the transformation looks incorrect, it might help to adjust the parameter 'w' within the function.
After the logicle transformation, per marker a plot window opens where you can select the peak of the negative population with the cursor. The blue line is the value of the reference negative population. Note that for CD81, CD38 and CD19 there is no true negative population. Here the peak of positive population is selected. For CD45 the mid peak needs to be selected.

After going through these steps, the normalized fcs files are to be found in your folder of choice.



Scripts from my internship @ BiGCaT, UM, Maastricht about network inference using breast cancer data and in silico data from HPN-DREAM Network Inference Challenge.

Based on the information on Wiki, relations are defined. Relations link three components: intervention case (stimuli and inhibitors), no-intervention case (only stimuli) and target node. Also, targets are defined according to inhibitors targeting them.

Time series data for each case (intervention/no-intervention) is obtained from the entire data. There are compared for each relation by calculating the area under time series curves. If intervention causes a decrese in the expression of child node, this shows the child node is activated by the target node. Or, if intervention causes an increase in the expression of child node, this shows the child node is inhibited by the target node.

Edges are relatively scored. The differences obtained from each case (intervention/no-intervention); the maximum difference is found among them; all values are divided by the maximum value and these scores (if they are over 0.1) are included in EDA file with edges and relations.

Notes:
- Missing points in experimental data are replaced with NA and if they are in the middle of the series, they are estimated by taking average of neighboring data values.
- AKT inhibitor (GSK690693) and MEK inhibitor (GSK1120212) are assumed to have the same effect on their phosphorylated forms. This is not considered for FGFR1/3 inhibitor (PD173074), its target is not found in the data.

The scripts require "simp.R" from StreamMetabolism R package (http://cran.r-project.org/web/packages/StreamMetabolism/index.html) for integration so make sure you have simp.R in the same directory with the scripts.

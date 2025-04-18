## Donahue_Lab_StandaloneTools
### Short, task-specific scripts  

A collection of small scripts that address individual tasks or problems. Scripts in this repository are standalone and not part of a larger pipeline.  
  
  

Scripts within this repository:  



<details>
<summary> KM_StratifiedPlot.R </summary>

**Description:** Generates Kaplan-Meier survival curves for a list of features stratified by a chosen percentile cutoff (e.g. Median, upper 25th, etc.). Useful for survival analysis comparing these two user-defined groups.



**The output .csv:**
| feature | p value | cutoff |
|-----:|-----:|-----:|


**Visualization:** Kaplan-Meier plot

</details>



<details>
<summary> Wilcoxon_Test.R </summary>

**Description:** Paired Wilcoxon Test (or unpaired) for multiple features between two groups. Useful for identifying significantly different distributions.

**The output .csv:**
| feature | p value | median1 | median2 |
|-----:|-----:|-----:|-----:|

**Visualization:** Box plot with p-values

</details>


<details>
<summary> TimeSeries_GridPlot.R </summary>

**Description:** Creates a 2x2 grid with group-based coloring and optional axis breaks. Useful for visualizing trends across treatment groups.

**Visualization:** Line graph (2x2 grid layout)

</details>




<details>
<summary> Heatmap_BubblePlot.R </summary>

**Description:** Generates heatmaps and bubble plots. Useful for exploring up/down regulation of pathways or other high-dimensional data. *Coming soon.*

**Visualizations:** Heatmap and bubble plots *Coming soon.*

</details>



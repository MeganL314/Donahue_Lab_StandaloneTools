## Donahue_Lab_StandaloneTools
### Short, task-specific scripts  

A collection of small scripts that address individual tasks or problems. Scripts in this repository are standalone and not part of a larger pipeline. Click on script name for more details.
<br>
<details>
<summary> KM_StratifiedPlot.R </summary>
<br>

**Description:** Generates Kaplan-Meier survival curves for a list of features stratified by a chosen percentile cutoff (e.g. Median, upper 25th, etc.). Useful for survival analysis comparing these two user-defined groups.  
  
**The output .csv:**
| feature | p value | cutoff |
|-----:|-----:|-----:|

**Visualization:** Kaplan-Meier plot

<br>
<br>
</details>


<details>
<summary> Wilcoxon_Test.R </summary>
<br>
  
**Description:** Paired Wilcoxon Test (or unpaired) for multiple features between two groups. Useful for identifying significantly different distributions.

**The output .csv:**
| feature | p value | median1 | median2 |
|-----:|-----:|-----:|-----:|

**Visualization:** Box plot with p-values

<br>
<br>
</details>


<details>
<summary> Wilcoxon_Paired_Complete_Cases.R </summary>
<br>
  
**Description:** Performs a paired Wilcoxon signed-rank test after removing any incomplete pairs (if either value in a pair is NA, the entire pair is excluded)

**The output .csv:**
| p value | feature | groups | sample size | median1 | median2 |
|-----:|-----:|-----:|-----:|-----:|-----:|

<br>
<br>
</details>






<details>
<summary> TimeSeries_GridPlot.R </summary>
<br>
  
**Description:** Creates a 2x2 grid with group-based coloring and optional axis breaks. Useful for visualizing trends across treatment groups.

**Visualization:** Line graph (2x2 grid layout)
<br>
<br>
</details>




<details>
<summary> Heatmap_BubblePlot.R </summary>
<br>
  
**Description:** Generates heatmaps and bubble plots. Useful for exploring up/down regulation of pathways or other high-dimensional data.

**Visualizations:** Heatmap and bubble plots
<br>
<br>
</details>




<details>
<summary> Correlation.R </summary>
<br>
  
**Description:** *Coming soon*
<br>
<br>
</details>




<details>
<summary> CommonDataWrangling.R </summary>
<br>
  
**Description:** Summary of common data wrangling solutions: cleaning, preprocessing, shaping, formatting, etc. Tidyverse and base R solutions.
<br>
<br>
</details>
<br>

Here you will find the repository for the Ritger et al 2020 paper, *"Diet choice in a generalist predator, the invasive lionfish (Pterois volitans/miles)"*. The repository is maintained by Amelia Ritger (GitHub: @ameliaritger) at the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.


**Abstract:** Diet choice in marine species is typically derived from indirect methods such as stomach contents and stable isotope analysis, while choice experiments in controlled laboratory settings are used to infer foraging decisions in the wild. However, these methods are limited in their capacity to make inferences about foraging decisions by predators in variable environments or recreate the array of factors (such as prey traits, predator condition, and environmental conditions) present in natural systems which may interact to affect diet decisions by predators. Recent work has provided evidence for selectivity in the invasive Indo-Pacific lionfish (Pterois volitans/miles) despite the predator’s apparent opportunistic, generalist feeding behavior. We directly tested diet choice by presenting wild-caught lionfish with multi-species prey assemblages in field enclosures. We offered lionfish equal biomasses of prey species sharing similar prey traits that are both highly abundant on coral reefs and prevalent in the lionfish diet across the invaded range. We then applied compositional analyses to determine relative prey consumption given prey availability. We observed lionfish selectively foraging on prey and manifesting strong consistent preferences for one prey species. Additionally, we observed condition-dependent foraging behavior, as lionfish with higher body conditions were more likely to exhibit selective foraging behavior. Our findings provide direct evidence for diet choice in an invasive generalist species and highlight the importance of preserving the ecological complexity of natural ecosystems in choice experiments, particularly when investigating predator-prey interactions in complex environments.


## Code
file name | file overview | description 
---|---|-----------
Lion_PrefIndex_cluster.R | Statistical analysis and Figures | compositional data analysis - application of ilr transformation and creation of index of selectivity; creation of Figure 2
Ternary_plots.R | Plotting ternary diagram | creation of Figure 3

## Metadata for "lionfishdata.csv"
*Data has been cleaned from /outdated-materials/All_data_Sept8 for analysis*

Each row represents one individual lionfish	

variable | description
---|---
Sample No. |	Lionfish number
Date captured	| Capture date of lionfish
Time captured	| Capture time of lionfish
Depth captured (m) | Depth where lionfish was captured, in meters
Quarantine No. |	Holding cage number where lionfish was placed between capture and trial
Date trial |	Date of the trial
Julian Date |	Julian date of the trial
Start time trial |	Start time of the trial
End time trial |	End time of the trial
Starvation time (hours) |	Number of hours between capture and trial start time when lionfish did not eat
Moon cycle |	Moon cycle on trial date
Cloud cover |	Cloud cover on trial date, averaged every 2 hours from trial start time until sunset
Current	| Ocean current at trial location during day of trial
Sunset time |	Sunset time on trial date
Enclosure No. |	Enclosure where lionfish was placed for trial
Sex |	Sex of lionfish
Lionfish total length (cm) |	Total length of lionfish in centimeters
Lionfish standard length (cm) |	Standard length of lionfish in centimeters
Lionfish Wet Weight (g) |	Wet weight of lionfish in grams
Body Condition | Total length  / wet weight ^3
Number chromis consumed |	Total number of brown chromis consumed
Proportion chromis consumed |	Number consumed / number provided
Proportion biomass chromis | Biomass* chromis consumed / total prey biomass consumed
Proportion diet chromis | Number chromis consumed / total number prey consumed
Number wrasse consumed | Total number of bluehead wrasse consumed
Proportion wrasse consumed | Number consumed / number provided
Proportion biomass wrasse | Biomass wrasse consumed / total prey biomass consumed
Proportion diet wrasse | Number wrasse consumed / total number prey consumed
Number goby consumed | Total number of glass gobies consumed
Proportion goby consumed | Number consumed / number provided
Proportion biomass goby | Biomass goby consumed / total prey biomass consumed
Proportion diet goby | Number goby consumed / total number prey consumed
Total fish consumed	| Total number fish consumed / total number fish provided

*Estimated biomass for each prey species based on established averages of weights and lengths

## Figures (extracted from Ritger et al, 2020)
![Alt text](/figures/fig1.png?raw=true)
**Figure 1** (A) Schematic diagram of enclosure dimensions and experimental design. (B) Holding cages were located proximal to enclosures and (C) maintained lionfish in isolation and without food prior to predation trials.

![Alt text](/figures/fig2.png?raw=true)
**Figure 2** Relationship between lionfish body condition and index of selectivity. Lionfish with higher body condition had a significantly higher mean index of selectivity than lionfish with lower body condition. Point shape indicates significantly different groups, separated by a body condition value of 0.0124. The dashed and dotted lines represent the mean index of selectivity of lionfish from each group.

![Alt text](/figures/fig3.png?raw=true)
**Figure 3** Ternary diagram showing index of selectivity values from predation trials. Most lionfish exhibited some or strong preference for brown chromis. Each point in the ternary diagram represents the index of selectivity of an individual lionfish. The size of each point indicates the number of lionfish sharing an index of selectivity. Disproportionate consumption of prey shifts a point away from the + symbol, representing prey consumption proportional to prey availability (in accordance with the 11:11:55 ratio of prey offered), towards a corner representing the targeted prey species. The ○ symbol represents the geometric mean of lionfish prey consumption across all predation trials. The dashed contour lines depict the index of selectivity at 0.1 increments.

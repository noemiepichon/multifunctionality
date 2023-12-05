MULTIFUNCTIONALITY DATASETS


Datasets used in "Nitrogen availability and plant functional composition modify biodiversity-multifunctionality relationships" 
from Noémie A. Pichon, Seraina L. Cappelli, Tosca Mannall, Thu Zar Nwe, Norbert Hölzel, Valentin H. Klaus, Till Kleinebecker, Santiago Soliveres, Hugo Vincent and Eric Allan
Initial manuscript can be found on bioRxiv: https://www.biorxiv.org/content/10.1101/2020.08.17.254086v1 

This experiment manipulates in a full factorial design species richness (1, 4, 8 species), nitrogen enrichment (0, 100 kg ha-1 year-1), fugicide spraying (unsprayed, sprayed) and initial community composition in SLA (gradient in growth strategies).

See methods for more info on the design. If you have further questions: noemie.pichon@wsl.ch

These datasets can be used with code on GitHub:
https://github.com/noemiepichon/multifunctionality

To calculate multifunctionality, please use the code in GitHub "Megamodel github.R"
To calculate species effect on functioning, please use the code in "Sp effect on function github.R"

Attention: the functions were transformed before calculating multifunctionality. See methods.




In the first dataset "All_BigData.txt", we put together the ten functions measured on each plot of the experiment, and each plot community metrics (SLA, richness etc). 

All_BigData.txt:

Plot_Nr
Block				1 to 4
Species_richness		1 to 20
Functional_composition		F, M, S, for Fast, mixed or slow initial sown composition
Nitrogen			0 or 1
Fungicide			0 or 1
Combination			unique combination of species, growing once in control conditions, once with nitrogen, once with fungicide, once with both
CWM_SLA				Community specific leaf area calculated with monoculture SLA values and plot-specific abundances
MPD_SLA_abundance		Mean pairwise distance in specific leaf area calculated with monoculture SLA values and plot-specific abundances
MPD_SLA_presence		Mean pairwise distance in specific leaf area calculated with monoculture SLA values equally weighted
Aboveground_biomass		(will be sqrt trandformed)
Herbivory			(will be sqrt trandformed)
Pathogens			(no transformation needed)
Soil_respiration		(will be log trandformed)
Plant_N_uptake			(will be log trandformed)
Plant_P_uptake			(will be log trandformed)
Belowground_biomass		(will be log trandformed)
BGlucosidase			(will be sqrt trandformed)
Phosphatase			(will be log trandformed)
Carbon_storage			(will be sqrt+1 trandformed)




The second dataset "Mean_SLA.txt" is the mean SLA per control monoculture, averaged over the years.

Plot_Nr				
CWM_SLA				Community specific leaf area calculated with monoculture SLA values and plot-specific abundances
Species				Abbreviation for species identity




The last dataset "species_presence.txt" indicates if a species is present (1) or absent (0) of a plot.

Am				Achillea millefolium
Ao				Anthoxanthum odoratum
As				Anthriscus sylvestris
Be				Bromus erectus
Cb				Crepis biennis
Cj				Centaurea jacea
Dc				Daucus carotta
Dg				Dactylis glomerata
Fr				Festuca rubra
Ga				Galium album
Hl				Holcus lanatus
Hp				Helictotrichon pubescens
Hs				Heracleum sphondylium
Lp				Lolium perenne
Pg				Prunella grandiflora
Pm				Plantago media
Pt				Poa trivialis
Ra				Rumex acetosa
Sp				Salvia pratensis
To				Taraxacum officinale
Block				1 to 4
Plot_Nr




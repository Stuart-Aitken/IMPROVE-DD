# IMPROVE
<table style="width:100%">
  <tr> <th style="width:30%">IMPROVE-DD:<br>
     Integrating Multiple Phenotype <br>Resources Optimises Variant Evaluation in genetically determined Developmental Disorders<br>
    Aitken, FitzPatrick, Semple et al. 2022 <br> <a href="https://doi.org/10.1101/2022.05.20.22275135">https://doi.org/10.1101/2022.05.20.22275135<a></th>
    <th style="width:70%"> 
  <img src="https://github.com/Stuart-Aitken/IMPROVE-DD/blob/main/improve.png" width="300" alt="IMPROVE-DD">
     </th>
  </tr>
  </table>

Code used for Figures 2 and 3. Depends on the Human phenotype ontology: https://hpo.jax.org/
Version 02/08/2021 available: https://bioportal.bioontology.org/ontologies/HP


To run either the TF IDF IC analysis or the IMPROVE HPO term classifier, a matrix of cases * HPO terms needs to be created (see IMPROVE_resource.r). Should a database of HPO annotations not be available, a small set of demo data can be created using IMPROVE_resource.r to enable the main code to be run.


All R code:
#Copyright (C) 2022 The University of Edinburgh 
#Author Stuart Aitken MRC HGU IGC s.aitken@ed.ac.uk
#All Rights Reserved.
#Funded by the Medical Research Council
#https://www.ed.ac.uk/mrc-human-genetics-unit**


License: GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

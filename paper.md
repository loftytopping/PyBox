---
title: ' PyBox: An automated box-model generator for atmospheric chemistry and aerosol simulations.'
tags:
  - Python
  - Atmosphere
  - Air-quality
  - Climate
  - Chemistry
  - Aerosol 
  - Fortran
authors:
  - name: David Topping
    orcid: 0000-0001-8247-9649
    affiliation: 1 
  - name: Paul Connolly
    orcid: 0000-0002-3294-7405
    affiliation: 1
  - name: Jonathan Reid
    orcid: 0000-0001-6022-1778
    affiliation: 2
affiliations:
 - name: School of Earth and Environmental Science, University of Manchester, M13 9PL
   index: 1
 - name: School of Chemistry, University of Bristol, BS8 1TS
   index: 2

date: 18 May 2018
bibliography: paper.bib
---

# Summary

Air pollution and climate change are two of the biggest multidisciplinary challenges in society today.  The need to understand the chemical and physical processes in the atmosphere that dictate the impacts of both has created a wide range of research platforms.  These include numerical models that aim to capture the aforementioned processes. Volatile organic compounds (VOCs), emitted from both natural and anthropogenic sources, are oxidised in the atmosphere to form lower-volatility species that form organic particulate matter in the atmosphere through gas-to-particle partitioning. Chemical mechanisms, such as the Master Chemical Mechanism (MCM) [@Jenkin2003], have been built to hold information of all relevant species and reactions in the atmosphere.  

However, it is becoming extremely hard to develop models that can respond to our growing body of knowledge since these mechanisms now treat thousands of species and tens of thousands of reactions. Researchers need to simulate the evolution of each individual species over a range of time-scales and ambient conditions. Manually setting up the relevant ordinary differential equations definitions [ODEs] and associated solvers would present a huge challenge given the number of compounds and equations involved. In addition, it is important to test sensitivity to rate coefficients or evaluate the impact of emerging laboratory data that might identify new reaction pathways or new, improved rate coefficients. Added to this complexity is the need to account for gas-to-particle partitioning of each compound. The atmosphere has varying concentrations of particulate matter that might act as a condensational sink for each compound. Predicting that partitioning requires calculations of molecular properties including saturation vapour pressures. Automation is essential, not least to ensure the gas phase chemistry is accounted for alongside the gas-to-particle partitioning and reproducibility is ensured. This is the driver behind ``PyBox``. 

``PyBox`` is a 0-D box model, where all species are homogeneously distributed, built around a chemical mechanism file; a file that represents all the individual chemical reactions taking place starting from a wide range of precursors. By parsing a chemical equation file, obtained from the MCM project, ``PyBox`` creates files that account for the gas phase chemistry as well as automatically calculating properties that dictate gas-to-particle partitioning through connection with the ``UManSysProp`` informatics suite [@Topping2016]. Written in Python, ``PyBox`` uses ``Numba`` [@Lam2015] or the Fortran-to-Python-Generator ``f2py`` [@Peterson2009] to perform calculations currently within a library of ODE solvers provided by the ``Assimulo`` package [@Andersson2015]. 

# Acknowledgements

This project was partly supported through NERC grant NE/N013794/1 which funded me to spend time developing new community models for sensitivity studies. 

# References

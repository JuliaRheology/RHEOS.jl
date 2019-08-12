---
title: 'RHEOS -- A Julia package for Rheology data analysis.'
tags:
  - Julia
  - Rheology
  - Fractional Rheology
  - Viscoelasticity
  - Fractional Viscoelasticity
  - Biomechanics
authors:
  - name: Jonathan Louis Kaplan
    affiliation: 1
  - name: Alessandra Bonfanti
    affiliation: 1
  - name: Alexandre J Kabla
    orcid: 0000-0002-0280-3531
    affiliation: 1
affiliations:
  - name: Engineering Department, Cambridge University, UK
    index: 1
date: 12 August 2019
bibliography: main.bib
---
# Summary
Rheology is generally understood to refer to the science of deformation and flow. It encompasses a broad range of experimental methodologies (such as macro-scale metal tensile testing and micro-scale indentation tests of biological samples) and a correspondingly broad range of theoretical tools. 

RHEOS (Rheology, Open-Source) is a software package that is designed to make the analysis of rheological data simpler, faster and more reproducible. RHEOS has a particular emphasis on rheological models containing fractional derivatives which have demonstrable utility for the modelling of biological materials [@bonfantiUnifiedRheologicalModel2019; kaplanPectinMethylesterificationImplications2019] but have hitherto remained in relative obscurity -- possibly due to their mathematical and computational complexity. RHEOS is written in Julia [@bezansonJuliaFreshApproach2017], which greatly assists achievement of our aims as it provides excellent computational efficiency and approachable syntax. RHEOS is fully documented and has extensive testing coverage.

RHEOS has already been used in one study [@bonfantiUnifiedRheologicalModel2019] and is currently being used in at least one other study and review that are in progress.

# Statement of Need


Indeed, possibly due to the computational simplicity of the above, many studies do simply assume that their loading is accurately described by a step function. However the validity of this approach is not always tenable.

Interdisciplinary nature

As mentioned above, the utility of fractional viscoelastic models has been well documented [@bonfantiUnifiedRheologicalModel2019]. However, attaining a well-rounded understanding of fractional viscoelastic theory can be difficult: of popular rheology textbooks published in 1989 [@barnesIntroductionRheology1989; @findleyCreepRelaxationNonlinear1989], 2008 [@brinsonPolymerEngineeringScience2008], 2009 [@lakesViscoelasticMaterials2009] and 2013 [@wardMechanicalPropertiesSolid2013], fractional viscoelasticity is only mentioned briefly in one of them [@lakesViscoelasticMaterials2009].

Furthermore, 

# Functionality and Implementation

# 

# Acknowledgements
JLK Would like thank the George and Lillian Schiff Foundation for the PhD funding which facilitated this project.

# References
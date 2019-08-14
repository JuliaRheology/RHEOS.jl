---
title: 'RHEOS.jl -- A Julia package for Rheology data analysis.'
tags:
  - Rheology
  - Fractional Rheology
  - Viscoelasticity
  - Fractional Viscoelasticity
  - Biomechanics
  - Julia
authors:
  - name: Jonathan Louis Kaplan
    orcid: 0000-0002-2700-5229
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
bibliography: paper.bib
---
# Summary
Rheology is generally understood to refer to the science of deformation and flow. It encompasses a broad range of experimental methodologies (such as macro-scale metal tensile testing and micro-scale indentation tests of biological samples) and a correspondingly broad range of theoretical tools.

RHEOS (Rheology, Open-Source) is a software package that is designed to make the analysis of rheological data simpler, faster and more reproducible. RHEOS has a particular emphasis on rheological models containing fractional derivatives which have demonstrable utility for the modelling of biological materials [@bonfantiUnifiedRheologicalModel2019; @kaplanPectinMethylesterificationImplications2019] but have hitherto remained in relative obscurity -- possibly due to their mathematical and computational complexity. RHEOS is written in Julia [@bezansonJuliaFreshApproach2017], which greatly assists achievement of our aims as it provides excellent computational efficiency and approachable syntax. RHEOS is fully documented and has extensive testing coverage.

RHEOS has already been used in one study [@bonfantiUnifiedRheologicalModel2019] and is currently being used in at least one other study and review that are in progress.

# Statement of Need
When fitting regular viscoelastic models to data under the assumption of step loading, the process is straightforward as it reduces to a direct optimisation of a function which consists of a small number of exponential terms. Indeed, possibly due to the computational convenience of the above, many studies do assume that their loading is accurately described by a step function and a traditional viscoelastic model. However the validity of this approach is not always tenable. With regard to the step assumption, a number of studies have noted the importance of small time-scale behaviour when fitting viscoelastic models [@bonfantiUnifiedRheologicalModel2019; @dipaolaInfluenceInitialRamp2014; @oyenSphericalIndentationCreep2005]. Taking into account a complete arbitrary loading history requires computing the viscoelastic hereditary integral:

$$ \sigma(t) = \int_{0}^t G(t - \tau) \frac{d \epsilon(\tau)}{d \tau} d \tau $$

Depending on the model kernel *G* and the sample rate(s) of the data, computation is not always straightforward and there are several approaches that can be used. With regard to the choice of model, the utility of fractional viscoelastic models has been well documented [@bonfantiUnifiedRheologicalModel2019]. Furthermore, these fractional viscoelastic models are often valuable in modelling biological materials, which may be studied by highly interdisciplinary research groups where there is neither the time or expertise required to utilise fractional viscoelastic modelling.

RHEOS reduces the two biggest complexities presented above, arbitrary loading history and fractional model kernels, to a simple, approachable syntax. It includes all of the most common traditional and fractional viscoelastic models and additional models can easily be added by users.

Obtaining intuition for fractional viscoelastic theory can be difficult, and learning material is sparse: of popular rheology textbooks published in 1989 [@barnesIntroductionRheology1989; @findleyCreepRelaxationNonlinear1989], 2008 [@brinsonPolymerEngineeringScience2008], 2009 [@lakesViscoelasticMaterials2009] and 2013 [@wardMechanicalPropertiesSolid2013], fractional viscoelasticity is only mentioned briefly in one of them [@lakesViscoelasticMaterials2009]. Given the above, RHEOS also includes signal generation features so that common loading patterns (e.g. step, ramp, stairs) can be applied to unfamiliar models to gain great qualitative intuition of their behaviour.

# Implementation Demonstration
In this section, the architecture of RHEOS is briefly discussed, accompanied by a demonstrative example. It should be noted that RHEOS is not, at its heart, an optimisation package. It builds on a another optimisation package, NLopt [@johnsonNLoptNonlinearoptimizationPackage] by adding a large number of abstractions specific to the exploration of viscoelastic data. Some of these abstractions are featured in the example below in which experimental viscoelastic data is fitted to a fractional viscoelastic model (spring-pot) and then this model is used to make a prediction of the behaviour so that its qualitative accuracy can be assessed.

 The RHEOS-specific commands used to load in the data, fit a model and predict using that model are the following:
- `data = importcsv("example1_data.csv", t_col=1, ϵ_col=2, σ_col=3)`
- `maxwell_model = modelfit(data, Maxwell, strain_imposed)`
- `maxwell_predict = modelpredict(extract(data, strain_only), maxwell_model)`

In brief: the `importcsv` command is a convenience function for loading in CSV data into a `RheoTimeData` struct; `modelfit` takes the data and a `RheoModelClass` which it should be fitted to; `modelpredict` takes the fitted model, and a partial form of the data and uses these arguments to predict the response according to the fitted `maxwell_model`. For a more detailed discussion of the syntax used, see the RHEOS documentation. This common workflow is shown schematically in Figure 1, along with the plotted data in Figure 2.

![RHEOS Workflow Diagram](diagram_resize_v2.png)

![Resultant Figure](predictfigure.png)

Note that although the above example dealt with time data, RHEOS handles frequency-based viscoelastic data equally well and models fitted on one type of viscoelastic moduli can easily be used to predict behaviour with a different moduli. (For instance, fitting against relaxation data and predicting the frequency response spectrum.)

# Acknowledgements
JLK Would like thank the George and Lillian Schiff Foundation for the PhD funding which facilitated this project.

# References

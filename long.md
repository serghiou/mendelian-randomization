# A method of data analysis you probably do not know but probably should: Mendelian Randomization

<div align="justify">

## TL;DR

* Methods of causal inference attempt to make causal statements using observational data
* Instrumental variables can be used for causal inference when not all confounders are known
* Mendelian randomization (MR) is an instrumental variable method using genetic instruments
* Advancements in MR provide an opportunity to probe previously impossible questions
* Estimates from MR analyses should be viewed with caution as they are hard to interpret

## Introduction

Observational data, unlike experimental data, are abundant, yet difficult to use in making causal statements. Mendelian randomization can harness what we know about the human genome to inform our understanding of the association between a putative cause and effect using observational data. It was initially proposed in the late 1980s (Katan, 1986; Katan, 2004), used in a handful of occasions thru the 1990s (Gray and Wheatley, 1991) and eventually formally established in the early 2000s (Smith and Ebrahim, 2003). Since then, it has gained substantial popularity and has been used to make several claims that were upheld by randomized controlled trials (RCTs) (Bennett and Holmes, 2017). However, very few data scientists are in fact aware of how to use this method or evaluate its results. I hereby aim to offer a basic introduction.

## Causal inference

Before we examine MR, we need to define causality. A causal factor is one in which a modification across at least one of its levels would modify the outcome, all else held constant. For example, if we take smoking as the exposure and cancer as the outcome, then smoking increases the probability of developing cancer. This is not the same as carrying a lighter, which albeit associated to developing cancer, it does not cause cancer. We care about causality because it implies that we can use this knowledge to intervene (Glass et al., 2013). For example, we know that an intervention to help someone give up smoking will reduce their probability of developing cancer, but an intervention to make someone not carry a lighter will not change the probability of developing cancer, unless through its impact on smoking. Such statements become more meaningful within diagrammatic acyclic graphs (DAGs).

The DAG associated to the example above is illustrated below. In this DAG, **X** refers to the exposure (smoking), **Y** refers to the outcome (cancer) and **W** refers to carrying a lighter. These are known as **nodes** and the directed arrows between nodes are the **edges**.  Altering X will have an effect on all nodes within its path - in this case, W and Y. Altering W will not have an effect on either X or Y because no edge from W feeds into any other node. 

(img)

_Tip. There is a lot more to learn about causal inference and [[this]] freely available book by Harvard professors Miguel Hernán and Jamie Robins is an excellent place to start._

## Instrumental Variables (IV)

To understand MR we need to first understand Instrumental Variables (IV) as MR is in fact a type of IV analysis. IV estimators were developed early on in the econometric literature (Wright, 1928; Theil, 1958) and have since been applied to measure many estimates of causal effects (e.g. estimation of effect of maternal smoking on birth weight, effect of children’s age at school entry on their eventual educational attainment, etc.). The popularity of such methods stems from their ability to circumvent the cardinal sin of observational data - confounding. 

Confounding is a form of bias that arises because the treatment and outcome share a common cause - even though traditional analyses try to use our knowledge of covariates to adjust for such confounding variables, we very rarely observe all confounding variables we need to adjust for. IV analyses use a clever trick to estimate the effect of interest without relying on potentially observed or unobserved confounders. 

### What do we mean by an IV?

Consider the following DAG, in which **X** is the exposure and **Y** is the outcome, as above. In this case, X is a cause of Y.

(img)

Now consider the following DAG, illustrating observed **O** and unobserved **U** confounders. As noted above, a confounder is a shared cause between exposure X and outcome Y and is the main difficulty in making causal statements with observational data. For example, take ingestion of red meat as a putative risk factor for cancer.  Then, consider that smokers tend to eat more red meat and that smokers have a high probability of cancer. In this scenario, it is very difficult to know whether red meat is a cause of cancer as any effect of red meat observed on cancer may in fact have been mediated by smoking. 

(img)

The following two DAGs illustrate the IV **Z**.  In the DAG on the left, Z is causally associated to X, but in the DAG on the right it is not - IVs do not have to be causally associated to exposure X, but if they are, it improves our ability to interpret our final estimate.

(img) (img)

For Z to be an IV, it has to abide by the following 3 assumptions:

**Assumption 1: Relevance.** There is an association (causal or non-causal) between Z and X. The stronger this association between Z and X, the better.

**Assumption 2: Exchangeability.** There is no line between Z and O or U. In other words, Z is not causally related to the confounders of the association between X and Y.

(img)

**Assumption 3: Exclusion restriction.** There is no line between Z and Y. In other words, Z does not cause Y. 

(img)

Even though the above assumptions are enough to call Z an IV, we need a further assumption to be able to estimate the effect of X on Y, known as the Complier Average Causal Effect (CACE).

**Assumption 4: Monotonicity.** This implies that there are no defiers. In other words, for a binary Z and X, no individual exposed under no treatment, i.e. X = 1 under Z = 0, would be unexposed under treatment, i.e. X = 0 under Z = 1.  For example, no person would quit smoking if smoking is cheap and not quit smoking if smoking is expensive. 

Even though we can use the monotonicity assumption to estimate an effect, we need another assumption to be able to interpret our estimate of effect as the **Average Causal Effect (ACE)** - this is usually the effect of interest and the one calculated by randomized controlled trials (RCTs).

**Assumption 5: Homogeneity.** This implies that the average effect is equal within all of compliers (i.e. those that comply with treatment), always takers (those that always receive the treatment whether they are assigned to treatment or not) and never takers (those that never receive the treatment, whether they have are assigned to treatment or not).

### How do we estimate the effect of interest using an IV?

Let us define,  ~$\beta = \text{effect}$~. Then, by the DAG above we can see that,

$$\beta_{zy} = \beta_{zx} \beta_{xy} \quad \Rightarrow \quad \beta_{xy} = \frac{\beta_{zy}}{\beta_{zx}}$$

This implies that, for a binary exposure and IV, the estimator of effect is,

$$\frac{\text{E}\{Y \mid Z = 1 \} - \text{E}\{Y \mid Z = 0 \}}{\text{E}\{X = 1 \mid Z = 1 \} - \text{E}\{X = 0 \mid Z = 0 \}}$$

The above estimator is known as the Wald estimator (or IV estimator). However, the most commonly used estimator is the **Two-Stage Least Squares (TSLS)** estimator. This procedure calculates the estimate of interest in two stages:

**Stage 1.** Fit X on Z and possible covariates C by ordinary least squares (OLS) and compute the predicted value ~$\hat{X}$~,

$$\hat{X} = \hat{\text{E}}\{A \mid Z, C \} = \hat{\alpha}_0 + \hat{\alpha}_1Z + \hat{\alpha}_2^TC$$

** Stage 2.** Fit Y on ~$\hat{X}$~,

$$\text{E}\{Y \mid \hat{X}, C \} = \mu_0 + \mu_1 \hat{X} + \mu_2^T C$$

Under assumptions 1-4, the coefficient of ~$\hat{X}$~ converges in probability to the causal effect of interest and with assumption 5, this is the desired ACE.

### How do we do this in R?

In R, it is possible to run an IV analysis using several packages. The most commonly used ones are `AER`, `ivmodel` and `ivpack`.  In the example below, I am using a well-known dataset to identify the causal effect of education on log-wage earned. To do so, I am using father education as an IV - this is in fact not the best IV because it is easy to consider how it may violate the exchangeability assumption. This code was taken from the STATS 266 course at Stanford University.

```
# Import packages
library(AER)       # ivreg 
library(ivpack)    # robust.se and anderson.rubin.ci

# Import data
dat <- read.table("http://statweb.stanford.edu/~rag/stat209/Mroz87.dat"
                  , header = T)

# Fit model use father education as IV
iv.fit1 <- ivreg(logWage ~ educ | fatheduc, dat)

# Observe model
summary(iv.fit1, diagnostics = T)

# Calculate robust standard error
robust.se(iv.fit1)

# Calculate the Anderson-Rubin confidence interval (robust to weak instruments)
anderson.rubin.ci(iv.fit1)
```

*Tip. For a great tutorial on IV methods for causal inference, have a look at [this] paper by Baiocchi et al. (2014).*

## Mendelian Randomization

MR is an IV method in which the IV of choice is genetic, most commonly what is known as a single nucleotide polymorphism (SNP). The correlation of SNPs to variables most commonly arises from Genome-Wide Association Studies (GWAS) and each SNP used in MR tends to be taken from a different gene. The term "Mendelian randomization" arises from  the "random assortment of alleles during meiosis where DNA is transferred from parent to offspring at the time of gamete formation, a process named Mendel’s second law" (Bennett and Holmes, 2017).

### What are the benefits of MR?

1. **Randomization.** As noted above, segregation and independent assortment during gametogenesis guarantees randomness of gene allocation.

2. **Reverse causality.** The genetic sequence is fixed from conception, which implies that the outcome could not have possibly caused the genetic instrument.

3. **GWAS.** Hundreds of genome-wide association studies (GWAS) have identified thousands of SNPs associated to many exposures. As such, there is an abundance of potential genetic instruments.

4. **Flexibility.** All of the already known associations between SNPs and exposures can be used to probe questions that were previously too expensive or unethical to probe. For example, it is possible to ask questions such as "does education increase myopia", without randomizing children to education or no education, probe treatment efficacy before running an RCT, probe target-mediated adverse effects and identify opportunities for drug repurposing (Bennett and Holmes, 2017).

5. **Automated causal inference.** Many are considering how recent developments in MR methods can make use of already known SNPs from large GWAS to probe questions in an automated or semi-automated fashion.

### What are the limitations of MR?

1. **Weak instruments.** This is a violation of the relevance assumption. It describes a situation in which the IV is only weakly associated with exposure. Weak instruments: (1) bias the estimate of effect towards the unadjusted estimate (i.e. the estimate we would have obtained had we not been using MR) (Burgess et al., 2011), (2) lead to a TSLS estimator with a non-Normal sampling distribution  regardless of sample size (Stock et al., 2002) and (3) reduce statistical power (Pierce et al., 2011). As a rule of thumb, when using a single IV with TSLS, seek an F-statistic greater than 10; loosely, the F-statistic arises from the null hypothesis test of no association between IV and exposure (refer to Stock et al. (2002) for a more accurate definition). This is not a panacea and many considerations have to be taken into account when dealing with weak instruments, including using estimators other than TSLS (Burgess et al., 2011; Stock et al., 2002). A relatively simple fix is to estimate the first stage regression into a separate sample from that used in the second stage of TSLS (the downside of this is loss in power). 

2. **Horizontal pleiotropy.** This is a violation of the exchangeability and restriction exclusion assumptions. Pleiotropy refers to a setting in which a genetic locus, such as a SNP, is related to more than one traits. This is a problem when these traits are either shared between exposure and outcome or between IV and outcome. Such pleiotropy can introduce substantial bias in the estimate of effect of X on Y. 

3. **Confounding by ancestry (= spurious pleiotropy).** This is a violation of the exclusion restriction assumption. It occurs when Z and Y share a common cause because of shared genetic ancestry (S). When this is the case, the observed effect of Z on Y does not only reflect its association through X, but also its association through genetic variants G (Pingault et al., 2018).

	(img)

### What can we do to mitigate the aforementioned limitations?

1. **Use multiple instruments.** This aims to reduce the weak instrument bias by using more than one IV. Most common practice is to create polygenic or allelic scores by weighted or unweighted linear combinations of more than one SNPs. For example, in a recent paper (Mountjoy et al., 2018), allele score for n alleles was calculated by,

	$$\text{Allele score} = \sum_{i=1}^n \text{weight}_i \times \text{dosage}_i
 
	The problem with introducing multiple instruments is that they exacerbate the horizontal pleiotropy problem, for which reason these methods are used in conjunction with the methods mentioned below. Another problem is that they reduce estimate interpretability (discussed later on).

2. **Use meta-analytic methods.** This aims to reduce the pleiotropy problem and offer a straight forward method to combine multiple instruments. The two most commonly used such methods are MR-Egger regression and inverse variance weighting (IVW). MR-Egger fits an ordinary least squares line through the IVs, where the X-axis refers to the strength of the Z-X association and the Y-axis refers to the strength of the Z-Y association. MR-Egger is robust to multiple invalid IVs (i.e. IVs that violate assumptions 1-3) and to directional pleiotropy (i.e. a pleiotropic effect that has a non-zero cumulative effect on the X-Y association). IVW weighs the effect identified through each IV by the proportion of the variance of its association with X to the total amount of such variance (after adding all variances together). This method has more statistical power than MR-Egger, but is not robust against directional pleiotropy, unless used within a random-effects model. A recently published paper offers a great overview of these methods (Pingault et al., 2018). The picture below was taken from this paper and illustrates IVs as blue points on the graph - the methods described here refer to the solid blue line (MR-Egger) and the solid red line (IVW). These methods can be fit using the `MendelianRandomization` package in R, an excellent vignette for which can be found [here].

	(img)

*Tip. Be wary of very strong instruments.*

### How are MR studies different from RCTs?

MR studies are oft compared to RCTs and called "nature’s RCTs". Indeed, many are looking forward to using MR studies to replace or reduce our reliance on RCTs (Gidding et al., 2012; Urbina et al., 2016). Even though in principle MR studies in which all aforementioned assumptions apply could provide accurate estimates of the true causal effect. This is rarely the case. An excellent review paper by Swanson et al. (2017) explores a few reasons why at the moment such comparisons are premature.

1. **Non-guaranteed exchangeability.** As noted above, assignment of genetic instruments can be confounded by ancestry. This is not the case in well done RCTs, where allocation of treatment assignment is independent of covariates or outcome. As such, even though aspects of population, such as ethnicity, can be used to mitigate the impact of such confounding, absence of residual confounding cannot be guaranteed. Furthermore, such "misclassifications" of treatment will reduce statistical power to observe an effect. 

2. **Unclear time zero.** In most MR studies, treatment assignment, eligibility criteria and ascertainment of outcome do not all begin at the same time, unlike RCTs, in which they do. As such, effects identified in MR studies may be biased by post-randomization events, such as in the case of genetic variants associated to being alive or eligible at the time events are recorded.

3. **Proxy randomizers.** In most MR studies, the association of genetic variant to exposure is not causal. This is unlike RCTs, in which all participants receive the treatment to which they were randomized. As such, (1) it is very difficult to know which subgroup of the population is characterized by the identified causal effect (Swanson and Hernán, 2017) and (2) a sufficiently large sample size is required to mitigate loss in power due to non-compliance (Hernán and Robins, 2006). To understand (1), consider that the effect identified in MR studies is "local"; this implies that it only characterizes a subgroup of "compliers" in which the monotonicity assumption applies. Unfortunately, this subgroup is unknown to us, but when the IV is causally associated with exposure, it is possible to characterize this subgroup and produce statements such as: "this result applies to 30% of the population" or "this result applies to older individuals". When it is not causally associated, this is much harder and often impossible.  This problem is further compounded when multiple genetic variants are used as IVs, where each of which may reflect a different subgroup of the population.

4. **Undefined adherence.** The definition of "adherence" to an exposure is often unclear in MR studies. This is not the case in RCTs, where a pre-registered protocol often makes such definitions explicit. Even though this does not affect our ability to estimate the effect of the genetic variant (i.e. the intention-to-treat analog), it limits its utility and interpretability as we cannot estimate how far this estimate may be from that of the treatment strategy. Measuring adherence becomes even more complicated when using sustained treatment strategies (rather than once-off treatments), where changes in adherence would have to be measured over time to provide what is known as a per protocol estimate (i.e. an estimate of effect had all participants been treated as per protocol).

5. **No treatment pathway.** As noted above, protocols are very rarely defined for MR studies, which means that there is no clear definition of the treatment strategies of interest. As such, it is difficult to know which one of these effects is the one found by an MR study. For example, as noted in the original paper, those who inherit some high-risk variant used in MR may be "assigned at conception to a treatment strategy of continuous obesity throughout the life-course", or "to a treatment strategy of progressive weight gain after a childhood of normal weight", etc. "Each of these treatment strategies will result in different effect estimates and it is unclear which one of them corresponds to the one targeted by the MR study. Because of this, it is all the more difficult to describe the resulting bias from an MR study using a classical instrumental variable analysis to study a sustained treatment strategy: without first addressing the ambiguity in the per-protocol effect of interest, we cannot logically address whether our estimates are valid or approximately valid."

6. **Time-varying treatments.** Even though most MR studies are of time-varying treatments, the techniques I am thus aware of are limited to settings for baseline treatments only.

7. **Per protocol.** "MR studies, are often not focused on the effect of the genetic variants themselves because the goal is not generally to inform decision-making that involves direct genetic manipulation. Instead, MR studies use genetic variation to study other treatment effects that are the analogue of per-protocol effects in randomized trials." "By precisely defining the causal effects being estimated, we underscore that MR estimates are often vaguely analogous to perprotocol effects in randomized trials, and that current MR methods for estimating analogues of per-protocol effects could be biased in practice."

### What is best practice in running a MR study?

First, a good MR study needs to specifically define its protocol, the main features of which should be (Swanson et al., 2017): eligibility criteria, treatment strategies under study, assignment procedures, follow-up period, outcome, causal contrasts of interest, and analysis plan. Second, all MR studies need to use sensitivity analyses to confirm the extent to which their results are robust to deviations from their assumptions. Even though this is usually done for the relevance and exclusion restriction assumptions, it is rarely done for the monotonicity assumption. However, the monotonicity assumption is a strong assumption that has to be justified in a case-by-case basis and explored in such a way as to make the estimate of effect more interpretable (Swanson and Hernán, 2017). Useful guidance in evaluating IV assumptions can be found in Glymour et al. (2012), even though these tests are not definitive and should supplement rather than replace scientific judgement (Burgess, 2012). 


### Should I believe analyses done by MR?

Theoretically, under the assumptions described above, MR should be able to identify the true causal effect. Even though the nature of genetic instruments comes with substantial limitations (weak instruments, horizontal pleiotropy), recent advancements in methods of producing MR studies has substantially reduced the impact of those limitations. Furthermore, so far, analyses done by MR seem to produce estimates that are fairly similar to those seen in RCTs (Voight et al., 2012; O’Donoghue et al., 2014; etc.). 

However, most MR studies produce estimates that are very hard to interpret. This is because those estimates represent the local per protocol effect in the unknown subgroup of “complier” participants "who would have complied with the assigned treatment level had they been assigned to either randomization arm", i.e. they only represent effect of treatment on individuals whose treatment status can be manipulated by the IV (Imbens and Angrist, 1994; Angrist et al., 1996), but it is very hard to know who these individuals are and when using multiple IVs this becomes very complicated.  Furthermore, changes in study participation after randomization can introduce selection bias that is not accounted for in most MR studies.

As such, estimates provided by MR studies done well and in a manner attempting to emulate an imaginary RCT can provide a substantial improvement over studies of observational data not using any methods found in causal inference. However, its estimates should be interpreted with caution and as suggestive of the possible causal effect, rather than a definitive account. Claims that MR studies can substitute RCTs are currently premature and unsubstantiated.

## What if not MR?

There are many methods of causal inference for observational studies. The most widely known ones are inverse propensity score weighting, matching, g-computation, difference in differences, and regression discontinuity designs. A fantastic introduction on most of these topics can be found in the Imbens and Rubin book [Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction].

## Acknowledgements

In preparing this document I have used knowledge and slides from the Causal Inference workshop at Harvard. I will keep updating this article on the basis of comments received and how my own knowledge and opinions change over time. I would also like to thank Constantinos Parisinos for introducing me to the topic, providing important bibliography and giving me the opportunity to work with him on an MR study of our own (Parisinos et al., 2018).

## Additional resources

 [Video by George Davey Smith.]

## References

Angrist JD, Imbens GW, Rubin DB. Identification of Causal Effects Using Instrumental Variables. JASA. 1996 Jun;91(434):444-55

Baiocchi M, Cheng J, Small DS. Instrumental variable methods for causal inference. Stat Med. 2014;33:2297–2340.

Burgess S, Thompson SG; CRP CHD Genetics Collaboration. Avoiding bias from weak instruments in Mendelian randomization studies. Int J Epidemiol. 2011 Jun;40(3):755-64.

Burgess S. Re:“Credible Mendelian Randomization Studies: Approaches For Evaluating The Instrumental Variable Assumptions”. Am J Epidemiol. 2012 Sept 01;176(5):456-7

Gidding SS, Daniels SR, Kavey RE; Expert Panel on Cardiovascular Health and Risk Reduction in Youth. Developing the 2011 Integrated Pediatric Guidelines for Cardiovascular Risk Reduction. Pediatrics. 2012;129:e1311–e1319. 

Glymour MM, Tchetgen Tchetgen EJ, Robins JM. Credible Mendelian randomization studies: approaches for evaluating the instrumental variable assumptions. Am J Epidemiol. 2012 Feb 15;175(4):332-9. 

Gray R, Wheatley K. How to avoid bias when comparing bone marrow transplantation with chemotherapy. Bone Marrow Transplant  1991;7 (Suppl.3):9–12.

Hernán MA, Robins JM. Instruments for causal inference: an epidemiologist’s dream? Epidemiology. 2006;17:360–372.

Hernán M., Robins J. Causal Inference Book. 2018 Edition. Access at: [https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/]

Imbens GW, Angrist JD. Identification and Estimation of Local Average Treatment Effects. Econometrica. 1994 Mar;62(1):467-75

Katan M, isoforms AE. Serum cholesterol, and cancer. Lancet 1986;1:507–8.

Katan MB1. Commentary: Mendelian Randomization, 18 years on. Int J Epidemiol. 2004 Feb;33(1):10-1.

O’Donoghue ML, Braunwald E, White HD, et al. Effect of darapladib on major coronary events after an acute coronary syndrome: the solid-TIMI 52 randomized clinical trial. JAMA 2014;312:1006–15.

Parisinos CA, Serghiou S, Katsoulis M, George MJ, Patel RS, Hemingway H, Hingorani AD. Variation in Interleukin 6 Receptor Gene Associates with Risk of Crohn's Disease and Ulcerative Colitis. Gastroenterology. 2018 May 15.

Pierce BL, Ahsan H, Vanderweele TJ. Power and instrument strength requirements for Mendelian randomization studies using multiple genetic variants. Int J Epidemiol. 2011 Jun;40(3):740-52.

Pingault JB, O'Reilly PF, Schoeler T, Ploubidis GB, Rijsdijk F, Dudbridge F. Using genetic data to strengthen causal inference in observational research. Nat Rev Genet. 2018 Jun 5.

Smith GD, Ebrahim S. 'Mendelian randomization': can genetic epidemiology contribute to understanding environmental determinants of disease? Int J Epidemiol. 2003 Feb;32(1):1-22.

Stock J, Yogo M, Wright J. A Survey of Weak Instruments and Weak Identification in Generalized Method of Moments. J Bus Econ Stat. 2002;20:518–529.

Swanson SA, Hernán MA. The challenging interpretation of instrumental variable estimates under monotonicity. Int J Epidemiol. 2017 Mar 30. 

Swanson SA, Tiemeier H, Ikram MA, Hernán MA. Nature as a Trialist?: Deconstructing the Analogy Between Mendelian Randomization and Randomized Trials. Epidemiology. 2017 Sep;28(5):653-659.

Urbina EM, de Ferranti SD. Lipid screening in children and adolescents. Lipid screening in children and adolescents. JAMA. 2016;316:589–591.

Voight BF, Peloso GM, Orho-Melander M, et al. Plasma HDL cholesterol and risk of myocardial infarction: a mendelian randomisation study. Lancet. 2012 Aug 11;380(9841):572-80. 


</div>





# Protocol v1 (FROZEN 14 October 2026, before the confirmatory run)

Title: Hypothesis Mismatch and Variance Mismatch in ANCOVA With Unequal Slopes: A Simulation Study
(working title). Status: exploratory run read (400 designs x 200 replications, base seed 7301,
201 null designs, 54 of them held out). Criteria below were set from that run with a margin, so a
criterion is met only if the confirmatory run shows a clearly stronger pattern than needed.

## Question, null, analyses
As in PROTOCOL_ancova_v0.md: eight analyses of the effect of a categorical factor A (2 to 4 levels) in a
model with a numeric covariate x and a possible A-by-x interaction. Null: equal adjusted means at the
population grand mean of x. A design is harmful for an analysis if its Type I error exceeds .10 (.075 and
.15 as sensitivity lines). Strata used below: interaction = slope range > 0; covariate imbalance = group
means of x differ; far = grand mean of x at least 2 SD from 0; unequal variances = residual SD ratio > 1.

## Checked before freezing
- validate_ancova.R (300 random designs from the study's generator): the hand-built F tests equal base R
  nested-model comparisons (largest absolute difference in p, 7.2e-14) and car::Anova for Type II, Type III
  (dummy and sum coding, raw and centered covariate) and the HC3 test (white.adjust = "hc3"; 4.8e-13).
- export_for_python.R and validate_python_ancova.py (100 data sets): statsmodels 0.15.0 (Python 3.13.16)
  reproduces all eight analyses, including HC3 (robust="hc3") and the drop-the-interaction rule, with a
  largest absolute difference in p of 6.7e-13.
- Dummy and sum coding give the same Type III test of A when the model has one covariate, so coding is not a
  factor here. These are agreement checks on the implementation, not tests of any hypothesis.

## Hypotheses and numeric criteria (exploratory value in brackets)
H1  Sequential default. anova(lm(y ~ A*x)) (R_default) is harmful in at least 50% of null designs and in
    at least 10 points more designs than any other analysis. [69%; next highest 47%]
H2  Covariate imbalance. Among null designs without an interaction and with covariate imbalance,
    Oneway and R_default are each harmful in at least 50%, whereas TypeII and ANCOVA are harmful in at
    most 20%. [68%, 93%; 12%, 12%]
H3  Raw covariate (hypothesis mismatch). Among null designs with an interaction and a far grand mean,
    TypeIII_raw is harmful in at least 60% and in at least 30 points more designs than TypeIII_ctr.
    [83%; 27%]
H4  Variance mismatch. TypeIII_ctr is harmful in at most 5% of null designs with equal variances and in at
    least 15% of those with unequal variances. [0%; 26%]
H5  HC3 remedy. TypeIII_ctr_HC3 is harmful in at most 10% of all null designs, and in fewer designs than
    each of the other seven analyses. [6%]
H6  Disagreement. The seven analyses 1 to 7 fail to reach a common decision in at least 25% of
    replications. [34%]
H7  Drop-the-interaction practice. Drop_int is harmful in at least 30% of null designs. [47%]
Each hypothesis counts as supported only if every criterion is met in the confirmatory run, with no
allowance for Monte Carlo error, and the same criteria are met in the held-out quarter of designs
(design_id divisible by 4). The Monte Carlo SE of a rate near .05 with 1,000 replications is about .007.

## Mechanism experiment (H8, scaling; criteria set after an exploratory scaling run, seed pool 9100000)
Designs are chosen by rule from a pool of generator draws (20 per set in the confirmatory run, 500
replications each, seed pool 9300000) and every group size is multiplied by 1, 2, 4, 8, 16. Sets: (1) null,
interaction, grand mean of x at least 2 SD from 0, equal variances; (2) null, no interaction, covariate
imbalance, equal variances; (3) null, residual SD ratio at least 2, small groups with large SD (pairing -1),
imbalance at least .4. Criteria on set means [exploratory value at 16x in brackets]:
 a. Set 1: TypeIII_raw error at 16x is at least .60 and rises at each step [.91]; TypeII error at 16x is at
    least .15 and at least .05 above its 1x value [.25 vs .07].
 b. Set 1: TypeIII_ctr error stays within .03 of .05 at every multiplier [.05 to .06].
 c. Set 2: Oneway and R_default error at 16x are each at least .60 [.82, .83]; TypeII, ANCOVA and TypeIII_ctr
    stay within .03 of .05 at every multiplier [.05].
 d. Set 3: TypeIII_ctr error is above .10 at 1x and 16x and differs by less than .05 between them
    [.16, .16]; TypeIII_ctr_HC3 error at 16x is at most .06 and below its 1x value [.053 vs .072].
Script: scaling_ancova.R.

## Runs
Confirmatory: 2,000 designs x 1,000 replications (SIM_MODE=full SIM_SEED=8419).
Extension (addendum written before it runs): 10,000 designs x 1,000 replications, new seed, adds a skewed
covariate and one more factor level layout (up to 5 levels); pooled results weighted toward common layouts.
No outside registration; dates rest on the deposited protocol and code.

## Theory section (planned, not yet built)
Conditional on the covariate values, each F test is a quadratic form in normal errors, so its exact Type I
error can be computed with Imhof's (1961) method for any group sizes, slopes, and variances. Two
propositions parallel the two-way paper: (1) noncentrality grows with N when the tested contrast is nonzero
under the null (hypothesis mismatch; error tends to 1); (2) with equal slopes the pooled-variance test
converges to a fixed error given by a weighted sum of chi-square variates (variance mismatch). To be built
and checked against the simulation before it is claimed.

## Fingerprints at freeze
explore_ancova.R md5 7444efa5d268d653ba5f77e2469a2c5b
scaling_ancova.R md5 7b61374f645a06bb75c7436e1390560a
Seeds: main run base seed 8419 (exploratory 7301); scaling seed pool 9300000 (exploratory 9100000).
Note: explore_ancova.R gained the SIM_SEED environment variable after the exploratory run (default 7301 keeps
the exploratory results reproducible); no test or generator code changed.
Open at freeze: the extension design (one covariate with a skewed covariate and up to 5 levels) is to be
written as a dated addendum before that run; the theory section is not built; the literature gap needs a
fuller search before it is claimed in the manuscript.


---
## Addendum A (14 October 2026, written before the extension run; the confirmatory results above are unchanged)

Purpose: check whether the confirmed pattern holds in a wider design space. Script: explore_ancova_ext.R
(same eight tests as explore_ancova.R, validated against base R and car::Anova on 150 random designs from
the new generator, largest absolute difference in p 7.7e-14). Run: SIM_MODE=ext SIM_SEED=7207, 10,000
designs x 1,000 replications. explore_ancova_ext.R md5 d783bf53055eec8b4c50f56d08b9ca2e.

Changes to the design space: (1) levels of A are 2, 3, 4 and 5, equally likely; (2) the covariate
distribution is normal, skewed or heavy-tailed (equally likely, standardized to mean 0 and SD 1); (3) the
SD of the covariate differs across groups with probability .5 (ratio log-uniform up to 2). Everything else
is as in the main run.

Weighting: pooled results are weighted so the layouts 2, 3 and 4 levels carry the shares of the main run
(.40, .35, .25) and the 5-level layout carries zero weight; the 5-level layout and each covariate shape are
reported separately.

Criteria (no allowance for Monte Carlo error):
E1  H1 to H7 hold with the same numeric criteria in the weighted extension results, and in its held-out quarter
    (design_id divisible by 4).
E2  In each of the four layouts, R_default is the analysis that is harmful in the most null designs, and
    TypeIII_ctr_HC3 is harmful in fewer null designs than TypeIII_ctr.
Descriptive only (no criterion): harmful shares by covariate shape and by smallest group size.
The scaling experiment is not repeated.

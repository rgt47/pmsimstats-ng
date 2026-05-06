## Detailed comments
 * Page #2 (Abstract): "context." -- add future directions of generalization?

 * Page #2 (1. Introduction): "The" -- Is this actually the division we want to highlight front and center?

 * Page #3 (1. Introduction): "not a research-program convention but a paper-specific design choice" -- would skip this part of framing

 * Page #3 (1. Introduction):
   > A comprehensive audit revealed that a parallel tidyverse reimplementation of the same simulation had inadvertently adopted the standard direct mean moderation approach.

   Would seem more straightforward to simply explain that the impact of carryover on power was greater than expected in this implementation vs other similar projects, an exploration of why identified the dynamic explored here

 * Page #3 (1. Introduction):
   > We formalize the distinction, demonstrate its consequences with simulation results, and provide guidance for investigators designing simulation studies for predictive biomarker-moderated treatment effects.

   Question of whether we 1) implement and include in comparisons the hybrid model, vs 2) identify the hybrid as future area of work, vs 3) only note in passing / as something less pragmatically relevant that's just a conceptual reference point...

 * Page #3 (2. Two DGP Architectures): "biomarker" -- baseline biomarker value, I assume?

 * Page #5 (2.2 Architecture B: Differential Correlation (MVN)): "The exponential-decay form is the default" -- Dumb question, but, is this just the default because we wrote it that way? In particular, is this what you added as part of fixing the positive definite problem? Is it actually structurally required by the rest of the setup, or could this decision be dissociated?

 * Page #5 (2.2 Architecture B: Differential Correlation (MVN)):
   > Carryover effects in the DGP inflate the BR means during off-drug periods

   the way this is written I don't understand how its specific to this model - isn't this what carryover of a drug effect always is?

 * Page #7 (3.2 Why the architectures diverge): "the noise structure is independent of drug state" -- this makes the most sense to me as a clear concrete difference - how much is lost if I think of this as the central difference?

 * Page #13 (6.4 Two-stage random slopes): "Two-stage random slopes" -- should we be considering this for our data?


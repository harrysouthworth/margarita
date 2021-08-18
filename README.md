# margarita #

Making texmex more digestible.

When doing extreme value modelling on clinical trial data, we have to account
for baseline effects, which can get messy. margarita provides support functions
to make it less messy and (hopefully) more reliable.

margarita also provides ggplot2 replacement functions for some of the plot
functions in texmex.

The development of margarita is partially funded by AstraZeneca.

**2021-06-15**
Removing lmr and associated functions. Only boxplot.lmr remains. The reason
is that lmr now has its own package.

**2021-08-18**
Fixing a longstanding bug in simulate.margarita.simple that caused those
observations with an actual residual beneath the threshold to have zero
probability of *ever* having a residual benath the threshold. I don't
think this has ever been an issue, but I'm not using the function in
cold blood and this popped up.

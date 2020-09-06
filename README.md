# plant-stan

_I ~~might be have been indoctrinated in~~ __fully embrace__ [Richard Mcelreath](https://xcelab.net/rm/statistical-rethinking/)'s ~~cult~~ __views on statistical inference__: we should not think of it as a set of pre-made tools but instead a set of strategies._

I personally don't like software such as [maxent](https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/maxent) or [HMSC](https://cran.r-project.org/web/packages/Hmsc/index.html). They could be (along with p-values deception) Richard's worst nightmare---black boxes that help you overlook basic aspects of statistical inference. What are my priors? Why did I choose a likelihood or prior distribution over another? Why is this model better than any other? These, among others, are basic questions that any scientist should be able to confidently answer about their models. Therefore, and for the sake of my own learning, I am trying here to compare the performance of the aforementioned software and similar models written in stan (sampled with the help of [rstan](https://cran.r-project.org/web/packages/rstan/index.html) and the R package [rethinking](https://github.com/rmcelreath/rethinking)).



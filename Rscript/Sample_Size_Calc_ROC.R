###  Using sensitivity of biomarkers to calculate
library(pwr)
#sensitivity of biomarker1: 0.52
#sensitivity of biomarker2: 0.56 
#sensitivity of biomarker3: 0.84

h1 <- ES.h(0.63, 0.52) #biomarker 2 vs biomarker 1
pwr.p.test(h = h1, n = NULL, sig.level = 0.05, power = 0.9, alt = "greater")

h1 <- ES.h(0.84, 0.52) #biomarker 3 vs biomarker 1
pwr.p.test(h = h1, n = NULL, sig.level = 0.05, power = 0.9, alt = "greater")


library(pROC)

###  Using full ROC of biomarkers to calculate
roc1 <- roc(<outcome>, <biomarker1>)#fill in as appropriate
roc2 <- roc(<outcome>, <biomarker2>)#fill in as appropriate
roc3 <- roc(<outcome>, <biomarker3>)#fill in as appropriate

## Sample size
# biomarker 2 vs biomarker 1
power.roc.test(roc1, roc2, power=0.9)
# biomarker 3 vs biomarker 1
power.roc.test(roc1, roc3, power=0.9)


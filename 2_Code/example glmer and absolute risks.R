library(lme4)

set.seed(123)
x <- rnorm(1000)
z <- rbinom(1000, 1, .5)
cluster <- rep(1:10, 100)
y <- metamisc:::inv.logit(-4 + x * 2 + z + cluster/2)
y <- rbinom(1000, 1, y)

glmer(y ~ x + z + (1 | cluster), family = binomial(link = "logit")) # random intercepts
glmer(y ~ x + z + (0+x | cluster), family = binomial(link = "logit")) # Random effects
glmer(y ~ x + z + (0+x | cluster), family = binomial(link = "logit")) # Random effects


fitlogit <- glm(y ~ x + z, family = binomial(link = "logit"))

newdata = data.frame(x = rnorm(1000),
                     z = rbinom(1000, 1, .5),
                     cluster = rep(1:10, 100))
newdata$y <- metamisc:::inv.logit(-4 + newdata$x * 2 + newdata$z + newdata$cluster/2)
newdata$y <- rbinom(1000, 1, newdata$y)


predict(fitlogit, newdata = newdata, type = "response") # absolute risks for each person

# Absolute risks for person with x = 0, and mean value for z:
p0 <- predict(fitlogit, newdata = data.frame(x = 0, z = mean(newdata$z)), type = "response")
p0

# Absolute risks for person with x = 1, and mean value for z:
p1 <- predict(fitlogit, newdata = data.frame(x = 1, z = mean(newdata$z)), type = "response") 
p1

# Risk difference:
p0 - p1

# standard errors are more difficult...


# how to get a log model to work (note I get 30 errors, one of which is 
# "algorithm did not converge", meaning that the result cannot be trusted at all)
fitlog_logitstart <- glm(y ~ x + z, family = binomial(link = "log"), start = c(coef(fitlogit)[1],0,0)) # 31 warnings


fitlog_null <- glm(y ~ 1, family = binomial(link = "log"))
fitlog_nullstart <- glm(y ~ x + z, family = binomial(link = "log"), start = c(coef(fitlog_null), 0, 0)) # warnings

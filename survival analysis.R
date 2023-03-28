
# import dataset

library(readr)
mortality <- read_csv("simulated HF mort data for GMPH (1K) final.csv")
View(mortality)


install.packages("survival")
library(ggplot2)
library(survival)


# putting columns into variables

gender <- factor(mortality$gender)
fu_time <- mortality$fu_time
death <- mortality$death


# Kaplan-Meier plot

km_fit <- survfit(Surv(fu_time, death) ~ 1)
plot(km_fit)

summary(km_fit, times = c(1:7,30,60,90*(1:10)))


# splitting the curve byb gender

km_gender_fit <- survfit(Surv(fu_time, death) ~ gender)
plot(km_gender_fit)


# logrank test to compare survaival by gender

survdiff(Surv(fu_time, death) ~ gender, rho = 0)


age <- mortality$age

age <- ifelse(age < 65, "Under 65", "over 65")
table(age)

km_age_fit <- survfit(Surv(fu_time, death) ~ age)
plot(km_age_fit)

survdiff(Surv(fu_time, death) ~ age, rho = 0)


# cox regression

install.packages("survminer")
library(survminer)

ethnicgroup <- mortality$ethnicgroup


# using ethnicgroup as continous
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup)
summary(cox)

# using ethnicgroup as categorical
ethnic_var <- factor(ethnicgroup)

cox2 <- coxph(Surv(fu_time, death) ~ ethnic_var, data = mortality)
summary(cox2)


# creating a new category for the missing values in ethnicgroup

levels(ethnic_var) <- c(levels(ethnic_var),"8")
ethnic_var[is.na(ethnic_var)] <- "8" # Change NA to None

cox3 <- coxph(Surv(fu_time, death) ~ ethnic_var, data = mortality)
summary(cox3)


# Descriptive statistics of some variables

summary(age)

t <- table(gender, exclude = NULL)
addmargins(t) # adds the total
round(100*prop.table(t), digits = 1) # gets % rounded to 1dp

copd <- mortality$copd
c <- table(copd, exclude = NULL)
addmargins(c) # adds the total
round(100*prop.table(c), digits = 1) # gets % rounded to 1dp


prior_dnas <- mortality$prior_dnas
p <- table(prior_dnas, exclude = NULL)
addmargins(p) # adds the total
round(100*prop.table(p), digits = 1) # gets % rounded to 1dp



e <- table(ethnic_var, exclude = NULL)
addmargins(e) # adds the total
round(100*prop.table(e), digits = 1) # gets % rounded to 1dp


# multiple cox model
cox3 <- coxph(Surv(fu_time, death) ~ age + gender + copd + prior_dnas + ethnic_var)
summary(cox3)


# non convergence
quintile <- factor(mortality$quintile)
cox4 <- coxph(Surv(fu_time, death) ~ age + gender + copd +  quintile + ethnic_var)
summary(cox4)

table(quintile, exclude = NULL)

q <- table(quintile, death)
q
round(100*prop.table(q,1),digits=1)



# Fixing non-converging model

# method 1 : relevel
quintile <- relevel(quintile, ref = 2) # quintile 1 as the reference  category
cox5 <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile + ethnic_var)
summary(cox5)

# method 2: combining categories
quintile_5groups <- mortality$quintile
quintile_5groups[quintile_5groups==0] <- 5
quintile_5groups <- factor(quintile_5groups)
table(quintile_5groups, exclude = NULL)
cox6 <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile_5groups + ethnic_var)
summary(cox6)


# method 3: Drop the category

drop_quan_cat <- mortality$quintile
drop_quan_cat[drop_quan_cat==0] <- NA
drop_quan_cat <- factor(drop_quan_cat)
table(drop_quan_cat, exclude = NULL)

cox7 <- coxph(Surv(fu_time, death) ~ age + gender + copd + drop_quan_cat + ethnic_var)
summary(cox7)


# method 4: Drop the variable
cox8 <- coxph(Surv(fu_time, death) ~ age + gender + copd + ethnic_var)
summary(cox8)


# checking proportionality assumptions
cox.zph(fit, transform = "km", global = TRUE)

fit <- coxph(Surv(fu_time, death) ~ gender) # fit the desired model
temp <- cox.zph(fit) # apply the cox.zph function to the desired model
fit2 <- coxph(Surv(fu_time, death) ~ copd)
temp2 <- cox.zph(fit2)
print(temp)
print(temp2)
plot(temp)


#km plot

km_fit <- survfit(Surv(fu_time, death) ~ gender)
plot(km_fit, xlab = "time", ylab = "survival probability")


# Deviance residuals

install.packages("survminer")
library(survminer)

res.cox <- coxph(Surv(fu_time, death) ~ age)
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

res.cox2 <- coxph(Surv(fu_time, death) ~ age)
ggcoxdiagnostics(res.cox2, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxfunctional(Surv(fu_time, death) ~ age + log(age) + sqrt(age))


# time transformation (tt) function
fit3 <- coxph(Surv(fu_time, death) ~ gender + tt(gender))
summary(fit3)

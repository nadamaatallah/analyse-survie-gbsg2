###############################################################
#  ANALYSE DE SURVIE – DONNÉES CANCER DU SEIN (GBSG2)
###############################################################

###############################
# 1. Packages
###############################
install.packages(c(
  "survival",
  "survminer",
  "TH.data",
  "ggplot2",
  "naniar",
  "dplyr",
  "tibble",
  "gridExtra",
  "broom"
))

library(survival)
library(survminer)
library(TH.data)
library(ggplot2)
library(naniar)
library(dplyr)
library(tibble)
library(gridExtra)
library(broom)

###############################
# 2. Chargement des données
###############################
data("GBSG2", package="TH.data")
df <- GBSG2



head_table <- tableGrob(head(df, 10))

str(df)
summary(df)

###############################
# 3. Préparation & Exploration
###############################
p_missing <- gg_miss_var(df) + ggtitle("Valeurs manquantes par variable")


df <- df %>%
  mutate(
    horTh = factor(horTh, levels=c("no","yes")),
    menostat = factor(menostat, levels=c("Pre","Post")),
    tgrade = factor(tgrade, ordered = TRUE)
  )

###############################
# 4. Création de l’objet Surv
###############################
surv_object <- Surv(time = df$time, event = df$cens)

###############################
# 5. Courbe KM Globale
###############################


km_global <- survfit(surv_object ~ 1, data=df)


p_km_global <- ggsurvplot(
  km_global,
  conf.int = TRUE,
  ggtheme = theme_minimal(),
  risk.table = TRUE,
  surv.median.line = "hv",
  palette = "#2E9FDF",
  title = "Courbe de survie Kaplan–Meier – GBSG2"
)


###############################
# 6. Courbes KM par traitement hormonal
###############################
km_group <- survfit(Surv(time, cens) ~ horTh, data=df)

p_km_group <- ggsurvplot(
  km_group,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  linetype = "strata",
  ggtheme = theme_minimal(),
  palette = c("#E64B35", "#4DBBD5"),
  title = "Survie selon traitement hormonal"
)


###############################
# 7. Test Log-Rank
###############################

logrank <- survdiff(Surv(time, cens) ~ horTh, data=df)
logrank
cat("\nP-value = ", 1 - pchisq(logrank$chisq, df=1))


###############################
# 8. Modèles de Cox
###############################
cox_uni <- coxph(Surv(time, cens) ~ horTh, data=df)
cox_multi <- coxph(Surv(time, cens) ~ age + horTh + tgrade + menostat, data=df)

summary(cox_uni)
summary(cox_multi)


###############################
# 9. Forest plot
###############################
p_forest <- ggforest(cox_multi, data=df)

###############################
# 10. Test des Hazards Proportionnels
###############################
cox_zph_test <- cox.zph(cox_multi)

cox_zph_test
print(cox_zph_test)



###############################
# 11. Courbes ajustées
###############################
p_adj <- ggadjustedcurves(
  cox_multi,
  data=df,
  variable="horTh",
  legend.title = "Traitement hormonal",
  palette=c("#E64B35","#4DBBD5")
) + ggtitle("Courbes ajustées – Modèle de Cox multivarié")




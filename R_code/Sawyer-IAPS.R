# R code for "Alcoholism Gender Differences in Brain Responsivity to Emotional Stimuli", 
# by Kayle S. Sawyer, Ph.D., Nasim Maleki, Ph.D., Trinity Urban, M.A., Ksenija Marinkovic, Ph.D., 
# Steven A. Karson, B.A., Susan M. Ruiz, Ph.D., 
# Gordon J. Harris, Ph.D., and Marlene Oscar-Berman, Ph.D.
# https://doi.org/10.1101/428565
# 2018-09-26
# This code is licened Creative Commons CC0: https://creativecommons.org/publicdomain/zero/1.0/

# Working directory should be where sawyer-iaps.Rproj is, one level above R_code, R_data, R_tables, and figures.

# Create R_tables and figures directories
if(!dir.exists(file.path("R_tables"))) { dir.create(file.path("R_tables")) }
if(!dir.exists(file.path("figures"))) { dir.create(file.path("figures")) }

# Load libraries (Pacman will install them automatically if necessary)
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(data.table, tidyverse, skimr, lme4, emmeans, lmerTest)

# Split violin function
# https://stackoverflow.com/a/45614547 (by https://stackexchange.com/users/1645427/yak )
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                               draw_quantiles = NULL, trim = TRUE, scale = "area", 
                               na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#
# Load participant characteristics data
#

taboo.dt <- fread("R_data/Sawyer-IAPS.csv", key = "ID")

# 
# Figure 2
# 

# Melt to long format for plotting
scores <- c("Age", "Education", "VIQ", "PIQ", "WMS_DMI", "HRSD", "DHD", "DD", "LOS")
taboo.dt.long <- melt(taboo.dt, id.vars = c("ID", "Group", "Gender"), measure.vars = scores, 
                      variable.name = "Measure", value.name = "Value")
# Split violin version
# ggplot(data = taboo.dt.long, aes(x = Gender, y = Value, fill = Group)) +
#   geom_split_violin() +
#   geom_point(position = position_jitterdodge(), size = 0.5) +
#   facet_wrap("Measure", scales = "free") + 
#   theme_minimal(base_size = 16, base_family = "Helvetica")

# Boxplot version
fig2 <- ggplot(data = taboo.dt.long, aes(x = Gender, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, colour="light blue", geom="point", 
               shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_wrap("Measure", scales = "free") +
  theme_minimal(base_size = 16, base_family = "Helvetica")
ggsave(filename = "figures/Sawyer-Figure-2.pdf", plot = fig2)


#
# Figure 2-S1 (Participant characteristics table)
#

# Select columns for Figure 2-S1 and make numeric
demog.taboo.dt <- taboo.dt[, .(Group, Gender, Age, Education, VIQ, PIQ, WMS_DMI, HRSD, DHD, DD, LOS)]
demog.taboo.dt[, (scores) := lapply(.SD, as.numeric), .SDcols = scores]

# Inspect data for all participants, then split by Group and Gender
skim(demog.taboo.dt)
demog.taboo.dt %>%
  group_by(Group, Gender) %>%
  skim()

# Make Figure 2-S1
mean.na.rm.t <- function (x) { mean(x, na.rm = T)}
sd.na.rm.t <- function (x) { sd(x, na.rm = T)}
skim_with(numeric = list(mean = mean.na.rm.t, sd = sd.na.rm.t), append = F)
skim_format(numeric = list(digits = 1)) # Doesn't work so used format(round()) as below.
demog.taboo.fig2s1 <- demog.taboo.dt %>%
  group_by(Group, Gender) %>%
  skim()
demog.taboo.fig2s1$value <- format(round(demog.taboo.fig2s1$value, 1), nsmall = 1)
skim_with_defaults()
skim_format_defaults()

# Cast wide and re-sort row order
demog.taboo.fig2s1.wide <- dcast(demog.taboo.fig2s1, variable ~ Group + Gender + stat, value.var="value")
setDT(demog.taboo.fig2s1.wide)[, variable := factor(variable, levels = scores)]
setorder(demog.taboo.fig2s1.wide, variable)

fwrite(demog.taboo.fig2s1.wide, "R_tables/Sawyer-Figure-2-S1.csv")


## T-tests Figure 2-S1 for group differences

#demog.taboo.dt[, lapply(.SD, function(y) t.test(data=.SD, y ~ Group)), .SDcols = scores]
group.ttests.men <- lapply(demog.taboo.dt[Gender=="Men", scores, with=F], 
                           function(y) t.test(y ~ demog.taboo.dt[Gender=="Men", Group]))
pander(group.ttests.men)
group.ttests.women <- lapply(demog.taboo.dt[Gender=="Women", scores, with=F], 
                             function(y) t.test(y ~ demog.taboo.dt[Gender=="Women", Group]))
pander(group.ttests.women)

gender.ttests.alc <- lapply(demog.taboo.dt[Group=="Alcoholic", scores, with=F], 
                            function(y) t.test(y ~ demog.taboo.dt[Group == "Alcoholic", Gender]))
pander(gender.ttests.alc)
gender.ttests.nc <- lapply(demog.taboo.dt[Group=="Nonalcoholic Control", scores, with=F], 
                           function(y) t.test(y ~ demog.taboo.dt[Group == "Nonalcoholic Control", Gender]))
pander(gender.ttests.nc)

interaction.lm <- lapply(demog.taboo.dt[, scores, with = FALSE], function(y) lm(y ~ Group*Gender, data = demog.taboo.dt))
pander(interaction.lm)

#
# Figure 2-S2 Participant characteristics for FSIQ, PIQ, IMI, DMI, WMI, POMS, MAACL
#

# Select columns for Figure 2-S2 and make numeric
demog.taboo.dt <- taboo.dt[, .(Group, Gender, FSIQ, PIQ, WMS_IMI, WMS_DMI, WMS_WMI, 
                               Ten, Dep, Ang, Vig, Fat, Con, 
                               A, D, H, PA, SS, Dys, PASS)]
scores.sup <- c("FSIQ", "PIQ", "WMS_IMI", "WMS_DMI", "WMS_WMI", 
                "Ten", "Dep", "Ang", "Vig", "Fat", "Con", 
                "A", "D", "H", "PA", "SS", "Dys", "PASS")
demog.taboo.dt[, (scores.sup) := lapply(.SD, as.numeric), .SDcols = scores.sup]

skim_with(numeric = list(mean = mean, sd = sd), append = F)
skim_format(numeric = list(digits = 1)) # Doesn't work so used format(round()) as below.
demog.taboo.fig2s2 <- demog.taboo.dt %>%
  group_by(Group, Gender) %>%
  skim()
demog.taboo.fig2s2$value <- format(round(demog.taboo.fig2s2$value, 1), nsmall = 1)
skim_with_defaults()
skim_format_defaults()

# Cast wide and re-sort row order
demog.taboo.fig2s2.wide <- dcast(demog.taboo.fig2s2, variable ~ Group + Gender + stat, value.var="value")
setDT(demog.taboo.fig2s2.wide)[, variable := factor(variable, levels = scores.sup)]
setorder(demog.taboo.fig2s2.wide, variable)
demog.taboo.fig2s2.wide$variable <- c("FSIQ", "PIQ", "WMS IMI", "WMS DMI", "WMS WMI",
                                       "POMS Tension", "POMS Depression", "POMS Anger", "POMS Vigor", "POMS Fatigue", "POMS Confusion",
                                       "MAACL Anxiety", "MAACL Depression", "MAACL Hostility", "MAACL Positive Affect", 
                                       "MAACL Sensation Seeking", "MAACL Dysphoria", "MAACL Positive Affect Sensation Seeking")

fwrite(demog.taboo.fig2s2.wide, "R_tables/Sawyer-Figure-2-S2.csv")


## T-tests Figure 2-S2 for group differences
group.ttests.men <- lapply(demog.taboo.dt[Gender=="Men", scores.sup, with=F], 
                           function(y) t.test(y ~ demog.taboo.dt[Gender=="Men", Group]))
pander(group.ttests.men)

group.ttests.women <- lapply(demog.taboo.dt[Gender=="Women", scores.sup, with=F], 
                             function(y) t.test(y ~ demog.taboo.dt[Gender=="Women", Group]))
pander(group.ttests.women)

gender.ttests.alc <- lapply(demog.taboo.dt[Group=="Alcoholic", scores.sup, with=F], 
                            function(y) t.test(y ~ demog.taboo.dt[Group == "Alcoholic", Gender]))
pander(gender.ttests.alc)

gender.ttests.nc <- lapply(demog.taboo.dt[Group=="Nonalcoholic Control", scores.sup, with=F], 
                           function(y) t.test(y ~ demog.taboo.dt[Group == "Nonalcoholic Control", Gender]))
pander(gender.ttests.nc)

interaction.lm <- lapply(demog.taboo.dt[, scores.sup, with=F], function(y) lm(y ~ Group*Gender, data = demog.taboo.dt))
pander(interaction.lm)


#
# Behavioral ratings (buttonbox) analyses
#

taboo.dt <- fread("R_data/Sawyer-IAPS.csv", key = "ID")
buttonbox.wide <- taboo.dt[behavioral_responses_included==1, grep("ID|Group|Gender|pics", names(taboo.dt)), with = FALSE]


# Separate and melt RT and percent to long format

# Select RT
buttonbox.RT <- buttonbox.wide[, grep("ID|Group|Gender|_RT_", names(buttonbox.wide)), with = FALSE] 

# Melt to long format for ratings
buttonbox.RT.rating <-
  melt.data.table(
    buttonbox.RT,
    id = c("ID", "Group", "Gender"),
    measure.vars = patterns(c("bad$", "good$", "neut$")),
    value.name = c("Bad", "Good", "Neutral"),
    variable.name = "Condition", 
    variable.factor = "true"
  )
buttonbox.RT.rating[, Condition := dplyr::recode_factor(Condition, "1" = "Happy", "2" = "Aversive", 
                                                        "3" = "Neutral", "4" = "Erotic", "5" = "Gruesome")]
# Melt to long format for condition
buttonbox.RT.long <-
  melt.data.table(
    buttonbox.RT.rating,
    id = c("ID", "Group", "Gender", "Condition"),
    measure.vars = c("Bad", "Good", "Neutral"),
    value.name = "RT",
    variable.name = "Rating", 
    variable.factor = "true"
  )


# Select percent
buttonbox.percent <- buttonbox.wide[,grep("ID|Group|Gender|_percent_", names(buttonbox.wide)), with = FALSE] 

# Melt to long format for ratings
buttonbox.percent.rating <-
  melt.data.table(
    buttonbox.percent,
    id = c("ID", "Group", "Gender"),
    measure.vars = patterns(c("bad$", "good$", "neut$")),
    value.name = c("Bad", "Good", "Neutral"),
    variable.name = "Condition", 
    variable.factor = "true"
  )
buttonbox.percent.rating[, Condition := dplyr::recode_factor(Condition, "1" = "Happy", "2" = "Aversive", 
                                                             "3" = "Neutral", "4" = "Erotic", "5" = "Gruesome")]
# Melt to long format for condition
buttonbox.percent.long <-
  melt.data.table(
    buttonbox.percent.rating,
    id = c("ID", "Group", "Gender", "Condition"),
    measure.vars = c("Bad", "Good", "Neutral"),
    value.name = "Percent",
    variable.name = "Rating", 
    variable.factor = "true"
  )

# Merge RT and percent
tb <- merge(buttonbox.percent.long, buttonbox.RT.long, by = c("ID", "Group", "Gender", "Condition", "Rating"))

# Set factors
factors <- c("ID", "Condition", "Rating", "Group", "Gender")
tb[, (factors) := lapply(.SD, as.factor), .SDcols = factors]

# Check
length(tb[,unique(ID)])
skim(tb)
tb[, .N, by = .(ID, Group, Gender)][, .N, by = .(Group, Gender)]
tb[, mean(Percent, na.rm=T), by = .(ID, Group, Gender)]


# Mixed model for percentages
tb.lmer.percent <- lmer(Percent ~ Condition * Rating * Group * Gender + (1|ID), data = tb)
tb.lmer.percent.anova <- anova(tb.lmer.percent)
tb.lmer.percent.emmeans <- emmeans(tb.lmer.percent, c("Rating", "Condition", "Group", "Gender"))

# Save ANOVA and emmeans as tables
sink("R_tables/Sawyer-Figure-3-S1.txt")
print(tb.lmer.percent.anova)
sink()
write.csv(tb.lmer.percent.emmeans, "R_tables/Sawyer-Figure-3-S1-emmeans.csv")


# Mixed model for reaction time
tb.lmer.RT <- lmer(RT ~ Condition * Rating * Group * Gender + (1|ID), data = tb)
tb.lmer.RT.anova <- anova(tb.lmer.RT)
tb.lmer.RT.emmeans <- emmeans(tb.lmer.RT, c("Rating", "Condition", "Group", "Gender"))

# Save ANOVA and emmeans as tables
sink("R_tables/Sawyer-Figure-3-S3.txt")
print(tb.lmer.RT.anova)
sink()
write.csv(tb.lmer.RT.emmeans, "R_tables/Sawyer-Figure-3-S3-emmeans.csv")


## Plot of Buttonbox Ratings (Figure 3)

# Make split violin for percent (Figure 3)
# fig3 <- ggplot(tb, aes(Gender, Percent, fill = Group)) +
#   geom_split_violin() +
#   stat_summary(fun.y=mean, colour="light blue", geom="point", 
#                shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
#   geom_point(position = position_jitterdodge(), size = 0.5) +
#   facet_grid(Condition ~ Rating) +
#   theme_minimal(base_size = 16, base_family = "Helvetica")
# ggsave("figures/Sawyer-Figure-3", fig3)

# Make split violin for reaction time (Figure 3-S2)
# fig3s2 <- ggplot(tb, aes(Gender, RT, fill = Group)) + 
#   geom_split_violin() +
#   stat_summary(fun.y=mean, colour="blue", geom="point", 
#                shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
#   geom_point(position = position_jitterdodge(), size = 0.5) + 
#   facet_grid(Condition ~ Rating) +
#   ylab("Reaction Time (ms)") +
#   theme_minimal(base_size = 16, base_family = "Helvetica")
# ggsave("figures/Sawyer-Figure-3-S2.pdf", fig3s2)

# Make boxplot for percent (Figure 3)
fig3 <- ggplot(tb, aes(Gender, Percent, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun.y=mean, colour="light blue", geom="point", 
               shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid(Condition ~ Rating) +
  theme_minimal(base_size = 16, base_family = "Helvetica")
ggsave(filename = "figures/Sawyer-Figure-3.pdf", plot = fig3)

# Make boxplot for reaction time (Figure 3-S2)
fig3s2 <- ggplot(tb, aes(Gender, RT, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, colour="light blue", geom="point", 
               shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid(Condition ~ Rating) +
  ylab("Reaction Time") +
  # coord_cartesian(ylim = c(500, 3500)) + #hides outliers!
  theme_minimal(base_size = 16, base_family = "Helvetica")
ggsave(filename = "figures/Sawyer-Figure-3-S2.pdf", plot = fig3s2)


# Split violin boxplots
# ggplot(tb, aes(Gender, Percent, fill = Group)) +
#   geom_split_violin(alpha = 0.5, color = "white") +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
#   stat_summary(fun.y=mean, colour="light blue", geom="point", 
#                shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
#   geom_point(position = position_jitterdodge(), size = 0.5) +
#   facet_grid(Condition ~ Rating) +
#   theme_minimal(base_size = 16, base_family = "Helvetica")
# 
# ggplot(tb, aes(Gender, RT, fill = Group)) +
#   geom_split_violin(alpha = 0.5, color = "white") +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#   stat_summary(fun.y=mean, colour="light blue", geom="point", 
#                shape="diamond", size=6, position = position_dodge(width = 0.75), show.legend = FALSE) + 
#   geom_point(position = position_jitterdodge(), size = 0.5) +
#   facet_grid(Condition ~ Rating) +
#   ylab("Reaction Time") +
#   # coord_cartesian(ylim = c(500, 3500)) + #hides outliers!
#   theme_minimal(base_size = 16, base_family = "Helvetica")


#
# Split violin plots of fMRI cluster activation
#

# Erotic vs Neutral, mni305 volume (Figure 4-S3)
# Select clusters and melt
erovneu <- taboo.dt[, .(ID, Group, Gender, `Limbic Structures`)]
erovneu[, GroupGender := paste(Group, Gender)]

erovneu.long <- melt.data.table(erovneu, variable.name = "Region", value.name = "CES", 
                                id.vars = c("ID", "Group", "Gender", "GroupGender"))
# Recode and re-order levels
erovneu.long[, Group := dplyr::recode_factor(Group, Alcoholic = "ALC", `Nonalcoholic Control` = "NC")]
erovneu.long[, Gender := dplyr::recode_factor(Gender, Men = "Men", Women = "Women")]

# Make plot
fig4s3 <- ggplot(erovneu.long, aes(Gender, `CES`, fill = Group)) +
  geom_split_violin() +
  geom_point(position = position_jitterdodge()) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Limbic Structures")

ggsave("figures/Sawyer-Figure-4-S3.pdf", fig4s3)

# Boxplot version
ggplot(erovneu.long, aes(Gender, `CES`, fill = Group)) +
  # geom_split_violin() +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, colour="yellow", geom="point", 
               shape=18, size=4, position = position_dodge(width = 0.75), show.legend = FALSE) + 
  ggtitle("Limbic Structures")


# Negative vs Neutral, left hemisphere surface (Figure 5-S4)
negvneu <- taboo.dt[, .(ID, Group, Gender, caudalmiddlefrontal1, caudalmiddlefrontal2, superiorfrontal, inferiorparietal)]
negvneu[, GroupGender := paste(Group, Gender)]

# Check N sample
negvneu[, .N, by = .(Group, Gender)]

negvneu.long <- melt.data.table(negvneu, variable.name = "Region", value.name = "CES", 
                                id.vars = c("ID", "Group", "Gender", "GroupGender"))

# Recode and re-order levels
negvneu.long[, Group := dplyr::recode_factor(Group, Alcoholic = "ALC", `Nonalcoholic Control` = "NC")]
negvneu.long[, Gender := dplyr::recode_factor(Gender, Men = "Men", Women = "Women")]


cluster_names <- c(
  'caudalmiddlefrontal1'="Caudal Middle\nFrontal 1",
  'caudalmiddlefrontal2'="Caudal Middle\nFrontal 2",
  'superiorfrontal'="Superior\nFrontal",
  'inferiorparietal'="Inferior\nParietal"
)
cluster_labeller <- function(variable,value){
  return(cluster_names[value])
}

# Make plot
fig5s4 <- ggplot(negvneu.long, aes(Gender, `CES`, fill = Group)) +
  geom_split_violin() +
  geom_point(position = position_jitterdodge()) +
  facet_grid( ~ Region, labeller = as_labeller(cluster_names)) +
  theme_minimal(base_size = 16, base_family = "Helvetica")

ggsave("figures/Sawyer-Figure-5-S4.pdf", fig5s4)

# Boxplot version
ggplot(negvneu.long, aes(Gender, `CES`, fill = Group)) +
  # geom_split_violin() +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid( ~ Region, labeller = as_labeller(cluster_names)) +
  stat_summary(fun.y=mean, colour="yellow", geom="point", 
               shape=18, size=4, position = position_dodge(width = 0.75), show.legend = FALSE) + 
  theme_minimal(base_size = 16, base_family = "Helvetica")

#
# Brain behavior analyses (relationships of fMRI with buttonbox data and participant characteristics)
#

# Merge datasets and remove missing data listwise (N=76)
tb[, ID:= as.factor(ID)]
setkeyv(tb, "ID")
taboo.dt[, ID := as.factor(ID)]
setkeyv(taboo.dt, "ID")

tb.brain <- tb[taboo.dt[, .(ID, `Limbic Structures`, 
                            caudalmiddlefrontal1, caudalmiddlefrontal2, 
                            superiorfrontal, inferiorparietal,
                            PIQ, VIQ, WMS_DMI, Dep, A, SS)]]
tb.brain <- tb.brain[!is.na(Group)]


# Mixed models

# Percent ratings for limbic regions cluster (erotic vs. neutral)
tb.lmer.percent.limbic <- lmer(Percent ~ `Limbic Structures` * Group * Gender * Rating + (1|ID), 
                               data = tb.brain[Condition=="Erotic"])
anova(tb.lmer.percent.limbic) # No significant interactions between fMRI activity and Group or Gender for percent
# emtrends(tb.lmer.percent.limbic, pairwise ~ Group*Gender|Rating, var = "Limbic Structures")
# 
# plot.percent.limbic <- ggplot(data = tb.brain[Condition=="Erotic"], 
#                       aes(x = `Limbic Structures`, y = Percent, color = Rating)) +
#   geom_point() + 
#   geom_smooth(method = 'lm') +
#   ggtitle("Erotic vs. Neutral") + 
#   facet_grid(Group~Gender)
# plot.percent.limbic

# Reaction times (RT) for limbic regions cluster (erotic vs. neutral)

tb.lmer.RT.limbic <- lmer(RT ~ `Limbic Structures` * Group * Gender * Rating + (1|ID), 
                          data = tb.brain[Condition=="Erotic"])
anova(tb.lmer.RT.limbic) 
# No significant interactions between fMRI activity and Group or Gender for RT
emtrends(tb.lmer.RT.limbic, pairwise ~ Group*Gender|Rating, var = "Limbic Structures")

plot.RT.limbic <- ggplot(data = tb.brain[Condition=="Erotic"],
                         aes(x = `Limbic Structures`, y = RT, color = Rating)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggtitle("Erotic vs. Neutral") +
  facet_grid(Group~Gender)
plot.RT.limbic

# Percent ratings for caudalmiddlefrontal1 cluster (aversive vs. neutral)
tb.lmer.percent.caudalmiddlefrontal1 <- lmer(Percent ~ caudalmiddlefrontal1 * Group * Gender * Rating + (1|ID), 
                                             data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.percent.caudalmiddlefrontal1) 
# Significant interaction: fMRI*Rating*Group*Gender for percent
emtrends(tb.lmer.percent.caudalmiddlefrontal1, pairwise ~ Group*Gender|Rating, var = "caudalmiddlefrontal1")
# However, slopes not significantly different between groups

plot.percent.caudalmiddlefrontal1 <- ggplot(data = tb.brain[Condition=="Aversive"],
                                            aes(x = caudalmiddlefrontal1, y = Percent, color = Rating)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggtitle("Aversive vs. Neutral") +
  facet_grid(Group~Gender)
plot.percent.caudalmiddlefrontal1

plot.percent.caudalmiddlefrontal1 <- ggplot(data = tb.brain[Condition=="Aversive", GroupGender := paste(Group, Gender)][Condition=="Aversive"],
                                            aes(x = caudalmiddlefrontal1, y = Percent, color = GroupGender)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggtitle("Aversive vs. Neutral") +
  facet_grid(~Rating)
plot.percent.caudalmiddlefrontal1
# Perhaps driven by two NCM outliers?

# RT ratings for caudalmiddlefrontal1 cluster (aversive vs. neutral)
tb.lmer.RT.caudalmiddlefrontal1 <- lmer(RT ~ caudalmiddlefrontal1 * Group * Gender * Rating + (1|ID), 
                                        data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.RT.caudalmiddlefrontal1) 
# No significant interactions between fMRI activity and Group or Gender for RT
# emtrends(tb.lmer.RT.caudalmiddlefrontal1, pairwise ~ Group*Gender|Rating, var = "caudalmiddlefrontal1")

# plot.RT.caudalmiddlefrontal1 <- ggplot(data = tb.brain[Condition=="Aversive"],
#                                        aes(x = caudalmiddlefrontal1, y = RT, color = Rating)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   ggtitle("Aversive vs. Neutral") +
#   facet_grid(Group~Gender)
# plot.RT.caudalmiddlefrontal1
# 
# plot.RT.caudalmiddlefrontal1 <- ggplot(data = tb.brain[Condition=="Aversive", GroupGender := paste(Group, Gender)][Condition=="Aversive"],
#                                        aes(x = caudalmiddlefrontal1, y = RT, color = GroupGender)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   ggtitle("Aversive vs. Neutral") +
#   facet_grid(~Rating)
# plot.RT.caudalmiddlefrontal1


# Percent ratings for caudalmiddlefrontal2 cluster (aversive vs. neutral)
tb.lmer.percent.caudalmiddlefrontal2 <- lmer(Percent ~ caudalmiddlefrontal2 * Group * Gender * Rating + (1|ID), 
                                             data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.percent.caudalmiddlefrontal2) 
# No significant interactions between fMRI activity and Group or Gender for percent
# emtrends(tb.lmer.percent.caudalmiddlefrontal2, pairwise ~ Group*Gender|Rating, var = "caudalmiddlefrontal2")

# RT ratings for caudalmiddlefrontal2 cluster (aversive vs. neutral)
tb.lmer.RT.caudalmiddlefrontal2 <- lmer(RT ~ caudalmiddlefrontal2 * Group * Gender * Rating + (1|ID), 
                                        data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.RT.caudalmiddlefrontal2) 
# No significant interactions between fMRI activity and Group or Gender for RT
# emtrends(tb.lmer.RT.caudalmiddlefrontal2, pairwise ~ Group*Gender|Rating, var = "caudalmiddlefrontal2")

# Percent ratings for superiorfrontal cluster (aversive vs. neutral)
tb.lmer.percent.superiorfrontal <- lmer(Percent ~ superiorfrontal * Group * Gender * Rating + (1|ID), 
                                        data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.percent.superiorfrontal) 
# Significant interaction: fMRI*Rating*Group*Gender for percent
emtrends(tb.lmer.percent.superiorfrontal, pairwise ~ Group*Gender|Rating, var = "superiorfrontal")
# However, slopes not significantly different between groups
plot.percent.superiorfrontal <- ggplot(data = tb.brain[Condition=="Aversive", GroupGender := paste(Group, Gender)][Condition=="Aversive"],
                                       aes(x = superiorfrontal, y = Percent, color = GroupGender)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggtitle("Aversive vs. Neutral") +
  facet_grid(~Rating)
plot.percent.superiorfrontal

# RT ratings for superiorfrontal cluster (aversive vs. neutral)
tb.lmer.RT.superiorfrontal <- lmer(RT ~ superiorfrontal * Group * Gender * Rating + (1|ID), 
                                   data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.RT.superiorfrontal) 
# No significant interactions between fMRI activity and Group or Gender for RT
# emtrends(tb.lmer.RT.superiorfrontal, pairwise ~ Group*Gender|Rating, var = "superiorfrontal")
# plot.RT.superiorfrontal <- ggplot(data = tb.brain[Condition=="Aversive", GroupGender := paste(Group, Gender)][Condition=="Aversive"],
#                                   aes(x = superiorfrontal, y = RT, color = GroupGender)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   ggtitle("Aversive vs. Neutral") +
#   facet_grid(~Rating)
# plot.RT.superiorfrontal

# Percent ratings for inferiorparietal cluster (aversive vs. neutral)
tb.lmer.percent.inferiorparietal <- lmer(Percent ~ inferiorparietal * Group * Gender * Rating + (1|ID), 
                                         data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.percent.inferiorparietal) 
# No significant interactions between fMRI activity and Group or Gender for percent
# emtrends(tb.lmer.percent.inferiorparietal, pairwise ~ Group*Gender|Rating, var = "inferiorparietal")

# RT ratings for inferiorparietal cluster (aversive vs. neutral)
tb.lmer.RT.inferiorparietal <- lmer(RT ~ inferiorparietal * Group * Gender * Rating + (1|ID), 
                                    data = tb.brain[Condition=="Aversive"])
anova(tb.lmer.RT.inferiorparietal) 
# No significant interactions between fMRI activity and Group or Gender for percent
# emtrends(tb.lmer.RT.inferiorparietal, pairwise ~ Group*Gender|Rating, var = "inferiorparietal")



#
# Relationships of affect and ability measures with buttonbox data
# 

# Create Z-score of each measure and check outliers
measures <- c("PIQ", "VIQ", "WMS_DMI", "Dep", "A", "SS")

taboo.dt[ID %in% tb[,ID]]
taboo.dt.long <- melt(data = taboo.dt[ID %in% tb[,ID]], id.vars = "ID", measure.vars = measures, 
                      variable.name = "measure", value.name = "value")
taboo.dt.long[, measureZ := scale(value), by=measure]
taboo.dt.outliers <- taboo.dt.long[(measureZ > 3) | (measureZ < -3)]
for (measure.i in unique(taboo.dt.outliers[,measure])) {
  assign(paste0(measure.i, ".out"), taboo.dt.outliers[measure==measure.i, ID])
}
taboo.dt[ID %in% A.out, .(Group, Gender)]
taboo.dt[ID %in% Dep.out, .(Group, Gender)]

# Mixed models for Percent

tb.lmer.percent.VIQ <- lmer(Percent ~ VIQ * Group * Gender * Condition * Rating + (1|ID), 
                            data = tb.brain)
anova(tb.lmer.percent.VIQ) 
# emtrends(tb.lmer.percent.VIQ, pairwise ~ Group*Gender|Rating, var = "VIQ")

tb.lmer.percent.PIQ <- lmer(Percent ~ PIQ * Group * Gender * Condition * Rating + (1|ID), 
                            data = tb.brain)
anova(tb.lmer.percent.PIQ) 
# emtrends(tb.lmer.percent.PIQ, pairwise ~ Group*Gender|Rating, var = "PIQ")

tb.lmer.percent.WMS_DMI <- lmer(Percent ~ WMS_DMI * Group * Gender * Condition * Rating + (1|ID), 
                                data = tb.brain)
anova(tb.lmer.percent.WMS_DMI) 
# emtrends(tb.lmer.percent.WMS_DMI, pairwise ~ Group*Gender|Rating, var = "WMS_DMI")

tb.lmer.percent.Dep <- lmer(Percent ~ Dep * Group * Gender * Condition * Rating + (1|ID), 
                            data = tb.brain[!(ID %in% Dep.out)])
anova(tb.lmer.percent.Dep) 
emtrends(tb.lmer.percent.Dep, pairwise ~ Group|Rating, var = "Dep")
plot.percent.Dep <- ggplot(data = tb.brain[!(ID %in% Dep.out)], 
                           aes(x = VIQ, y = RT, color = Rating)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(Group~Gender)
plot.percent.Dep

tb.lmer.percent.A <- lmer(Percent ~ A * Group * Gender * Condition * Rating + (1|ID), 
                          data = tb.brain[!(ID %in% A.out)])
anova(tb.lmer.percent.A) 
# emtrends(tb.lmer.percent.A, pairwise ~ Group*Gender|Rating, var = "A")

tb.lmer.percent.SS <- lmer(Percent ~ SS * Group * Gender * Condition * Rating + (1|ID), 
                           data = tb.brain)
anova(tb.lmer.percent.SS) 
# emtrends(tb.lmer.percent.SS, pairwise ~ Group*Gender|Rating, var = "SS")


# Mixed models for RT

tb.lmer.RT.VIQ <- lmer(RT ~ VIQ * Group * Gender * Condition * Rating + (1|ID), 
                       data = tb.brain)
anova(tb.lmer.RT.VIQ) # Significant interaction: VIQ*Group*Gender
emtrends(tb.lmer.RT.VIQ, pairwise ~ Group*Gender, var = "VIQ") # No comparisons significant
plot.RT.VIQ <- ggplot(data = tb.brain, 
                      aes(x = VIQ, y = RT)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(Group~Gender)
plot.RT.VIQ

tb.lmer.RT.PIQ <- lmer(RT ~ PIQ * Group * Gender * Condition * Rating + (1|ID), 
                       data = tb.brain)
anova(tb.lmer.RT.PIQ) 
# emtrends(tb.lmer.RT.PIQ, pairwise ~ Group*Gender|Rating, var = "PIQ")

tb.lmer.RT.WMS_DMI <- lmer(RT ~ WMS_DMI * Group * Gender * Condition * Rating + (1|ID), 
                           data = tb.brain)
anova(tb.lmer.RT.WMS_DMI) # Significant interaction: WMS_DMI*Group*Gender*Rating
emtrends(tb.lmer.RT.WMS_DMI, pairwise ~ Group*Gender|Rating, var = "WMS_DMI")
plot.RT.WMS_DMI <- ggplot(data = tb.brain, 
                          aes(x = WMS_DMI, y = RT, color = Rating)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(Group~Gender)
plot.RT.WMS_DMI

tb.lmer.RT.Dep <- lmer(RT ~ Dep * Group * Gender * Condition * Rating + (1|ID), 
                       data = tb.brain[!(ID %in% Dep.out)])
anova(tb.lmer.RT.Dep) # Significant interaction: Dep*Group*Gender*Rating
emtrends(tb.lmer.RT.Dep, pairwise ~ Group|Condition*Rating, var = "Dep")
plot.RT.Dep <- ggplot(data = tb.brain, 
                      aes(x = Dep, y = RT, color = Rating)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(Group~Condition)
plot.RT.Dep
# Looks heavily influenced by outliers.

tb.lmer.RT.A <- lmer(RT ~ A * Group * Gender * Condition * Rating + (1|ID), 
                     data = tb.brain)
anova(tb.lmer.RT.A) # Significant interaction: A*Group*Gender*Rating
emtrends(tb.lmer.RT.A, pairwise ~ Group*Gender|Rating, var = "A")
plot.RT.A <- ggplot(data = tb.brain, 
                    aes(x = A, y = RT, color = Rating)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(Group~Gender)
plot.RT.A

tb.lmer.RT.SS <- lmer(RT ~ SS * Group * Gender * Condition * Rating + (1|ID), 
                      data = tb.brain)
anova(tb.lmer.RT.SS) 
emtrends(tb.lmer.RT.SS, pairwise ~ Group*Gender|Rating, var = "SS")

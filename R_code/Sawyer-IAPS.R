# R code for "Alcoholism Gender Differences in Brain Responsivity to Emotional Stimuli", 
# by Kayle S. Sawyer, Ph.D., Nasim Maleki, Ph.D., Trinity Urban, M.A., Ksenija Marinkovic, Ph.D., 
# Steven A. Karson, B.A., Susan M. Ruiz, Ph.D., 
# Gordon J. Harris, Ph.D., and Marlene Oscar-Berman, Ph.D.
# https://doi.org/10.1101/428565
# 2018-09-26
# This code is licened Creative Commons CC0: https://creativecommons.org/publicdomain/zero/1.0/

# Working directory should be where sawyer-iaps.Rproj is, one level above R_code, R_data, R_tables, and figures.

# Load libraries (Pacman will install them automatically if necessary)
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(data.table, tidyverse, skimr, lme4, emmeans, lmerTest)

#
# Load participant characteristics data
#

taboo.dt <- fread("R_data/Sawyer-IAPS.csv", key = "ID")

#
# Table 1
#

# Create R_tables directory
if(!dir.exists(file.path("R_tables"))) { dir.create(file.path("R_tables")) }

# Select columns for Table 1 and make numeric
demog.taboo.dt <- taboo.dt[, .(Group, Gender, Age, Education, VIQ, PIQ, WMS_DMI, HRSD, DHD, DD, LOS)]
scores <- c("Age", "Education", "VIQ", "PIQ", "WMS_DMI", "HRSD", "DHD", "DD", "LOS")
demog.taboo.dt[, (scores) := lapply(.SD, as.numeric), .SDcols = scores]

# Inspect data for all participants, then split by Group and Gender
skim(demog.taboo.dt)
demog.taboo.dt %>%
  group_by(Group, Gender) %>%
  skim()

# Make Table 1
skim_with(numeric = list(mean = mean, sd = sd), append = F)
skim_format(numeric = list(digits = 1)) # Doesn't work so used format(round()) as below.
demog.taboo.table1 <- demog.taboo.dt %>%
  group_by(Group, Gender) %>%
  skim()
demog.taboo.table1$value <- format(round(demog.taboo.table1$value, 1), nsmall = 1)
skim_with_defaults()
skim_format_defaults()

# Cast wide and re-sort row order
demog.taboo.table1.wide <- dcast(demog.taboo.table1, variable ~ Group + Gender + stat, value.var="value")
setDT(demog.taboo.table1.wide)[, variable := factor(variable, levels = scores)]
setorder(demog.taboo.table1.wide, variable)

fwrite(demog.taboo.table1.wide, "R_tables/table1.csv")


## T-tests Table 1 for group differences

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

## MAKE Table 1-S1 FSIQ, PIQ, IMI, DMI, WMI, POMS, MAACL

# Select columns for Table 1-S1 and make numeric
demog.taboo.dt <- taboo.dt[, .(Group, Gender, FSIQ, PIQ, WMS_IMI, WMS_DMI, WMS_WMI, 
                               Ten, Dep, Ang, Vig, Fat, Con, 
                               A, D, H, PA, SS, Dys, PASS)]
scores.sup <- c("FSIQ", "PIQ", "WMS_IMI", "WMS_DMI", "WMS_WMI", 
                "Ten", "Dep", "Ang", "Vig", "Fat", "Con", 
                "A", "D", "H", "PA", "SS", "Dys", "PASS")
demog.taboo.dt[, (scores.sup) := lapply(.SD, as.numeric), .SDcols = scores.sup]

skim_with(numeric = list(mean = mean, sd = sd), append = F)
skim_format(numeric = list(digits = 1)) # Doesn't work so used format(round()) as below.
demog.taboo.tables1 <- demog.taboo.dt %>%
  group_by(Group, Gender) %>%
  skim()
demog.taboo.tables1$value <- format(round(demog.taboo.tables1$value, 1), nsmall = 1)
skim_with_defaults()
skim_format_defaults()

# Cast wide and re-sort row order
demog.taboo.tables1.wide <- dcast(demog.taboo.tables1, variable ~ Group + Gender + stat, value.var="value")
setDT(demog.taboo.tables1.wide)[, variable := factor(variable, levels = scores.sup)]
setorder(demog.taboo.tables1.wide, variable)
demog.taboo.tables1.wide$variable <- c("FSIQ", "PIQ", "WMS IMI", "WMS DMI", "WMS WMI",
                                       "POMS Tension", "POMS Depression", "POMS Anger", "POMS Vigor", "POMS Fatigue", "POMS Confusion",
                                       "MAACL Anxiety", "MAACL Depression", "MAACL Hostility", "MAACL Positive Affect", 
                                       "MAACL Sensation Seeking", "MAACL Dysphoria", "MAACL Positive Affect Sensation Seeking")
  
fwrite(demog.taboo.tables1.wide, "R_tables/tables1-s1.csv")


## T-tests Table 1-S1 for group differences
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
anova(tb.lmer.percent)
emmeans(tb.lmer.percent, c("Rating", "Condition", "Group", "Gender"))

# Mixed model for reaction time
tb.lmer.RT <- lmer(RT ~ Condition * Rating * Group * Gender + (1|ID), data = tb)
anova(tb.lmer.RT)
emmeans(tb.lmer.RT, c("Rating", "Condition", "Group", "Gender"))

# # Five mixed models, one for each condition
# conditions <- c("Aversive", "Erotic", "Gruesome", "Happy", "Neutral")
# for(c in conditions){
#   print(paste0("Percent analysis. Condition = ", c))
#   tb.lmer.percent <- lmer(Percent ~ Rating * Group * Gender + (1|ID), data = tb[Condition==c])
#   print(anova(tb.lmer.percent))
#   print(emmeans(tb.lmer.percent, c("Group", "Gender")))
#   
#   print(paste0("RT analysis. Condition = ", c))
#   tb.lmer.RT <- lmer(RT ~ Rating * Group * Gender + (1|ID), data = tb[Condition==c])
#   print(anova(tb.lmer.RT))
#   print(emmeans(tb.lmer.RT, c("Group", "Gender")))
# }


## Plot of Buttonbox Ratings (Figure 2)

# https://stackoverflow.com/a/45614547 (by https://stackoverflow.com/users/1870254/jan-glx )
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
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
                               draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# Make Split Violin for percent (Figure 2)
fig2 <- ggplot(tb, aes(Gender, Percent, fill = Group)) +
  geom_split_violin() +
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid(Condition ~ Rating) +
  theme_minimal(base_size = 16, base_family = "Helvetica")

ggsave("figures/Sawyer-Figure-2.pdf", fig2)

# Make Split Violin for reaction time (Figure 2-S2)
fig2s2 <- ggplot(tb, aes(Gender, RT, fill = Group)) + 
  geom_split_violin() +
  geom_point(position = position_jitterdodge(), size = 0.5) + 
  facet_grid(Condition ~ Rating) +
  ylab("Reaction Time (ms)") +
  theme_minimal(base_size = 16, base_family = "Helvetica")

ggsave("figures/Sawyer-Figure-2-S2.pdf", fig2s2)

# Boxplot versions
ggplot(tb, aes(Gender, Percent, fill = Group)) +
  # geom_split_violin() +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid(Condition ~ Rating) +
  theme_minimal(base_size = 16, base_family = "Helvetica")

ggplot(tb, aes(Gender, RT, fill = Group)) +
  # geom_split_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid(Condition ~ Rating) +
  ylab("Reaction Time") +
  # coord_cartesian(ylim = c(500, 3500)) + #hides outliers!
  theme_minimal(base_size = 16, base_family = "Helvetica")

#
# Split violin plots of cluster activation
#

# Erotic vs Neutral, mni305 volume (Figure 3-S3)
# Select clusters and melt
erovneu <- taboo.dt[, .(ID, Group, Gender, `Limbic Structures`)]
erovneu[, GroupGender := paste(Group, Gender)]

erovneu.long <- melt.data.table(erovneu, variable.name = "Region", value.name = "CES", 
                                id.vars = c("ID", "Group", "Gender", "GroupGender"))
# Recode and re-order levels
erovneu.long[, Group := dplyr::recode_factor(Group, Alcoholic = "ALC", `Nonalcoholic Control` = "NC")]
erovneu.long[, Gender := dplyr::recode_factor(Gender, Men = "Men", Women = "Women")]

# Make plot
fig3s3 <- ggplot(erovneu.long, aes(Gender, `CES`, fill = Group)) +
  geom_split_violin() +
  geom_point(position = position_jitterdodge()) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Limbic Structures")

ggsave("figures/Sawyer-Figure-3-S3.pdf", fig3s3)

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


# Negative vs Neutral, left hemisphere surface (Figure 4-S4)
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
fig4s4 <- ggplot(negvneu.long, aes(Gender, `CES`, fill = Group)) +
  geom_split_violin() +
  geom_point(position = position_jitterdodge()) +
  facet_grid( ~ Region, labeller = as_labeller(cluster_names)) +
  theme_minimal(base_size = 16, base_family = "Helvetica")

ggsave("figures/Sawyer-Figure-4-S4.pdf", fig4s4)

# Boxplot version
ggplot(negvneu.long, aes(Gender, `CES`, fill = Group)) +
  # geom_split_violin() +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  facet_grid( ~ Region, labeller = as_labeller(cluster_names)) +
  stat_summary(fun.y=mean, colour="yellow", geom="point", 
               shape=18, size=4, position = position_dodge(width = 0.75), show.legend = FALSE) + 
  theme_minimal(base_size = 16, base_family = "Helvetica")


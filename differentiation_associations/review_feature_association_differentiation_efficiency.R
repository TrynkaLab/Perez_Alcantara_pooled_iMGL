#Association of cell line features with neuronal differentiation efficiency
library(lme4)
library(lmerTest)
library(tidyverse)

dir = "review_feature_association_differentiation_efficiency/"
dir.create(paste0("../../data/results/",dir))

donor_metadata = read_csv("../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv")
metadata_info = read_csv("../../../OTAR2065_differentiation_efficiency/data/metadata_info_hipsci_IPMAR.csv") %>%
  mutate(Disease.Status = case_when(Disease.Status == "Control" ~ "Normal",
                                    Disease.Status == "EOAD" ~ "AD",
                                    Disease.Status == "LOAD" ~ "AD",
                                    .default = "Normal"))
merged_metadata = donor_metadata %>%
  left_join(metadata_info) %>% 
  distinct() %>%
  mutate(line = case_when(donor == "Arlene" ~ "Arlene-003",
                          donor == "Bertha" ~ "Bertha-004",
                          donor == "Cindy" ~ "Cindy-005",
                          donor == "Dexter" ~ "Dexter-006",
                          donor == "Fiona" ~ "Fiona-010",
                          donor == "Gilma" ~ "Gilma-009",
                          donor == "Hector" ~ "Hector-011",
                          donor == "Imani" ~ "Imani-012",
                          donor == "Javier" ~ "Javier-013",
                          donor == "Keoni" ~ "Keoni-014",
                          donor == "Mindy" ~ "Mindy-016",
                          donor == "Nestor" ~ "Nestor-017",
                          donor == "Olaf" ~ "Olaf-018",
                          donor == "Qiana" ~ "Qiana-022",
                          .default = line),
         Disease.Status = case_when(Disease.Status == "Control" ~ "Normal",
                                    Disease.Status == "EOAD" ~ "AD",
                                    Disease.Status == "LOAD" ~ "AD",
                                    .default = "Normal"))

rm(metadata_info)
rm(donor_metadata)
## differentiation efficiency
ipsc_premac =  read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_premac_iPSC.csv")
premac_microglia = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv")
old_young_premac = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_old_vs_young_premac.csv")
## aggregate metadata
ipsc_premac =  ipsc_premac %>%
  left_join(merged_metadata) %>%
  group_by(pool) %>%
  mutate(n_lines_per_pool = n_distinct(line)) %>%
  ungroup()

premac_microglia =  premac_microglia %>%
  left_join(merged_metadata) %>%
  group_by(pool, treatment) %>%
  mutate(n_lines_per_pool = n_distinct(line)) %>%
  ungroup()

old_young_premac = old_young_premac %>%
  left_join(merged_metadata) %>%
  group_by(pool) %>%
  mutate(n_lines_per_pool = n_distinct(line)) %>%
  ungroup()

## sanity checks

length(unique(premac_microglia$donor)) # 221
length(unique(premac_microglia$line)) # 241
length(unique(old_young_premac$line)) # 220
length(unique(premac_microglia$line)) # 241

table(premac_microglia$data_source)

# HipSci  IPMAR 
# 2094     130


### Measuring drivers of differentiation efficiency

ipsc_premac$line<- as.factor(ipsc_premac$line)
ipsc_premac$sex<- as.factor(ipsc_premac$sex)
ipsc_premac$pool<- as.factor(ipsc_premac$pool)

old_young_premac$line = as.factor(old_young_premac$line)
old_young_premac$sex = as.factor(old_young_premac$sex)
old_young_premac$pool = as.factor(old_young_premac$pool)

premac_microglia$line = as.factor(premac_microglia$line)
premac_microglia$sex = as.factor(premac_microglia$sex)
premac_microglia$pool = as.factor(premac_microglia$pool)

model <- lmer(
  scaled_log_fraction ~ preMAC_age +  
    #Disease.Status + 
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    n_lines_per_pool +
    (1 | pool) + (1 | line),
  data = ipsc_premac
)

plot(model)

model_no_pool <- lmer(
  scaled_log_fraction ~ preMAC_age + Disease.Status + sex   + n_lines_per_pool + (1 | line),
  data = ipsc_premac,
  REML = FALSE
)

model_with_pool <- lmer(
  scaled_log_fraction ~ preMAC_age + Disease.Status + sex  + n_lines_per_pool+ (1 | pool) + (1 | line),
  data = ipsc_premac,
  REML = FALSE
)

anova(model_no_pool, model_with_pool)
summary(model_with_pool)
# all the pool variation is captured in the lines, because there are many that are unique per pool
# subsetting to shared lines only

shared_ipsc_premac = ipsc_premac %>%
  group_by(line) %>%
  filter(n_distinct(pool) >= 2) %>%
  ungroup()

length(unique(shared_ipsc_premac$line)) # 89

model <- lmer(
  scaled_log_fraction ~ preMAC_age +  
    #Disease.Status + # can't fit this (all IPMAR that remain are normal)
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + (1 | pool) + (1 | line),
  data = shared_ipsc_premac
)

plot(model) # much better than with non-shared lines due to reduced collinearity between lines and pools

model_no_pool <- lmer(
  scaled_log_fraction ~ preMAC_age + 
    #Disease.Status + 
    sex   + n_lines_per_pool +
    (1 | line),
  data = shared_ipsc_premac,
  REML = FALSE
)

model_with_pool <- lmer(
  scaled_log_fraction ~ preMAC_age + 
    n_lines_per_pool +
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    (1 | line) + 
    pool,
  data = shared_ipsc_premac,
  REML = FALSE
)

drop1(model_with_pool, test = "Chisq")



# scaled_log_fraction ~ preMAC_age + n_lines_per_pool + Pluritest.pluripotency.score + Pluritest.novelty.score + sex + (1 | line) + pool
# Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# preMAC_age                    0.0369 0.03692     1 149.889  0.1112 0.739293   
# n_lines_per_pool                                                              
# Pluritest.pluripotency.score  1.8198 1.81983     1  78.895  5.4787 0.021777 * 
#   Pluritest.novelty.score       0.5465 0.54646     1  79.411  1.6451 0.203353   
# sex                           0.1849 0.18487     1  78.239  0.5566 0.457886   
# pool                         12.1740 0.93647    13 175.444  2.8193 0.001082 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#removes each term one at a time
#performs likelihood ratio tests
#gives a table of p-values


model_with_line <- lmer(
  scaled_log_fraction ~ preMAC_age + 
    n_lines_per_pool +
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    (1 | pool) + 
    line,
  data = shared_ipsc_premac,
  REML = FALSE
)

drop1(model_with_line, test = "Chisq")
# Model:
# scaled_log_fraction ~ preMAC_age + Pluritest.pluripotency.score + Pluritest.novelty.score + sex + n_lines_per_pool + (1 | pool) + line
# Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)    
# preMAC_age                     0.047 0.04663     1  12.528  0.2042 0.6591    
# Pluritest.pluripotency.score                                                 
# Pluritest.novelty.score                                                      
# sex                                                                          
# n_lines_per_pool               0.188 0.18848     1  12.205  0.8253 0.3812    
# line                         124.974 1.60223    78 218.383  7.0157 <2e-16 ***

summary(model_with_pool)
# changing direction of estimates for pools means that the effect on the proliferation is not constant
# other significant is Pluritest.pluripotency.score  with Estimate =  0.016239, so it increases proliferation
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                   -1.693234   6.843766 154.121737  -0.247 0.804918    
# preMAC_age                    -0.070322   0.210919 149.889151  -0.333 0.739293    
# n_lines_per_pool               0.050704   0.020349 160.280855   2.492 0.013731 *  
#   Pluritest.pluripotency.score   0.016239   0.006938  78.895187   2.341 0.021777 *  
#   Pluritest.novelty.score        0.612396   0.477454  79.410562   1.283 0.203353    
# sexMale                        0.131051   0.175665  78.239237   0.746 0.457886    
# poolpool11                     1.713487   0.503868 174.774258   3.401 0.000833 ***
#   poolpool13                     2.124678   0.711051 160.774556   2.988 0.003248 ** 
#   poolpool14                     1.553628   0.666303 164.721107   2.332 0.020926 *  
#   poolpool15                     1.273726   0.862289 156.973175   1.477 0.141640    
# poolpool16                     1.928331   0.620149 164.174933   3.109 0.002210 ** 
#   poolpool17                     1.391342   0.456586 164.678030   3.047 0.002691 ** 
#   poolpool2                      1.143231   1.150836 155.850529   0.993 0.322059    
# poolpool3                      1.929750   0.567809 170.133840   3.399 0.000844 ***
#   poolpool4                      1.799931   0.470621 165.432443   3.825 0.000185 ***
#   poolpool5                      1.555273   0.593223 167.929676   2.622 0.009552 ** 
#   poolpool6                      1.289847   0.510244 171.536688   2.528 0.012378 *  
#   poolpool7                      1.156441   0.578837 159.894958   1.998 0.047426 *  
#   poolpool8                      1.534171   0.660300 161.899504   2.323 0.021399 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## old vs young premacs #######


shared_old_young_premac = old_young_premac %>%
  group_by(line) %>%
  filter(n_distinct(pool) >= 2) %>%
  ungroup()

length(unique(shared_old_young_premac$line)) # 73


model <- lmer(
  scaled_log_fraction ~ n_lines_per_pool +
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + (1 | pool) + (1 | line),
  data = shared_old_young_premac
)

plot(model) # much better than with non-shared lines due to reduced collinearity between lines and pools



model_with_pool <- lmer(
  scaled_log_fraction ~ 
    n_lines_per_pool +
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    pool + 
    (1 | line),
  data = shared_old_young_premac,
  REML = FALSE
)

drop1(model_with_pool, test = "Chisq")
#Single term deletions using Satterthwaite's method:
# 
# Model:
# scaled_log_fraction ~ n_lines_per_pool + Pluritest.pluripotency.score + Pluritest.novelty.score + sex + pool + (1 | line)
# Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# n_lines_per_pool                                                           
# Pluritest.pluripotency.score 0.7367 0.73669     1  62.675  1.6379 0.20532  
# Pluritest.novelty.score      0.0014 0.00139     1  64.572  0.0031 0.95579  
# sex                          0.2895 0.28947     1  61.989  0.6436 0.42547  
# pool                         9.7765 0.75204    13 146.693  1.6721 0.07264 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(model_with_pool)
# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: scaled_log_fraction ~ Pluritest.pluripotency.score + Pluritest.novelty.score +      sex + pool + (1 | line)
#    Data: shared_old_young_premac
# 
#       AIC       BIC    logLik -2*log(L)  df.resid 
#     528.6     593.8    -244.3     488.6       173 
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -3.03087 -0.52562  0.03308  0.54091  2.09926 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  line     (Intercept) 0.5469   0.7395  
#  Residual             0.4498   0.6706  
# Number of obs: 193, groups:  line, 67
# 
# Estimate Std. Error         df t value Pr(>|t|)  
# (Intercept)                   -0.631973   2.313921 177.616109  -0.273   0.7851  
# n_lines_per_pool              -0.009140   0.074189 141.772224  -0.123   0.9021  
# Pluritest.pluripotency.score   0.011293   0.008824  62.674799   1.280   0.2053  
# Pluritest.novelty.score        0.032853   0.590279  64.572127   0.056   0.9558  
# sexMale                       -0.178391   0.222362  61.988639  -0.802   0.4255  
# poolpool11                     0.308746   0.584178 148.880301   0.529   0.5979  
# poolpool13                     0.756547   0.791353 146.584669   0.956   0.3406  
# poolpool14                     0.185761   0.564036 153.829391   0.329   0.7423  
# poolpool15                    -0.185647   1.090125 138.495847  -0.170   0.8650  
# poolpool16                     0.513441   0.448885 147.009342   1.144   0.2546  
# poolpool17                     0.659392   0.286016 166.451978   2.305   0.0224 *
#   poolpool2                      0.208680   0.765158 143.059718   0.273   0.7855  
# poolpool3                      0.290457   0.790324 142.717581   0.368   0.7138  
# poolpool4                      0.250670   0.522201 152.290096   0.480   0.6319  
# poolpool5                      0.353663   0.804428 145.502984   0.440   0.6608  
# poolpool6                      0.193523   0.528541 149.613444   0.366   0.7148  
# poolpool7                      0.856692   0.524954 145.701743   1.632   0.1049  
# poolpool8                      0.357009   0.576107 143.326278   0.620   0.5364  
# # ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

model_with_line <- lmer(
  scaled_log_fraction ~ n_lines_per_pool + 
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    (1 | pool) + 
    line,
  data = shared_old_young_premac,
  REML = FALSE
)

drop1(model_with_line, test = "Chisq")
# Model:
# scaled_log_fraction ~ n_lines_per_pool + Pluritest.pluripotency.score + Pluritest.novelty.score + sex + (1 | pool) + line
# Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)    
# n_lines_per_pool               0.205 0.20455     1  14.461  0.6556 0.4312    
# Pluritest.pluripotency.score                                                 
# Pluritest.novelty.score                                                      
# sex                                                                          
# line                         117.755 1.86913    63 177.165  5.9905 <2e-16 ***

### microglia vs premac


shared_premac_microglia = premac_microglia %>%
  group_by(line) %>%
  filter(n_distinct(pool) >= 2) %>%
  ungroup()

table(shared_premac_microglia$pool,shared_premac_microglia$treatment)
length(unique(shared_premac_microglia$line)) # 89

model <- lmer(
  scaled_log_fraction ~ preMAC_age + 
    n_lines_per_pool+
    treatment + 
    #Disease.Status + # can't fit this (all IPMAR that remain are normal)
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + (1 | pool) + (1 | line),
  data = shared_premac_microglia
)

plot(model) # much better than with non-shared lines due to reduced collinearity between lines and pools
length(unique(shared_premac_microglia$line)) # 89

model_with_pool <- lmer(
  scaled_log_fraction ~ preMAC_age + 
    treatment + 
    n_lines_per_pool+
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    pool + 
    (1 | line),
  data = shared_premac_microglia,
  REML = FALSE
)

drop1(model_with_pool, test = "Chisq")

# scaled_log_fraction ~ preMAC_age + treatment + Pluritest.pluripotency.score + Pluritest.novelty.score + sex + pool + (1 | line)
# Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# preMAC_age                    0.0172  0.0172     1 1234.50  0.0285 0.8659595    
# treatment                     1.0543  0.5272     2 1219.53  0.8748 0.4172210    
# n_lines_per_pool              4.5282  4.5282     1 1233.23  7.5140 0.0062107 ** 
#   Pluritest.pluripotency.score  0.0656  0.0656     1   76.91  0.1089 0.7423315    
# Pluritest.novelty.score       0.1561  0.1561     1   79.06  0.2591 0.6121556    
# sex                           1.5750  1.5750     1   77.81  2.6135 0.1100061    
# pool                         25.6331  1.8309    14 1290.90  3.0382 0.0001217 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(model_with_pool)
# Fixed effects:
# Estimate Std. Error         df t value Pr(>|t|)   
# (Intercept)                  -1.197e+00  6.949e-01  1.789e+02  -1.722  0.08676 . 
# preMAC_age                   -5.478e-04  3.244e-03  1.235e+03  -0.169  0.86596   
# treatmentLPS                 -6.586e-02  4.983e-02  1.219e+03  -1.322  0.18646   
# treatmentuntreated           -3.044e-02  6.205e-02  1.220e+03  -0.491  0.62375   
# n_lines_per_pool              1.926e-02  7.026e-03  1.233e+03   2.741  0.00621 **
#   Pluritest.pluripotency.score  1.777e-03  5.386e-03  7.691e+01   0.330  0.74233   
# Pluritest.novelty.score       1.898e-01  3.730e-01  7.906e+01   0.509  0.61216   
# sexMale                      -2.219e-01  1.373e-01  7.781e+01  -1.617  0.11001   
# poolpool11                    3.720e-01  2.529e-01  1.295e+03   1.471  0.14158   
# poolpool13                    8.318e-01  2.698e-01  1.282e+03   3.083  0.00210 **
#   poolpool14                    6.598e-01  2.694e-01  1.286e+03   2.449  0.01445 * 
#   poolpool15                    2.104e-01  3.409e-01  1.248e+03   0.617  0.53720   
# poolpool16                    5.728e-01  2.392e-01  1.292e+03   2.395  0.01676 * 
#   poolpool17                    3.433e-01  1.889e-01  1.302e+03   1.817  0.06945 . 
# poolpool2                     5.468e-01  3.017e-01  1.260e+03   1.812  0.07019 . 
# poolpool3                     5.779e-01  2.887e-01  1.269e+03   2.002  0.04553 * 
#   poolpool4                     3.580e-01  2.500e-01  1.279e+03   1.432  0.15239   
# poolpool5                     2.638e-01  2.766e-01  1.281e+03   0.954  0.34047   
# poolpool6                     3.489e-01  2.534e-01  1.289e+03   1.377  0.16881   
# poolpool7                     5.889e-01  2.365e-01  1.283e+03   2.490  0.01292 * 
#   poolpool8                     7.212e-01  2.516e-01  1.271e+03   2.867  0.00421 **
#   poolpool9                     5.083e-01  1.591e-01  1.270e+03   3.195  0.00143 **
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

model_with_line <- lmer(
  scaled_log_fraction ~ preMAC_age + 
    n_lines_per_pool + 
    treatment + 
    Pluritest.pluripotency.score + Pluritest.novelty.score + 
    sex  + 
    (1 | pool) + 
    line,
  data = shared_premac_microglia,
  REML = FALSE
)

drop1(model_with_line, test = "Chisq")

# Model:
#   scaled_log_fraction ~ preMAC_age + n_lines_per_pool + treatment + Pluritest.pluripotency.score + Pluritest.novelty.score + sex + (1 | pool) + line
# Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)    
# preMAC_age                     0.00  0.0001     1 1301.21  0.0002 0.98983    
# n_lines_per_pool               3.50  3.4958     1   29.45  6.1189 0.01937 *  
#   treatment                      1.05  0.5225     2 1287.28  0.9146 0.40093    
# Pluritest.pluripotency.score                                                 
# Pluritest.novelty.score                                                      
# sex                                                                          
# line                         386.26  4.9520    78 1203.12  8.6677 < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
### fitting a t-test for the AD lines vs the rest, because they can't be included in the joint model
## same for data_source (HipSci or IPMAR)
t.test(scaled_log_fraction ~ Disease.Status, data = ipsc_premac)
length(unique(ipsc_premac$line)) # 238


# data:  scaled_log_fraction by Disease.Status
# t = -0.15627, df = 10.638, p-value = 0.8787
# alternative hypothesis: true difference in means between group AD and group Normal is not equal to 0
# 95 percent confidence interval:
#   -0.6891501  0.5981315
# sample estimates:
#   mean in group AD mean in group Normal 
# -0.02734717           0.01816210 

t.test(scaled_log_fraction ~ Disease.Status, data = old_young_premac)
length(unique(old_young_premac$line)) # 220

# data:  scaled_log_fraction by Disease.Status
# t = -0.39814, df = 10.511, p-value = 0.6985
# alternative hypothesis: true difference in means between group AD and group Normal is not equal to 0
# 95 percent confidence interval:
#   -0.9127354  0.6344503
# sample estimates:
#   mean in group AD mean in group Normal 
# -0.136061685          0.003080835 

t.test(scaled_log_fraction ~ Disease.Status, data = premac_microglia)
length(unique(premac_microglia$line)) # 241

# 
# data:  scaled_log_fraction by Disease.Status
# t = 0.32553, df = 90.054, p-value = 0.7455
# alternative hypothesis: true difference in means between group AD and group Normal is not equal to 0
# 95 percent confidence interval:
#   -0.1887326  0.2627043
# sample estimates:
#   mean in group AD mean in group Normal 
# 0.032183276         -0.004802581 


## data source (HipSci vs IPMAR)

t.test(scaled_log_fraction ~ data_source, data = ipsc_premac)
# data:  scaled_log_fraction by data_source
# t = -0.5794, df = 17.741, p-value = 0.5696
# alternative hypothesis: true difference in means between group HipSci and group IPMAR is not equal to 0
# 95 percent confidence interval:
#   -0.6132528  0.3483403
# sample estimates:
#   mean in group HipSci  mean in group IPMAR 
# 0.01095698           0.14341322 

t.test(scaled_log_fraction ~ data_source, data = old_young_premac)
# data:  scaled_log_fraction by data_source
# t = -0.20076, df = 17.841, p-value = 0.8432
# alternative hypothesis: true difference in means between group HipSci and group IPMAR is not equal to 0
# 95 percent confidence interval:
#   -0.5567738  0.4597026
# sample estimates:
#   mean in group HipSci  mean in group IPMAR 
# -0.003888611          0.044646989 

t.test(scaled_log_fraction ~ data_source, data = premac_microglia)
# data:  scaled_log_fraction by data_source
# t = -1.566, df = 143.9, p-value = 0.1195
# alternative hypothesis: true difference in means between group HipSci and group IPMAR is not equal to 0
# 95 percent confidence interval:
#   -0.32564374  0.03774303
# sample estimates:
#   mean in group HipSci  mean in group IPMAR 
# -0.01180337           0.13214699 






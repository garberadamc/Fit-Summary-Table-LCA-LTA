
# Creating Summary Fit Tables for LCA and LTA Analyses Using `MplusAutomation`
# Adam Garber

# This syntax automates the process of generating a model fit summary table for comparing a series of LCA models generated during enumeration. The code produces a publication read table that approximately adheres to APA formatting guidelines. Included in the table is the exhaustive set of fit indices recommended as current best practice for making enumeration decisions based on simulation results (Nylund et al., 2007). A separate example is provided below to demonstrate the procedure for an LTA analysis with two time points. 

# If using this tutorial to produce tables for publication it would be greatly appreciated if you cite this resource using the citation provided here:
  
# Garber, A. C. (2021). Creating Summary Fit Tables for LCA and LTA Analyses Using MplusAutomation. Retrieved from psyarxiv.com/uq2fh

# --------------------------------------------------------------------------------------

library(MplusAutomation) # Conduit between R & Mplus
library(tidyverse)       # Data manipulation
library(here)            # Location, location, location
library(gt)              # Tables with seamless knitting

# --------------------------------------------------------------------------------------

## Create Model Fit Summary Table for Latent Class Analysis (LCA)

# --------------------------------------------------------------------------------------

output_enum <- readModels(here("enum_mplus"), quiet = TRUE)

# --------------------------------------------------------------------------------------

enum_extract <- LatexSummaryTable(output_enum,                                 
                keepCols=c("Title", "Parameters", "LL", "BIC", "aBIC",
                           "BLRT_PValue", "T11_VLMR_PValue","Observations"))

# --------------------------------------------------------------------------------------

allFit <- enum_extract %>% 
  mutate(aBIC = -2*LL+Parameters*log((Observations+2)/24)) %>% 
  mutate(CIAC = -2*LL+Parameters*(log(Observations)+1)) %>% 
  mutate(AWE = -2*LL+2*Parameters*(log(Observations)+1.5)) %>%
  mutate(SIC = -.5*BIC) %>% 
  mutate(expSIC = exp(SIC - max(SIC))) %>% 
  mutate(BF = exp(SIC-lead(SIC))) %>% 
  mutate(cmPk = expSIC/sum(expSIC)) %>% 
  select(1:5,9:10,6:7,13,14) %>% 
  arrange(Parameters)

# --------------------------------------------------------------------------------------

allFit %>% 
  mutate(Title = str_remove(Title, " LCA Enumeration - Youth Coping Strategies")) %>% 
  gt() %>%
  tab_header(
    title = md("**Model Fit Summary Table**"), subtitle = md("&nbsp;")) %>% 
  cols_label(
    Title = "Classes",
    Parameters = md("Par"),
    LL = md("*LL*"),
    T11_VLMR_PValue = "VLMR",
    BLRT_PValue = "BLRT",
    BF = md("BF"),
    cmPk = md("*cmP_k*")) %>%
  tab_footnote(
    footnote = md(
    "*Note.* Par = parameters; *LL* = log likelihood;
      BIC = bayesian information criterion;
      aBIC = sample size adjusted BIC;
      CAIC = consistent Akaike information criterion;
      AWE = approximate weight of evidence criterion;
      BLRT = bootstrapped likelihood ratio test p-value;
      VLMR = Vuong-Lo-Mendell-Rubin adjusted likelihood ratio test p-value;
      cmPk = approximate correct model probability."), 
    locations = cells_title()) %>% 
  tab_options(column_labels.font.weight = "bold") %>% 
  fmt_number(10,decimals = 2,
             drop_trailing_zeros=TRUE,
             suffixing = TRUE) %>% 
  fmt_number(c(3:9,11),
             decimals = 1) %>% 
  fmt_missing(1:11,
              missing_text = "--") %>% 
  fmt(c(8:9,11),
    fns = function(x) 
    ifelse(x < 0.001, "<.001",
           scales::number(x, accuracy = 0.01))) %>%
  fmt(10, fns = function(x) 
    ifelse(x>100, ">100",
           scales::number(x, accuracy = .1))) 

# --------------------------------------------------------------------------------------

# Create Model Fit Summary Table for Latent Transition Analysis (LTA)

# --------------------------------------------------------------------------------------

# TIME 1
output_enum_t1 <- readModels(here("enum_LCA_time1"), quiet = TRUE)

# TIME 2
output_enum_t2 <- readModels(here("enum_LCA_time2"), quiet = TRUE)

# --------------------------------------------------------------------------------------

enum_extract1 <- LatexSummaryTable(output_enum_t1,                                 
                keepCols = c("Title", "Parameters", "LL", "BIC", "aBIC",
                           "BLRT_PValue", "T11_VLMR_PValue","Observations"))   

enum_extract2 <- LatexSummaryTable(output_enum_t2,                                 
                keepCols = c("Title", "Parameters", "LL", "BIC", "aBIC",
                           "BLRT_PValue", "T11_VLMR_PValue","Observations")) 

# --------------------------------------------------------------------------------------
                           
allFit1 <- enum_extract1 %>% 
  mutate(aBIC = -2*LL+Parameters*log((Observations+2)/24)) %>% 
  mutate(CIAC = -2*LL+Parameters*(log(Observations)+1)) %>% 
  mutate(AWE = -2*LL+2*Parameters*(log(Observations)+1.5)) %>%
  mutate(SIC = -.5*BIC) %>% 
  mutate(expSIC = exp(SIC - max(SIC))) %>% 
  mutate(BF = exp(SIC-lead(SIC))) %>% 
  mutate(cmPk = expSIC/sum(expSIC)) %>% 
  select(1:5,9:10,6:7,13,14) %>% 
  arrange(Parameters)

allFit2 <- enum_extract2 %>% 
  mutate(aBIC = -2*LL+Parameters*log((Observations+2)/24)) %>% 
  mutate(CIAC = -2*LL+Parameters*(log(Observations)+1)) %>% 
  mutate(AWE = -2*LL+2*Parameters*(log(Observations)+1.5)) %>%
  mutate(SIC = -.5*BIC) %>% 
  mutate(expSIC = exp(SIC - max(SIC))) %>% 
  mutate(BF = exp(SIC-lead(SIC))) %>% 
  mutate(cmPk = expSIC/sum(expSIC)) %>% 
  select(1:5,9:10,6:7,13,14) %>% 
  arrange(Parameters)

allFit <- full_join(allFit1,allFit2)

# --------------------------------------------------------------------------------------

allFit %>% 
  mutate(Title = str_remove(Title, "_Time*")) %>% 
  gt() %>%
  tab_header(
    title = md("**Model Fit Summary Table**"), subtitle = md("&nbsp;")) %>% 
  cols_label(
    Title = "Classes",
    Parameters = md("Par"),
    LL = md("*LL*"),
    T11_VLMR_PValue = "VLMR",
    BLRT_PValue = "BLRT",
    BF = md("BF"),
    cmPk = md("*cmPk*")) %>%
  tab_footnote(
    footnote = md(
    "*Note.* Par = Parameters; *LL* = model log likelihood;
      BIC = Bayesian information criterion;
      aBIC = sample size adjusted BIC; CAIC = consistent Akaike information criterion;
      AWE = approximate weight of evidence criterion;
      BLRT = bootstrapped likelihood ratio test p-value;
      VLMR = Vuong-Lo-Mendell-Rubin adjusted likelihood ratio test p-value;
      *cmPk* = approximate correct model probability."), 
    locations = cells_title()) %>% 
  tab_options(column_labels.font.weight = "bold") %>% 
  fmt_number(10,decimals = 2,
             drop_trailing_zeros=TRUE,
             suffixing = TRUE) %>% 
  fmt_number(c(3:9,11), 
             decimals = 1) %>% 
  fmt_missing(1:11,
              missing_text = "--") %>% 
  fmt(c(8:9,11),
    fns = function(x) 
    ifelse(x<0.001, "<.001",
           scales::number(x, accuracy = 0.01))) %>%
  fmt(10,
      fns = function(x) 
      ifelse(x>100, ">100",
             scales::number(x, accuracy = .1))) %>%
  tab_row_group(
    group = "Time-1",
    rows = 1:6) %>%
  tab_row_group(
    group = "Time-2",
    rows = 7:12) %>% 
  row_group_order(groups = c("Time-1","Time-2"))

# --------------------------------------------------------------------------------------

## References

# Hallquist, Michael N., and Joshua F. Wiley. 2018. "MplusAutomation: An R Package for Facilitating Large-Scale Latent Variable Analyses in Mplus." Structural Equation Modeling, 1--18. <https://doi.org/10.1080/10705511.2017.1402334>.

# Nylund, K. L., Asparouhov, T., & Muthén, B. O. (2007). Deciding on the number of classes in latent class analysis and growth mixture modeling: A Monte Carlo simulation study. Structural equation modeling: A multidisciplinary Journal, 14(4), 535-569.

# R Core Team. 2019.R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. <https://www.R-project.org/>.

# Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). "Welcome to the tidyverse." Journal of Open Source Software, 4(43), 1686. doi: 10.21105/joss.01686.

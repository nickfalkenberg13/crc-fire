
# function for printing model results with robust SEs
model_sum_rb <- function(pm, ci_round = 4){
  # pm is a poisson model fit using glm
  model_tib <- broom::tidy(coeftest(pm, vcov. = vcovHC(pm, type = 'HC3')), 
                           conf.int = T)
  model_tib$estimate_exp <- exp(model_tib$estimate)
  model_tib$lower_exp <- exp(model_tib$conf.low)
  model_tib$upper_exp <- exp(model_tib$conf.high)
  model_tib$ci <- paste("(", round(model_tib$conf.low, ci_round), ", ", round(model_tib$conf.high, ci_round), ")",
                        sep = "")
  model_tib$ci_exp <- paste("(", round(model_tib$lower_exp, 3), ", ", round(model_tib$upper_exp, 3), ")",
                            sep = "")
  
  return(model_tib[, c("term",
                       "estimate",
                       "ci",
                       "estimate_exp",
                       "ci_exp",
                       "p.value")])
}

# function for creating model formulas
mod_formula <- function(trms, resp_var = "count"){
  str_trms <- paste(trms, collapse = " + ")
  str_form <- paste(resp_var, str_trms, sep = " ~ ")
  
  return(formula(str_form))
}

# function to run models
get_effect_estimate <- function(A, covars, d){
  mod_terms <- c(A, covars)
  frmla <- mod_formula(mod_terms)
  
  psn_mod <- glm(frmla, data = d, family = poisson)
  print(paste("Model results for exposure", A))
  effect_est <- model_sum_rb(psn_mod) %>% filter(grepl(A, term))
  print("--------------------------")
  
  return(effect_est)
}
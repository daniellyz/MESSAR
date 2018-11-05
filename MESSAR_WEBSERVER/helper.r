mcc_calculator<-function(rules){
  P_A_B = rules$support/rules$Py
  tpr = rules$support
  fpr = (1 - P_A_B)*rules$Py
  fnr = (1 - rules$confidence)*rules$Px
  tnr = 1 - tpr - fpr - fnr
  MCC = (tpr*tnr - fpr*fnr)/sqrt((tpr+fpr)*(tpr+fnr)*(tnr+fpr)*(tnr+fnr))
  return(MCC)
}

search_rules<-function(ref_feature, ref_feature_type, mass, mass_diff, ppm_search){
  
  # Combine masses:
  
  test_feature = c(mass, mass_diff)
  test_feature_type = c(rep('mass', length(mass)), rep("mass_diff", length(mass_diff)))
  NF = length(ref_feature)
    
  # Search ref features:
  
  matched_feature = c()
  for (i in 1:NF){

    valid = 1
    feature = ref_feature[[i]]
    feature_type = ref_feature_type[[i]]
    
    for (f in 1:length(feature)){
      errors = abs((feature[f]-test_feature)/feature[f])*1000000
      min_error = min(errors)
      index_error = which.min(errors)

      if ((min_error <= ppm_search) & (test_feature_type[index_error] == feature_type[f])){
      } else {
        valid = 0
        break
      }
    }
    if (valid == 1){ # If always valid
        matched_feature = c(matched_feature, i)
    }
  }
  return(matched_feature)}
  
ppm_calc <- function(x, ref){
  errors = abs(x-ref)/ref*1000000
  return(errors)
}

compute_similariy<-function(candidates,true_structure){
  
  # Compare list of candidates (SMILE code) with the true structure
  candidates_finger <- lapply(parse.smiles(candidates), get.fingerprint, type='maccs')
  true_structure_finger <- get.fingerprint(parse.smiles(true_structure)[[1]],type="maccs")
  fp.distance=sapply(candidates_finger,function(x) 1-distance(x,true_structure_finger))
  return(fp.distance)
}

fn <- function(x, y){
  
  # Compare list of candidates 
  x_finger <- try(get.fingerprint(parse.smiles(x)[[1]],type='maccs'),silent=T)
  y_finger <-  try(get.fingerprint(parse.smiles(y)[[1]],type='maccs'),silent=T)
  if ((class(x_finger)!="try-error") && (class(y_finger)!="try-error")){
    fp.distance=1-distance(x_finger,y_finger)} else {
      fp.distance = 1}
  return(fp.distance)
}

numextract <- function(string){as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))} 

aggregate_rules <- function(rules_extracted){

  # Group rules by their substructure
  
  combined_type = aggregate(rules_extracted$SPECTRAL_FEATURE_TYPE, list(rules_extracted$SUBSTRUCTURE), toString)[,2]
  combined_features = aggregate(rules_extracted$SPECTRAL_FEATURE, list(rules_extracted$SUBSTRUCTURE), toString)[,2]
  combined_type = sapply(combined_type, clean_feature_text)
  combined_features = sapply(combined_features, clean_feature_text)

  sum_confidence = aggregate(rules_extracted$CONFIDENCE, list(rules_extracted$SUBSTRUCTURE), sum)
  sum_lift = aggregate(rules_extracted$LIFT, list(rules_extracted$SUBSTRUCTURE), sum)[,2]
  sum_mcc = aggregate(rules_extracted$MCC, list(rules_extracted$SUBSTRUCTURE), sum)[,2]
  sum_aggregated = cbind.data.frame(sum_confidence,sum_lift,sum_mcc)
  colnames(sum_aggregated) = c("SUBSTRUCTURE", "CONFIDENCE", "LIFT","MCC")

  median_confidence = aggregate(rules_extracted$CONFIDENCE, list(rules_extracted$SUBSTRUCTURE), median)
  median_lift = aggregate(rules_extracted$LIFT, list(rules_extracted$SUBSTRUCTURE), median)[,2]
  median_mcc = aggregate(rules_extracted$MCC, list(rules_extracted$SUBSTRUCTURE), median)[,2]
  median_aggregated = cbind.data.frame(median_confidence,median_lift,median_mcc)
  colnames(median_aggregated) = c("SUBSTRUCTURE", "CONFIDENCE", "LIFT","MCC")
  
  max_confidence = aggregate(rules_extracted$CONFIDENCE, list(rules_extracted$SUBSTRUCTURE), max)
  max_lift = aggregate(rules_extracted$LIFT, list(rules_extracted$SUBSTRUCTURE), max)[,2]
  max_mcc = aggregate(rules_extracted$MCC, list(rules_extracted$SUBSTRUCTURE), max)[,2]
  max_aggregated = cbind.data.frame(max_confidence,max_lift,max_mcc)
  colnames(max_aggregated) = c("SUBSTRUCTURE", "CONFIDENCE", "LIFT","MCC")
  
  return(list(type=combined_type, features=combined_features, 
              sum_aggregated=sum_aggregated, median_aggregated=median_aggregated, max_aggregated=max_aggregated))
}

eval_rules_aggregated <- function(rule_aggregated, aggregate_type, score_type, pb_threshold){

  # Pick a type of aggregated rule evaluation
  
  if (aggregate_type=="SUM"){score = rule_aggregated$sum_aggregated[, score_type]}
  if (aggregate_type=="MEDIAN"){score = rule_aggregated$median_aggregated[, score_type]}
  if (aggregate_type=="MAX"){score = rule_aggregated$max_aggregated[, score_type]}
  
  aggregated = cbind.data.frame(rule_aggregated$sum_aggregated$SUBSTRUCTURE, score)
  colnames(aggregated) = c("SUBSTRUCTURE", score_type)

  threshold = quantile(aggregated[,2], probs=pb_threshold/100, na.rm=T)
  aggregated =  aggregated[aggregated[,2]>=threshold, ]
  aggregated =  aggregated[order(aggregated[,2], decreasing = T),]

  return(aggregated)  
}

img_uri <- function(x) {sprintf('<img src="%s"/>', knitr::image_uri(x))}

clean_feature_text <- function(feature){
  feature = str_replace_all(feature, "\\[", "")
  feature = str_replace_all(feature, " ", "")
  feature = str_replace_all(feature, "\\]", "")
  feature = str_replace_all(feature, "\\'","")
  return(feature)
}

fdr_compute <- function(rules, decoy, fdr_thr){
  
  combined_label = c(rep("Target", nrow(rules)), rep("Decoy", nrow(decoy)))
  combined_mcc = c(rules$MCC, decoy$MCC)
  combined_lift = c(rules$LIFT, decoy$LIFT)
  NC = length(combined_label)
  NG = 100 # Number of grids
  grid_mcc = seq(min(combined_mcc), max(combined_mcc), length.out = NG)
  grid_lift = seq(min(combined_lift), max(combined_lift), length.out = NG)
  fdr_matrix = matrix(1, NG, NG)
  valid_matrix  = matrix(0, NG, NG) # Nb of valid (target) rules 
  
  for (i in 1:NG){
    for (j in 1:NG){
      valid = which(combined_mcc>=grid_mcc[i] & combined_lift>=grid_lift[j])
      if (length(valid)>0){
        selected_label = combined_label[valid]
        fdr = sum(selected_label=="Decoy")/length(valid)
        fdr_matrix[i,j] = fdr
        valid_matrix[i,j]  = sum(selected_label=="Target")
      }
    }
  }
  
  if (min(fdr_matrix)<=fdr_thr){
    index_fdr = which(fdr_matrix<=fdr_thr, arr.ind=T) # Index that have a fdr smaller than fdr_thr
    LS = valid_matrix[index_fdr] # Nb of rules of all minimal fdrs
    max_size = which(LS == max(LS)) # Maximal number of rules that allow minimal fdr
    grid_fdr = index_fdr[max_size[length(max_size)],] # Apply highest possibe threshold to maximize
  
    mcc_min = grid_mcc[grid_fdr[1]]
    lift_min = grid_lift[grid_fdr[2]]
    fdr = min(fdr_matrix)
  } else { # No filter applied if FDR optimization not possible
    max_size =  nrow(rules)
    mcc_min = 0
    lift_min = 1
    fdr = -1 # Failed!
  }
  return(list(max_size = max_size, mcc_min = mcc_min, lift_min = lift_min, fdr = fdr))
}

rule_mass_match<-function(rule_types, rule_features, masslist, ppm_search){
  
  # The function match back comobined rules to raw masslist
  
  rule_types = strsplit(rule_types,",")[[1]]
  rule_features = as.numeric(strsplit(rule_features,",")[[1]])
  
  # Match to raw masslist:
  
  index1 = c() # Index of matched mass in raw masslist
  rule_fragments = rule_features[rule_types=="mass"]
  if (length(rule_fragments)>0){
    for (fragment in rule_fragments){
      errors = ppm_calc(masslist, fragment)
      valid = which(errors<=ppm_search)
      if (length(valid)>0){index1 = c(index1, valid[1])}
  }}
  
  # Match to raw mass diff:
  
  index2 = list() # List of matched mass difference (from...to...)
  rule_mdiff = rule_features[rule_types=="mass_diff"]
  rule_mdiff = unique(rule_mdiff)
  dist_mass = data.matrix(dist(masslist))
  k = 0
  
  if (length(rule_mdiff)>0){
    for (mdiff in rule_mdiff){
      errors = abs((dist_mass-mdiff)/mdiff*1000000)
      valid = which(errors<ppm_search, arr.ind = T)
      if (nrow(valid)>0){
       k = k+1 
       index2[[k]] = as.numeric(valid[1,])
      }}
  }
  
  return(list(index1=index1, index2=index2))
}
  
  
  

options(stringsAsFactors = F)
options(warn=-1)
options(shiny.maxRequestSize=60*1024^2)

library(shiny)
library("V8")
library(shinyjs)
#library(MSnbase)
library(formattable)
library(stringr)
require(DT, quietly = TRUE)
library(prozor)

load("rule_db_multiple_sub_raw.RData")
source('helper.r')
colnames(decoy) = toupper(colnames(decoy))
colnames(rules) = toupper(colnames(rules))

#library(rsconnect)
#deployApp(server="shinyapps.io",appName="MESSAR")

shinyServer(function(input, output,clientData, session) {
  
  observeEvent(input$killButton,{
    shinyjs::js$refresh()})
  
  check_input <-eventReactive(input$goButton,{
  
    masslist = NULL
    massdiff = NULL
    intlist = NULL # Peak intensities
    valid = 1
    mms = ""
    
    inFile1=input$file1
    inFile2=input$file2
    
    if (input$blank_file1==""){
      mms="Please paste your mass peaks!"
      valid = 0}
    
    if (input$prec_mz!=""){  
      if (is.na(as.numeric(input$prec_mz))){
        mms="Precursor mass must be a numeric value!"
        valid = 0} else {
          prec_mz = as.numeric(input$prec_mz)}}
    
    if (input$blank_file1!=""){
      input_str = input$blank_file1
      input_str = strsplit(input_str,"\n")[[1]]
  
      input_str = lapply(input_str, function(x) strsplit(x, " ")[[1]])
      input_str = lapply(input_str, function(x) x[x!="\t"])

      if (all(sapply(input_str,length)==1)){ # One column situation
        masslist=as.numeric(unlist(input_str))
        masslist=masslist[!is.na(masslist)]
        intlist = rep(100, length(masslist))}
      
      if (all(sapply(input_str,length)==2)){ # Two column situation
        masslist = as.numeric(sapply(input_str,function(x) x[1]))
        intlist = as.numeric(sapply(input_str,function(x) x[2]))
        valid_peaks = which(!is.na(masslist) & !is.na(intlist))
        masslist = masslist[valid_peaks]
        intlist = intlist[valid_peaks]}
      
      if (length(masslist)==0){
        mms = "Input mass format not valid! It must be either a single column mass values or two columns <Mass, intensity> separated by "
        valid = 0}}
    
    if (!is.null(masslist) && valid==1){ 
      
      # Add precursor mass
      
      if (input$prec_mz!=""){
      
      valid1 = which(masslist<(prec_mz+1))
      masslist = masslist[valid1]
      intlist  = intlist[valid1]
      ppm_prec = min(ppm_calc(masslist, prec_mz))
      
      if (ppm_prec > input$ppm_search){
        masslist= c(masslist, prec_mz) # Add precursor peak by exact mass
        intlist = c(intlist, max(intlist)*1.1)
      }}
      
      # Filter intensity
      
      kept = which(intlist/max(intlist)>=input$Relative/100)
      intlist = intlist[kept]  
      masslist = masslist[kept]}
    
    if (input$blank_file2!=""){
      input_str=input$blank_file2
      massdiff=as.numeric(strsplit(input_str,"\n")[[1]])
      massdiff=massdiff[!is.na(massdiff)]
      if (length(massdiff)==0){
        mms = "Input mass difference format not valid! It must be a single column numeric vector!"
        valid = 0}}
  
    if (!is.null(masslist) && is.null(massdiff) && valid==1){
      massdiff = unique(as.numeric(dist(masslist)))
      massdiff = massdiff[massdiff>=30]}

    list(masslist=masslist, intlist=intlist, massdiff=massdiff, message=mms,valid=valid)
  })
  
observeEvent(input$exampleButton, {
  
  fileText <- paste(readLines("https://raw.githubusercontent.com/daniellyz/MESSAR/master/example_casmi_2017_203.txt"), collapse = "\n")
  updateTextAreaInput(session, "blank_file1", value = fileText)
})  
  
find_rules <- eventReactive(input$goButton,{

    rules_extracted = NULL
    decoy_extracted = NULL
    selected_index = NULL
    mms = ""
    
    if (check_input()$valid==1){

      masslist=check_input()$masslist
      massdiff=check_input()$massdiff

      # Target rules
      selected_index = search_rules(rules_feature, rules_type, masslist, massdiff, input$ppm_search)
      
      if (length(selected_index)>0){
        
        rules_extracted = rules[selected_index,]
        if (length(selected_index)==1){
          rules_extracted = data.frame(matrix(rules_extracted,nrow=1))
          colnames(rules_extracted) = colnames(rules)}
        
        rules_extracted$MCC = round(rules_extracted$MCC, 2) # Round to improve display
        rules_extracted$CONFIDENCE = round(rules_extracted$CONFIDENCE, 2) # Round to improve display
        rules_extracted$LIFT = round(rules_extracted$LIFT, 2)

        if (input$fdr_control){
        decoy_index = search_rules(decoy_feature, decoy_type, masslist, massdiff, input$ppm_search*2)
        if (length(decoy_index)>0){
          decoy_extracted = decoy[decoy_index,]
          if (length(decoy_index)==1){
            decoy_extracted = data.frame(matrix(decoy_extracted,nrow=1))
            colnames(decoy_extracted) = colnames(decoy)
          }
        }}
  
      mms = "Annotation succeeded! Please check panel B) then C)!"}
        
      if (length(selected_index)==0){
      mms = "No substructure recommended!"} 
    }
    
    list(selected_index=selected_index,rules_extracted=rules_extracted,
        decoy_extracted = decoy_extracted, message=mms)
})
  
observeEvent(input$goButton,{

  withProgress({
    setProgress(message="Check data format...")
    Sys.sleep(1)
    setProgress(message=check_input()$message)
    if (check_input()$valid==1){
      setProgress(message="Annotating substructures by target and decoy rule databases...")
      Sys.sleep(1)
      setProgress(message=find_rules()$message)
    }
  })

  if (check_input()$valid==0){
    updateActionButton(session, "goButton",label = "Try again")
    output$blank_message1<-renderText({check_input()$message})
  } else {
    output$blank_message1<-renderText({find_rules()$message})
  }
})

process_rules <- reactive({
  
  # Depend on fdr threshold, filter rules extracted and calculate aggregated rules
  
  rules_extracted = find_rules()$rules_extracted
  decoy_extracted = find_rules()$decoy_extracted
  
  NR0 = nrow(rules_extracted)
  NS0 = length(unique(rules_extracted$SUBSTRUCTURE))
  ms = ""
  thr1 = 0
  thr2 = 0
  
  if (!is.null(rules_extracted) & !is.null(decoy_extracted)){ 
  
  if ((nrow(rules_extracted)>0) & (input$fdr_control)){
    
      fdr_summmary = fdr_compute(rules_extracted, decoy_extracted, 0.05)
      
      thr1 = fdr_summmary$mcc_min
      thr2 = fdr_summmary$lift_min
      
      valid=which(rules_extracted$MCC>=thr1 & rules_extracted$LIFT>=thr2)
      if (length(valid)>=1){
        rules_extracted = rules_extracted[valid,]}
      
  NR = nrow(rules_extracted)
  NS = length(unique(rules_extracted$SUBSTRUCTURE))
  
  if (fdr_summmary$fdr!=-1){
   ms1= paste0(NR0, " rules and ", NS0, " substructures before")
   ms2= paste0(NR, " rules and ", NS, " substructures kept after FDR control")
   ms = paste("", ms1, ms2, sep="<br/>")} else {
   ms = "No FDR filtering Possible! All rules and substructures are kept!"}}}
  
  rules_aggregated = aggregate_rules(rules_extracted)
  
  output$blank_message2 = renderUI({
    HTML(ms)
  })
  list(rules_extracted=rules_extracted, rules_aggregated = rules_aggregated, thr1 = thr1, thr2 = thr2)
})

output$plot_fdr <- renderPlot({
  
  # Plot raw fdr distribution:
  if (input$fdr_control){
  
  rules_extracted = find_rules()$rules_extracted
  decoy_extracted = find_rules()$decoy_extracted
  thr1 = process_rules()$thr1
  thr2 = process_rules()$thr2
  
  if (!is.null(rules_extracted) & !is.null(decoy_extracted)){
  if (nrow(rules_extracted)>0){
  
  par(mfrow=c(1,3))
  
  TMP1 <- hist(rules_extracted$MCC, plot= FALSE)$counts
  TMP2 <- hist(decoy_extracted$MCC, plot= FALSE)$counts
  ylim = c(0, round(max(c(TMP1, TMP2))*1.05))
  xlim1 = c(0, round(max(c(rules_extracted$MCC, decoy_extracted$MCC))*1.05,2))
  hist(rules_extracted$MCC, col=rgb(0,0,1,0.5),xlim=xlim1, ylim=ylim, xlab="MCC", ylab="Number of rules", main = "", font=2, cex.lab=1.7)
  hist(decoy_extracted$MCC, col=rgb(1,0,0,0.5), add=T)
  
  abline(v=thr1, col="red", lwd=2.5)
  
  legend("topright", c("Target", "Decoy"), col=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), lwd=10, text.font=2, cex=1.5)
  box()
  
  TMP1 <- hist(rules_extracted$LIFT, plot= FALSE)$counts
  TMP2 <- hist(decoy_extracted$LIFT, plot= FALSE)$counts
  ylim = c(0, round(max(c(TMP1, TMP2))*1.05))
  xlim2 = c(1, round(max(c(rules_extracted$LIFT, decoy_extracted$LIFT))*1.05,2))
    hist(rules_extracted$LIFT, col=rgb(0,0,1,0.5),xlim=xlim2, ylim=ylim, xlab="LIFT", ylab="Number of rules", main = "", font=2, cex.lab=1.7)
  hist(decoy_extracted$LIFT, col=rgb(1,0,0,0.5), add=T)
  abline(v=thr2, col="red", lwd=2.5)
  
  legend("topright", c("Target", "Decoy"), col=c( rgb(0,0,1,0.5), rgb(1,0,0,0.5)), lwd=10, text.font=2, cex=1.5)
  box()
  
  plot(rules_extracted$MCC, rules_extracted$LIFT, xlim=xlim1, ylim=xlim2, pch=20, col=rgb(0,0,1,0.5), xlab = "MCC", ylab= "LIFT", main="", font=2, cex.lab=1.7)
  points(decoy_extracted$MCC, decoy_extracted$LIFT, pch=20, col=rgb(1,0,0,0.5))
  abline(v=thr1, h = thr2, col="red", lwd=2.5)
  legend("topright", c("Target", "Decoy"), col=c( rgb(0,0,1,0.5), rgb(1,0,0,0.5)), lwd=10, text.font=2, cex=1.5)
  box()
  }}}
})


output$table1 <- renderDataTable({
  
  rule_table=NULL

  withProgress({
    setProgress(message="Generating annotation results...")
    Sys.sleep(1)
    
    columns = c('SPECTRAL_FEATURE_TYPE', 'SPECTRAL_FEATURE', 'SUBSTRUCTURE', 'LIFT', 'MCC')
    new_columns = c("TYPE","FEATURE","SUBSTRUCTURE", "LIFT", "MCC")
    
    rules_extracted = process_rules()$rules_extracted

    if (nrow(rules_extracted)>0){

      rules_extracted = rules_extracted[,columns]

      if (nrow(rules_extracted)==1){rules_extracted=matrix(rules_extracted, 1, dimnames=list(NULL,new_columns))
      } else {colnames(rules_extracted)=new_columns}
    
    rule_table=datatable(rules_extracted,escape=c(TRUE,TRUE,TRUE,TRUE,TRUE), rownames = F)}
  
    return(rule_table)
    })
  })

RER_generator <-reactive({
  
  rules_aggregated = process_rules()$rules_aggregated
  RER = NULL
  
  if (!is.null(rules_aggregated)){
    
    if (input$score_type=="L1"){RER = eval_rules_aggregated(rules_aggregated, "SUM", "LIFT", 0)}
    if (input$score_type=="L2"){RER = eval_rules_aggregated(rules_aggregated, "MEDIAN", "LIFT", 0)}
    if (input$score_type=="M1"){RER = eval_rules_aggregated(rules_aggregated, "SUM", "MCC", 0)}
    if (input$score_type=="M2"){RER = eval_rules_aggregated(rules_aggregated, "MEDIAN", "MCC", 0)}
    if (input$score_type=="C1"){RER = eval_rules_aggregated(rules_aggregated, "SUM", "CONFIDENCE", 0)}
    if (input$score_type=="C2"){RER = eval_rules_aggregated(rules_aggregated, "MEDIAN", "CONFIDENCE", 0)}}
  
  return(RER)
})


output$table2 <- renderDataTable({
  
  rule_table=NULL
  
  withProgress({
    
    setProgress(message="Generating annotation results...")
    Sys.sleep(1)
    RER = RER_generator()
    
    # Add image:
    
    setProgress(message="Generating substructure visualization...")
    id_substructure = match(RER$SUBSTRUCTURE, substructure_db$substructure)
    id_substructure = substructure_db$ID[id_substructure]
    RER$IMG = paste0('<img src=','"', id_substructure,'.png" height="150"></img>') 
    RER = RER[,c(1,3,2)]

    if (nrow(RER)>0){
      if (nrow(RER)==1){
        RER=matrix(RER, 1 , dimnames=list(NULL,colnames(RER)))} 
    
    rule_table=datatable(RER,escape=c(TRUE,FALSE,TRUE), rownames = F,  selection = "single", options = list(pageLength=5))}
    
    return(rule_table)
    })
  })

selected_substructures <- eventReactive(input$table2_rows_selected,{
  RER_generator()$SUBSTRUCTURE[c(input$table2_rows_selected)]
})

output$plot_selected <- renderPlot({

  RER = RER_generator()
  selected = which(RER$SUBSTRUCTURE == selected_substructures())[1]
  
  # Plot distributions:
  par(mfrow = c(2,1))
  TMP1 <- hist(RER[,2], plot= FALSE)$counts
  
  hist(RER[,2], col="blue", xlab="Score", ylab="Frequency", main = "", font=2, cex.lab=1.7, font.lab=2)
  arrows(RER[selected,2],0, RER[selected,2], max(TMP1)*0.8,col="red", lwd = 2)
  box()
  
  raw_mass_list = check_input()$masslist
  raw_int_list = sqrt(check_input()$intlist)
  orders = order(raw_mass_list, decreasing = F)
  raw_mass_list= raw_mass_list[orders]
  raw_int_list = raw_int_list[orders]
  
  # Plot features behind each substructure:
  
  rules_aggregated = process_rules()$rules_aggregated
  selected = which(rules_aggregated$sum_aggregated$SUBSTRUCTURE==selected_substructures())[1]
  rule_types = rules_aggregated$type[selected] # Features behind the selected rule
  rule_features = rules_aggregated$features[selected]
  matches = rule_mass_match(rule_types, rule_features, raw_mass_list, input$ppm_search)
  
  # Plot fragments:
  collist = rep("grey", length(raw_mass_list))
  if (length(matches$index1)>0){
    collist[matches$index1] = "red"
  }
  ylim = c(0, max(raw_int_list)*1.1)
  plot(raw_mass_list,raw_int_list, type = "h", xlab = "m/z", ylab = "Intensity", ylim = ylim,
       font=2, font.lab=2, cex.lab=1.7, col=collist, lwd = 2)
  if (length(matches$index1)>0){
   text(raw_mass_list[matches$index1], raw_int_list[matches$index1]*1.1, round(raw_mass_list[matches$index1],2), col="blue", cex = 0.9)
  }
  
  # Plot mass differences:
  if (length(matches$index2)>0){
    for (i in 1:length(matches$index2)){
     frag_from = raw_mass_list[matches$index2[[i]][1]]
     frag_to = raw_mass_list[matches$index2[[i]][2]]
     int_base = raw_int_list[matches$index2[[i]][1]]/2
     arrows(frag_from,int_base, frag_to, int_base, col="red", lwd = 2)
     
     frags = unique(c(frag_from, frag_to))
     ints = raw_int_list[match(frags, raw_mass_list)]
     text(frags, ints*1.1, round(frags,2), col="blue", cex = 0.9)
     
    }
  }
})

})
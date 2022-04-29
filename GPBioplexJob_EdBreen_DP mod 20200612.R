
# zipFile
# DataofInterest FI, "FI,Obs Conc"
# CreatePID "no"
# adjustment "none" # fdr, holm
# exclude
# normaliseRoutine "none" "row,col"
# Kit Stats "(1|ID)"

# setwd("C:\\Users\\Dana\\Documents\\0-WORK\\Bioplex")

# setwd("V:/Projects/Bioplex/P31251 - Fatima Nasrallah")


args = c("zplates.zip","Design-remove-missing-samples.xlsx","DFI")

# bioplexGP("zJune27_2017.zip", "MDesign.xlsx", "DFI", "mDescription")
# bioplexGP("zRESMED.zip", "DFI", "pYes", "MDesign.xlsx", 'mpid:Type', "anone", "eDescription", "cNo")
# bioplexGP("zplates.zip","MDesign-remove-missing-samples.xlsx","DFI")
# bioplexGP("zplates.zip","MDesign-remove-missing-samples.xlsx","DFI, Obs Conc", "afdr", "Nrow,col")
# bioplexGP("zplates.zip","MDesign-remove-missing-samples.xlsx","DFI, Obs Conc", "afdr", "Nkit:row,kit:col")

bioplexGP<-function(...){
  
  # options( scipen = 0 )
  # options( digits = 3 )
  
  
  bplxInstalls()     
  trim <- function (x) gsub("^\\s+|\\s+$", "", x) # trim white space from ends
  
  makeIDfromBaseName = function(filename, ...)
  {
    sprintf("p%03d", as.numeric(gsub("[^0-9]", '', basename(filename))));
  }
  stars.pval <- function(p.value)
  {
    unclass(
      symnum(p.value, corr = FALSE, na = FALSE,
             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
             symbols = c("***", "**", "*", ".", " "))
    )
  }
  addTable=function(wb, sheetName, df, hstyle, tstyle,rowNames=F, 
                    withFilter=F,
                    name=NULL,
                    tableStyle = "TableStyleLight1") 
  {
    addWorksheet(wb, sheetName)
    if(!is.null(name))
      writeData(wb, sheetName, name, #headerStyle=hstyle, 
                     rowNames=rowNames)
    
    startrow = 2
    writeDataTable(wb, sheetName, as.data.frame(df), headerStyle=hstyle, 
              rowNames=rowNames, withFilter=withFilter,
              tableStyle=tableStyle)
    
    startrow + nrow(as.data.frame(df))
    
  }
  
  addBoxPlot = function(wb,sheet, name, form, data, ylab, xlab, main, las = 1, 
                        startCol=1, startRow = 1, width = 6, height = 5, dpi =600,
                        mar = c(4,4,3,1))
  {
    png(name, res=dpi, units = "in", width=width, height=height)
    par(mar=mar)
    boxplot(form, data, main=main, xlab=xlab, ylab=ylab, main=main, las = las)
    dev.off()
    insertImage(wb, sheet = sheet,
                file = name,
                startCol = startCol,
                startRow = startRow, 
                width = width, height = height,
                dpi=dpi)
  }
  
  print("BioPlex APAF Version 1.0.0")
  
  hstyle <- createStyle(fontSize = 11, fontName = "Arial",
                       textDecoration = "bold", 
                       halign = "center" #fgFill="gray95", 
                       #fontColour ="black"
                       )
  
  tstyle <- createStyle(halign = "center", numFmt = "0.000")
  
  
   
  args<-list(...)
  metaFile = NA
  print(args)
  exclude = NA
  
  doplateID = 'No'
  adjustment = 'none'
  normalise = 'none'
  withinGrp = NA
  for(i in 1:length(args)){
    flag<- substring(args[[i]],0,1) ## -x
    value<-substring(args[[i]],2, nchar(args[[i]]))
    
    if(flag == "z") zipFile <- value
    if(flag == "D")  Sheetvalue <-trim(unlist(strsplit(value,',')))  # a single addon lib file or a zip file
    if(flag == 'M') metaFile <- value
   
    if(flag == 'p') doplateID = value
    if(flag == 'a') adjustment = value
    if(flag == 'e') exclude = trim(unlist(strsplit(value,',')))
    if(flag == 'N') normalise = trim(unlist(strsplit(value,',')))
    if(flag == 'S') withinGrp = value
    
    
  }
  
  zipFiles = unzip(zipFile)
  #print(Sheetvalue)
  if(doplateID == 'Yes')
    rmd = bplxGetData(zipFiles,Sheetvalue, exclude = exclude, FUN = makeIDfromBaseName)
  else
    rmd = bplxGetData(zipFiles,Sheetvalue,exclude = exclude)
  
  file="BioPlexAnalaysis.xlsx"
  
  bp= rmd$bp
  
  nbp = NA
  bp = droplevels(bp[bp$grp == 'X' & bp$Type2 == 'raw',])
  if(normalise[1] != 'none') {
    nstr = paste("(1|", normalise, ")", collapse ='+') 
    nbp = bplxCorrectSysmBias(bp, str = nstr)
    bp = nbp$bp
  } else {
    nstr = 'none'
    bp$nFI = bp$FI
  }
  
  if(!is.na(metaFile)) {
    md = read.xlsx(metaFile, sheet = 1)
    md = convertChar2Factor(md)
    
    by = intersect(colnames(bp), colnames(md))
    bp = merge(bp,md, by = by, all.x = T)
    
    dropOn = colnames(md)[1]
    bp = droplevels(bp[ !is.na(bp[,dropOn]),])
  }
  
  wb <- createWorkbook()
  tableStyle="TableStyleLight1"
  
  # write header sheet contain rmd's kits, Readings and Header information
  addTable(wb, "Header", rmd$kits,hstyle,tstyle,withFilter=F)
  startRow = 3 + nrow(rmd$kits)
  writeDataTable(wb, "Header", as.data.frame(rmd$Readings), #headerStyle=hstyle, 
                 tableStyle=tableStyle,
                 startRow =startRow,
                 withFilter=F,
                 rowNames=T)
  startRow = startRow + nrow(rmd$Readings) + 2
  writeDataTable(wb, "Header", rmd$Header, #headerStyle=hstyle, 
                 tableStyle=tableStyle,
                 withFilter=F,
                 startRow =startRow,rowNames = F, tableName ="Yes")
  
  addTable(wb, "Control-Data", rmd$bp[rmd$bp$grp != 'X', ], hstyle, tstyle)
  addTable(wb, "BP-Data", bp, hstyle,tstyle)
  
  
 
  # Kit.Stats
  #str=paste('(1|', withinGrp, ')'
  if(!is.na(withinGrp)) {
    icc = bplxICC(bp, clevel = 'Cytokine', str=withinGrp, DV = 'FI')
    lod = bplxGetLODs(droplevels(rmd$bp[ rmd$bp$Type2 == 'raw',]))
    colnames(lod) = paste("LoD.", colnames(lod), sep='')
    colnames(lod)[1] = "Cytokine"
    startRow = 1
    addWorksheet(wb, sheetName="Q.Stats")
    writeData(wb,"Q.Stats", paste("Within", withinGrp),  startCol =1, startRow=startRow)
    # ct = as.data.frame(VarCorr(nbp$fit))[,-(2:3)]
    # writeDataTable(wb,"Q.Stats", ct,  startCol =2, startRow=startRow+1,
    #                tableStyle=tableStyle)
    # startRow = nrow(ct) + 1
    
    startRow = startRow + 2
    for(ic in names(icc)) {
          
      #merge(lod, icc[[ic]])
      m = merge(icc[[ic]], lod, by = "Cytokine")
      writeData(wb, "Q.Stats", ic, headerStyle=hstyle, startRow=startRow)
      writeDataTable(wb, "Q.Stats", m, headerStyle=hstyle,
                     tableStyle=tableStyle,
                startRow = startRow +1)
      
      startRow<- startRow + nrow(m) + 3
      
      
    }
  }
  # collect individual concentration curves
  # sh = sheetNames(zipFiles[1])
  sh = names(loadWorkbook(zipFiles[1]))
  
  if(length(grep("Curve",sh))>0) {
    df = data.frame()
    for(i in 1:length(zipFiles)) {
      curves = bplxgetCurves(zipFiles[i])
      d = data.frame(File = zipFiles[i],  curves )
      df = rbind(df,d)
      
    }
    
    addTable(wb, "Con. Curves", df,hstyle, tstyle)
    
    df = bplxAllCurves(rmd)
    addTable(wb, "Curves", df,hstyle, tstyle)
  }
  
  
  # do standard plots
  qualitySheet = "Plots Quality"
  
  addWorksheet(wb, qualitySheet)

 
  if(nlevels(factor(bp$kit:bp$pid)) > 10)
        width = 14
  else
        width = 6
  
  startRow = 1
  cexs = 0.22
  addBoxPlot(wb, qualitySheet, "Quality1.png",log2(FI)~sample, data = rmd$bp, 
             startCol = 1, startRow = startRow,
             main = "Dynamic Range", xlab ="Sample Type", ylab="log2(FI)", 
             las=1, width=6, height=5, mar = c(4,4,2,1))
  
   addBoxPlot(wb, qualitySheet, "kitPid.png",log2(FI)~factor(bp$kit:bp$pid), data = bp, 
             startCol = 10,
             main = "kit:pid", xlab ="", ylab="log2(FI)", 
             las=2, 
             width=max(4, cexs * nlevels(factor(bp$kit:bp$pid))), 
             height=5, mar = c(6,4,2,1))
  
  
  startRow = 27
  addBoxPlot(wb, qualitySheet, "Quality2.png",log2(FI)~Cytokine, data = bp, 
             startCol = 1, startRow = startRow,
             main = "Cytokine Levels", xlab ="", ylab="log2(FI)", 
             las=2, 
             width= max(6, cexs * nlevels(bp$Cytokine)), 
             height=5, mar = c(10,4,2,1))
  startRow <-startRow + 24
  addBoxPlot(wb, qualitySheet, "Quality3.png",log2(FI)~row, data = bp, 
             startCol = 1, startRow = startRow,
             main = "Test Sample Rows", xlab ="Row", ylab="log2(FI)", 
             las=2, width=6, height=4, mar = c(4,4,2,1))
  if(normalise[1] != 'none') {
       addBoxPlot(wb, qualitySheet, "Quality4.png",log2(nFI)~row, data = bp, 
             startCol = 10, startRow = startRow,
             main = "Normalised: Rows", xlab ="Row", ylab="log2(nFI)", 
             las=2, width=6, height=4, mar = c(4,4,2,1))
  }
  startRow <-startRow + 22
  addBoxPlot(wb, qualitySheet, "Quality5.png",log2(FI)~col, data = bp, 
             startCol = 1, startRow = startRow,
             main = "Test Sample Columns", xlab ="Column", ylab="log2(FI)", 
             las=1, width=6, height=4, mar = c(4,4,2,1))
  if(normalise[1] != 'none') {
    addBoxPlot(wb, qualitySheet, "Quality6.png",log2(nFI)~col, data = bp, 
             startCol = 10, startRow = startRow,
             main = "Normalised: Columns", xlab ="Column", ylab="log2(nFI)", 
             las=1, width=6, height=4, mar = c(4,4,2,1))
  }
    
  # construct formulas and do Statistics
  
  if(!is.na(metaFile)) { 
    comps = read.xlsx(metaFile, sheet=2)
    d = dim(comps)
    cnt = 2
    for(i in 1:d[1]) {
      DV = "nFI"
      bp[,  comps[i,'Comparison']] = factor(bp[,  comps[i,'Comparison']])
      
      if(is.na(comps[i,'wrt'])) 
        fstr = paste("log2(", DV, ")~", comps[i,'Comparison'], sep='')
      else 
        fstr = paste("log2(", DV, ")~", comps[i,'Comparison'], ':', comps[i,'wrt'],  sep='')
      
      
      pfile = sprintf("Q%03d.png", cnt <- cnt +1)
      
      png(pfile)
      par(mar=c(10,4,2,1))
      boxplot(formula(fstr), bp,
              main = fstr,
              ylab="log2(nFI)", xlab="", las = 2)
      dev.off()
      insertImage(wb, sheet = qualitySheet,
                  file = pfile,
                  startCol = 1,
                  startRow = startRow <- startRow + 25,
                  width = 6, height = 5,
                  dpi=600)
      
      
      #build model
      f = formula(paste("log2(", DV, ")~", comps[i,'Comparison'], '*Cytokine', sep=''))
      str =paste("log2(", DV, ")~", comps[i,'Comparison'], '*Cytokine', sep='');
      
      fixed =  "Cytokine"
      paired = comps[i,'Comparison']
      mixed = F
      
      if(!is.na(comps[i,"wrt"])) {
        minus = paste("- Cytokine", comps[i,'Comparison'], comps[i,"wrt"], sep=':')
        str = paste("log2(", DV, ")~Cytokine *", comps[i,'Comparison'],"*", comps[i,"wrt"], minus, sep='')
        fixed = c(fixed, comps[i,"wrt"])
      }
      
      
      rstr = NULL
      fstr = NULL
      covariates = NULL
      if(!is.na(comps[i,"Covariate"])) {
            covariates = unlist(strsplit(comps[i,"Covariate"],','))
            
            fstr = paste("+", covariates, collapse = '')
            str = paste(str, fstr,sep=' ')
            
            # for(k in covariates) {
            #       if(is.numeric(bp[,k])) 
            #             fstr = c(fstr,k)
            #       else
            #             rstr = c(rstr,k)
            # }   
            # if(!is.null(fstr)) {
            #       fstr = paste("+", fstr, collapse = '')
            #       str = paste(str, fstr,sep=' ')
            # }
            # if(!is.null(rstr)) {
            #       rstr = paste('+(1|', rstr, ')', collapse = '')
            #       str = paste(str, rstr, sep = " ")
            # }
            
      } 
      
      if(!is.na(comps[i,"pairedOn"])) { 
            str = paste(str,"+(1|", comps[i,"pairedOn"], ")", sep='')
            mixed = T
            
      } 
      
      f = formula(str)
      if(mixed)
            fit = lmer(f, bp)
      else
            fit = lm(f, bp)
      
      covariates = paste(covariates, collapse ="")
      sheetName = paste(comps[i,'Comparison'], comps[i,'wrt'], covariates)
      addWorksheet(wb, sheetName)
      
      #setColWidths(wb, sheetName, widths = "auto")
      
      t = testInteractions(fit, fixed = fixed, pairwise= paired, adjust = adjustment)
      t = data.frame(t)
      #t = round(t,4)
      t$Sig = stars.pval(t[,4])
      
      s = strsplit(rownames(t), ':')
      d = lapply(s, trim)
      d1 = data.frame(matrix(unlist(d), nrow = length(d), byrow=T))
      if(dim(d1)[2] == 3)
            colnames(d1) = c("Contrast", "Cytokine", "wrt")
      else
            colnames(d1) = c("Contrast", "Cytokine")
      
      t = cbind(d1,t)
      rownames(t) = NULL
      
      
      a = Anova(fit,type = 3)   
      outFit = capture.output(a)
      mlines = length(outFit)
      writeData(wb, sheetName, outFit[1], startCol =1, startRow=1)
      
      d = data.frame(as.matrix(a))
      d[ d < 1e-100] = 0.00
      # is.num <- sapply(d, is.numeric)
      # d[is.num] <- lapply(d[is.num], formatC, format = "e", digits = 3)
      d$Sig = stars.pval(d$Pr..Chisq.)
      mergeCells(wb, sheetName, cols = 1:(dim(d)[2]+1), rows =1)
      
      writeDataTable(wb,sheetName, d, startCol =1, startRow=2, rowNames=T,
                     tableStyle=tableStyle,
                     withFilter=F)
      
      
      writeData(wb,sheetName, outFit[mlines], startCol =1, startRow=(dim(d)[1]+3))
      mergeCells(wb, sheetName, cols = 1:(dim(d)[2]+1), rows =(dim(d)[1]+3))
      
      start = dim(d)[1] + 6
      writeData(wb,sheetName, "Model",  startCol =1, startRow=start)
      writeData(wb,sheetName, str,  startCol =2, startRow=start)
      mergeCells(wb, sheetName, cols = 2:10, rows =start)
      
      start = start + 2
      writeData(wb,sheetName, "Normalized By",  startCol =1, startRow=start)
      writeData(wb,sheetName, nstr,  startCol =2, startRow=start)
      mergeCells(wb, sheetName, cols = 2:10, rows =start)
      
      start = start + 1
      
      if(normalise[1] != "none") {
            ct = as.data.frame(VarCorr(nbp$fit))[,-(2:3)]
            ct[,2:3] = round(ct[,2:3], 4)
            writeDataTable(wb,sheetName, ct,  startCol =2, startRow=start,
                           tableStyle=tableStyle,
                           withFilter=F)
            start = start + nrow(ct) + 1
      }
      
      # is.num <- sapply(t, is.numeric)
      # t[is.num] <- lapply(t[is.num], formatC, format = "e", digits = 3)
      start = start + 1
      writeData(wb,sheetName, "Multiple test correction",  startCol =1, startRow=start)
      writeData(wb,sheetName, adjustment,  startCol =2, startRow=start)
      
      start = start + 2
      writeData(wb,sheetName, "Comparison",  startCol =1, startRow=start)
      writeDataTable(wb,sheetName, 
                     data.frame(PairWise = paired, Within = paste(fixed, collapse=':')),
                     startCol =2, startRow=start,
                     tableStyle=tableStyle,
                     withFilter=F)
      
      writeDataTable(wb, sheetName, t[order(t$Pr..Chisq.),], rowNames = F, startRow =start + 5,
                     headerStyle = hstyle,
                     tableStyle=tableStyle)
    }
    
  }
  
  
  saveWorkbook(wb, file, overwrite=TRUE)
  
  
}

convertChar2Factor = function(sheet) {
  i <- sapply(sheet, is.character)
  sheet[i] <- lapply(sheet[i], as.factor)
  
  
  sheet
}
convertFactor2Char = function(sheet) {
  i <- sapply(sheet, is.factor)
  sheet[i] <- lapply(sheet[i], as.character)
  
  
  sheet
}

bplxCorrectSysmBias = function(bp, DV="log2(FI)", str="(1|Description) + (1|Analyte) + (1|Analyte:Description)")
{
  bp = droplevels(bp[ bp$Type2 == 'raw',])
   
  f = as.formula(paste(DV, "~", str, sep='' ))
  
  fit = lmer(f, bp)
  #fit = lmer(log2(FI)~(1|pid) +  (1|Well:sample:Analyte), bp)
  
  a = fixef(fit)
  
  bp$nFI = 2^(resid(fit) + a)  
  
   list(bp = bp, fit = fit)
}

bplxGetLODs = function(bp)
{
  x = droplevels(bp[bp$grp == 'B',])
  
  m = aggregate(FI~Cytokine, data=x, mean, na.rm=T)
  s = aggregate(FI~Cytokine, data=x, sd, na.rm=T)
  l = aggregate(FI~Cytokine, data=x, length)
  
  data.frame(Analyte = m$Cytokine, N = l[,2], mean = m$FI, sd = s$FI, LOD = m$FI+2*s$FI)
  
  
}

bplxAllCurves = function(bpl)
{
  if(!is.list(bpl))
    stop("Require list input")
  #analytes = bpl$kits
  
  df = NULL
  for(i in 1:nrow(bpl$Header)) {
    
    kit = bpl$Header$kit[i]
    curves = bplxgetCurves(as.character(bpl$Header$file[i]))
    # a = unlist(strsplit(as.character(bpl$kits$Analytes[ bpl$kits$kit == kit]), ','))
    # curves$Analyte = unique(a)
    #curves$pid = bpl$Header$pid[i]
    
    df = rbind(df,curves)
    
  }
  a = aggregate(df[, c('d','a', 'c', 'b', 'g')], list(Analyte = df$Analyte), FUN=mean)
  a
  
}

bplxgetCurves = function(file){
  
  curvefmt = function(x) {
    ss = unlist(strsplit(x, " ")[[1]])
    n = 7 - length(ss)
    if(n > 0) ss = c(ss, rep(NA,n))
    ss
  }
  
  bplxInstalls('openxlsx')
  wb = loadWorkbook(file)
  
  sheet = read.xlsx(wb, sheet="STD Curve", colNames=F)
  
  cols = 1:dim(sheet)[2]
  t = c(t(sheet[,cols]))
  Type = t[ grep("Regression", t)]
  Curve = t[ grep("Std[\\.]? Curve", t)]
  Fit = t[ grep("FitProb.", t)]
  #s = gsub("[^[:digit:].]"," ", x)
  
  Type = unlist(lapply(Type, function(x) strsplit(as.character(x),": ")[[1]][2]))
  Curve =unlist(lapply(Curve, function(x) strsplit(as.character(x)," = ")[[1]][2]))
  Fit1 = unlist(lapply(Fit, function(x) strsplit(as.character(x),". = ")[[1]][2]))
  Fit2 = unlist(lapply(Fit1, function(x) strsplit(as.character(x),",")[[1]][1]))
  FitProb = as.numeric(Fit2)
  ResVar = as.numeric(unlist(
    lapply(Fit, function(x) strsplit(as.character(x),". = ")[[1]][3])))
  
  #Curve Type
  # watch for 4PL and and 5PL 
  
  sheet = read.xlsx(wb, sheet=1, colNames=F)
  row = grep('Type', sheet[,1])[1] - 1
  analytes = sheet[row,!is.na(sheet[row,])]
  
  ss = gsub("[- |+ ] ", "", gsub("x+", " ", gsub("[^0-9+.Ee-]", "x", Curve)))
  df  = data.frame(lapply(ss, FUN=curvefmt))
  colnames(df) = NULL;
  
  
  
  df = t(df)[,c(1:2, 5:7)]
  if(is.null(dim(df))) {
        df = data.frame(t(as.numeric(df)))
  } else {
      nr = nrow(df)
      df = matrix(as.numeric(df), nrow=nr, byrow=F)
  }
  colnames(df) = c('d', 'a', 'c', 'b', 'g')
  
  # get Analytes names 
  

  
  df = data.frame(Analyte = t(analytes)[,1], Type,df, FitProb , ResVar)
  rownames(df) = NULL  
  df
}

bplxColors = function(bp=NULL) {
  
  if(is.null(bp)) {
    bp = data.frame( grp = c('C', 'S', 'B', 'X') ,
                     samples = paste('S', 1:8,sep=''),
                     kit = 'k001')
  }
  
  mypal = c("grey80", "grey40", "hotpink1", "blue3", "tomato1", 
            "purple1","cyan1")
  
  
  num.stds = length(unique(grep('S', levels(bp$sample))))
  bplxColor = list()
  bplxColor$plate = "palegreen"
  
  
  if(!is.null(bp$kit))
    bplxColor$kit = mypal[1:(nlevels(bp$kit))]
  
  
  std.cols = heat.colors(num.stds)
  if(num.stds > 0) # assume decreasing
    std.cols = std.cols[num.stds:1]
  
  exp = c('C', 'S', 'B', 'X')
  id = match(exp, levels(bp$grp))
  
  std.col = NULL
  if(!is.na(id[2])) std.col = std.cols[as.integer(length(std.cols)/2)]
  
  
  blnk.col = NULL
  if(!is.na(id[3]))
    blnk.col =  "forestgreen"
  
  tst.col = NULL
  if(!is.na(id[4]))
    tst.col =  "steelblue"
  
  cntrl.cols = NULL
  if(!is.na(id[1])) cntrl.cols = c("turquoise1", "turquoise3")
  if(is.na(id[2])) std.cols = NULL
  
  
  bplxColor$std.col = std.col
  bplxColor$standards = std.cols
  bplxColor$blank = blnk.col
  bplxColor$test = tst.col
  bplxColor$control = cntrl.cols
  bplxColor$group = c( cntrl.cols[1],std.col,
                       blnk.col,tst.col)
  
  sample.cols =  c(cntrl.cols, std.cols,   blnk.col, tst.col)
  bplxColor$sample = factor(sample.cols, levels = sample.cols)
  
  
  return (convertFactor2Char(bplxColor))
}

bplxInstalls =function(pkgs = c('openxlsx','phia', 'lme4', 'lmerTest', 'nnet', 
                                'gdata', 'lattice', 'latticeExtra', 'reshape', 'visreg', 
                                'gplots', 'e1071', 'drc', 'magic')) {
  dynamic.load = function(package)
  {
    #,'sjmisc', 'sjPlot'
    if(eval(parse(text=paste("require(",package,")")))) return(T)
    install.packages(package,repos = "http://cran.us.r-project.org")
    return (eval(parse(text=paste("require(",package,")"))))
  }
  
  
  for(i in pkgs) 
    dynamic.load(i)
}

bplxFactorNos = function(bp, factors, rtn=F) {
  df = data.frame(table(bp[, factors]))
  df = droplevels(df[ df$Freq != 0,])
  if(!rtn) {
    nms = colnames(df)
    df = bplxFactorNos(df, nms[2: (length(nms)-1)], rtn=T) 
    colnames(df)[1] = nms[2]
  }
  df
}
bpxlAnalyteOrder = function(files, sheetName="FI")
{
  
}
bplxGetAnalytes = function(file,sheetName="FI", toNumeric=T){
  
  
  if(length(file) > 1) stop("Input file must be singular entity")
  bplxInstalls('openxlsx')
  bplxInstalls('gdata')
 
  sheet = read.xlsx(file, sheet=sheetName, colNames =F)
  
  
  xp = grep("^X[0-9]+$", sheet[,1])
  bp = grep("^B[0-9]*$", sheet[,1])
  cp = grep("^C[0-9]+$", sheet[,1])
  sp = grep("^S[0-9]+$", sheet[,1])
  
  if(length(sp)) {
    # catch non-standard standard names
    y = gsub("[^S0-9]","", sheet[sp,1])
    sheet[sp,1] = y
    
  }
  
  
  
  ap = sort(c(bp,sp,cp,xp))
  
  n = sapply(1:ap[1], function(x) {   length(grep(sheetName, sheet[x,], fixed=T))})
  a = which(n == max(n))
  cols = grep(sheetName, sheet[a,])
  end.c = cols[1]-1
  
  #     an = length(ap)
  #     cn = length(cols)
  x = sheet[ap, c(1:end.c,cols)]
  
  # create column header
  
  h1 = as.character(sheet[ a, 1:end.c])
  h2 = as.character(sheet[ a-1, cols])
  
  colnames(x) = append(h1,h2)
  if(toNumeric)
    suppressWarnings(ifelse(length(cols) == 1, 
                            x[,h2] <- sapply(x[,h2], as.numeric),
                            x[,h2] <- lapply(x[,h2], as.numeric)))
  
  # get Type variable components
  t = x$Type
  cb = grep('B', t)
  cs = grep('S', t)
  cc = grep('C', t)
  cx = grep('X', t)
  
  
  # order tests and standards
  xod = order( suppressWarnings(as.numeric(gsub("[A-Z]", "", t[cx]))))
  sod = order( suppressWarnings(as.numeric(gsub("[A-Z]", "", t[cs]))))
  idx = c(cc, cs[sod], cb, cx[xod])
  x$Type = factor(t,  levels=unique(t[idx]))
  
  
  # create type2 variable
  type2 = rep("raw", nrow(x))
  n = sapply(x$Well, function(x) {nchar(x)})
  id = n > 4
  type2[id] = "mfi"
  
  # identify and create replicate data
  rip = bplxReplicates(x$Type)

  # order Well variable
  x$Well = factor(x$Well, levels=unique(x$Well[idx]))
  
  #make grp
  grp = t
  grp[cx] = 'X'
  grp[cb] = 'B'
  grp[cs] = 'S'
  grp[cc] = 'C'
  
  # make group variable
  grp = factor(grp, levels = unique(c(grp[cc], grp[cs],grp[cb],grp[cx])))
  
  # make sample variable
  t[cx] = 'X'
  sample = factor(t, levels=unique(c(t[cc], t[cs][sod], t[cb], "X")))
  
  # create and order column variable Well
  col = gsub("[[:alpha:]]","",x$Well)
  cod = order(as.numeric(gsub(",",'0', col)))
  col = factor(col, levels = unique(col[cod]))
  
  # create and order row variable
  row =gsub("[[:digit:]]","",x$Well)
  od = data.frame( as.character(row),nchar(row), stringsAsFactors = F)
  od = od [ order( od[,2], od[,1]),]
  row = factor(row, levels = unique(od)[,1])
  
  # construct return object
  x = cbind(sample= sample,  grp = grp, row=row , Rep = rip$Rep,
            col = col, Type2 = type2, x)
  unik = !duplicated(x[ , c("Type", "Well")])
  
  
  
  return(x[ unik, ])
  
}


bplxGetData = function(files, sheetName="FI", 
                       exclude = NA,
                       FUN = makePlateId, ...){
  # nested functions
  makePlateId = function(filename, cnt) sprintf("p%03d", cnt)
  stackIt = function(d, description = "Description", parName="Cytokine", parVal="Fl"){
    c = which(colnames(d) == description) 
    nms = colnames(d)[1:c]
    ncyt = ncol(d) - c
    if(ncyt == 0)
          return (NULL)
    if(ncyt == 1) {
      s = data.frame(d[, c+1], rep(colnames(d)[c+1], nrow(d)))    
    } else {
      s = stack(d[,(c+1):ncol(d)])  
    }
    colnames(s) = c(parVal, parName)  
    if(length(nms) == 1) {
      x = as.data.frame(rep(d[,nms], ncyt))
      colnames(x)[1] = nms
    } else
      x = as.data.frame(sapply(d[,nms], rep, ncyt, simplify=F))
    x$plex = ncyt
    return(cbind(x,s[, 2:1]))
    
  }
  getBplxHeader = function(file, sheetName="FI"){
    bplxInstalls('openxlsx')
    rtn=data.frame()
    #analytes = NULL
    #    for(i in files) {
    #       print(i)
    
    lst = c("Date: ","Number: ","Plate ID: ","RP1 PMT \\(Volts\\): ", "RP1 Target: ")
    
    sheet = read.xlsx(file, sheet=sheetName,colNames=F)
    
    date = NULL;time = NULL;date=NULL;Reader=NULL;Plate =NULL;Volts=NULL;Target=NULL;
    for(j in 1:length(lst)) {
      sp = grep(lst[j], sheet[,1])[1]
      if(length(sp)) {
        wd = strsplit(as.character(sheet[sp,1]), lst[j])[[1]][2]
        switch(j,  
               {date = as.Date(substr(wd,1,11),format = "%d-%b-%Y"); time = substr(wd, 14, 21)},
               {Reader = wd},
               {Plate = wd},
               {Volts = wd},
               {Target=wd})
      }
    }   
    
    
    id = grep("Type", sheet[,1])
    cols = ""
    if(length(id)) {
      cols = sheet[id[1]-1,] 
      cols = paste(cols[!is.na(cols)], collapse='\t')
      
    }
    
    rtn = data.frame(Reader = Reader, 
                     Read.Date = date,
                     Read.Time = time,
                     Plate = Plate,
                     Volts = Volts,
                     Target = Target)
    
    return (rtn)
    
  }  
  bplxOrderFactors = function(bp) {
    # get Type variable components
    type = as.character(levels(bp$Type))
    cb = grep('B', type)
    cs = grep('S', type)
    cc = grep('C', type)
    cx = grep('X', type)
    
    
    # order tests and standards
    xod = order( suppressWarnings(as.numeric(gsub("[A-Z]", "", type[cx]))))
    sod = order( suppressWarnings(as.numeric(gsub("[A-Z]", "", type[cs]))))
    idx = c(cc, cs[sod], cb, cx[xod])
    bp$Type = factor(bp$Type,  levels=type[idx])
    
    
    y = data.frame(type = bp[,"Type"], well = bp[, "Well"])
    bp$Well = factor(bp$Well, levels = unique(y[order(y[,1]),2]))
    
    #make grp
    grp = type
    grp[cx] = 'X'
    grp[cb] = 'B'
    grp[cs] = 'S'
    grp[cc] = 'C'
    
    # make group variable
    bp$grp = bp$Type
    levels(bp$grp) = c(grp[cc], grp[cs],grp[cb],grp[cx])
    #bp$grp = factor(grp, levels = unique(c(grp[cc], grp[cs],grp[cb],grp[cx])))
    
    # make sample variable
    type = levels(bp$Type)
    type[grep('X', type)] = 'X'
    bp$sample = bp$Type
    levels(bp$sample) = type
    #bp$sample = factor(type, levels=unique(c(type[cc], type[cs][sod], type[cb], "X")))
    
    # create and order column variable Well
    col = gsub("[[:alpha:]]","",bp$Well)
    cod = order(as.numeric(gsub(",",'0', col)))
    bp$col = factor(col, levels = unique(col[cod]))
    
    bp  
    
    
  }
  processbplxInput = function(bpl, FUN=NULL, ...) {
    if(!is.list(bpl))
      stop("Require list input")
    nCytNme = function(cytokines, remleader=F) {
      # normalize cytokine names of the form  'XX name(region.No)'
      #for example: cytokines = levels(bp$Cytokine)
      if(remleader) cytokines = gsub("^.. ","", cytokines)
      gsub(" \\([0-9]+\\)$","", cytokines)
    }
    bplxAssignCytSet = function(bpl) {
      if(!is.list(bpl)) stop("Require bpl to be a list!")
      kit = list()
      bpl$Header$kit = ""
      bpl$Header$plex = ""
      bpl$Header$standards = ""
      bpl$Header$blanks = ""
      bpl$Header$controls = ""
      bpl$Header$tests = ""
      
      #bp$Header$Analytes=""
      
      for( i in levels(bpl$bp$pid)){ 
        id = bpl$bp$pid == i
        
        bp = droplevels(bpl$bp[id,])
        cytokines = unique( bp$Cytokine)
        nm = paste(cytokines, collapse=',')
        #if(nchar(nm) > 30) nm = substr(nm,1,30)
        kit[[nm]] = cytokines
        bpl$Header$kit[bpl$Header$pid == i] = sprintf("k%03d",which(names(kit) == nm))
        bpl$Header$plex[bpl$Header$pid == i]  = length(cytokines)
        bpl$Header$standards[bpl$Header$pid == i]  = length( grep('S', levels(bp$sample)))
        bpl$Header$blanks[bpl$Header$pid == i]  = length( grep('B', levels(bp$sample)))
        bpl$Header$controls[bpl$Header$pid == i]  = length( grep('C', levels(bp$sample)))
        bpl$Header$tests[bpl$Header$pid == i]  = length( grep('X', levels(bp$Type)))
        
      }
      #   ab = abbreviate(bpl$Header$file,8)
      #   bpl$files = data.frame(file = bpl$Header$file, fid = ab, row.names=NULL)
      
      #   bpl$Header$file = ab 
      #   y = which(colnames(bpl$Header) == 'file')
      #   colnames(bpl$Header)[y] ="fid"
      
      bpl$bp$kit = ""
      
      for(i in 1:length(kit)) {
        kid = sprintf("k%03d",i)
        id = bpl$bp$Cytokine %in% kit[[i]]
        bpl$bp$kit[ id] = kid
      }
      
      bpl$Header$kit = as.factor(bpl$Header$kit)
      bpl$bp$kit = as.factor(bpl$bp$kit)
      
      bpl$kits = data.frame(Analytes = names(kit),
                            kit = sprintf("k%03d",1:length(kit)))
      
      str = c( 'standards', 'blanks', 'controls')
      bpl$kits = merge(bpl$kits, unique(bpl$Header[, c('kit', 'plex', str)]), by='kit')[,c('kit', 'plex', 'Analytes', str)]
      return(bpl)
      
    }
    #levels(bpl$bp$Cytokine) = nCytNme(levels(bpl$bp$Cytokine))
    t = rbind(Readings=t1<-table(bpl$bp[,"sample"]),
              Missing=t2<-table(bpl$bp[is.na(bpl$bp$FI),"sample"]), For.Analysis=t1-t2)
    #print(t)
    
    bpl$Readings = t
    
    bpl$bp = bpl$bp[!is.na(bpl$bp$FI),]
    bpl = bplxAssignCytSet(bpl)
    
    bpl$bp = bplxOrderFactors(droplevels(bpl$bp))
    
    
    return(bpl)
  }
  
  #code starts here
  bplxInstalls('openxlsx')
  rtn = data.frame(); cnt = 0
  
  header = data.frame()
  for(i in files) {
    #print(i)
    x = NULL
    wb = loadWorkbook(i)
    for(j in sheetName) {
      #print('Here')
      s = bplxGetAnalytes(wb, sheetName=j)
      
      #colnames(s) = make.names(colnames(s))
      
      if(!is.na(exclude[1])) {
        if(length(k <- which(colnames(s) %in% exclude)))
          s = s[, -k]
      }
      n = which(colnames(s) %in% c("Type", "Well", "Description"))
      s= stackIt(s, description=colnames(s)[n[length(n)]], parName = "Cytokine", parVal = j)
      if(!is.null(s)) {
            if(is.null(x)) x = s
            else  x[,j] = s[,j]
      }
    }
    # also generate a plate identifier
    pid = FUN(i, cnt <-cnt +1)
    if(!is.null(x)) {
          rtn = rbind(rtn, cbind(pid = pid, x))
          header = rbind(header, data.frame(getBplxHeader(wb), file=i,pid=pid))
    }
  }
  
  header$pid = factor(as.character(header$pid))
  rtn$pid = factor(as.character(rtn$pid))
  #colnames(rtn) = make.names(colnames(rtn))
  
  
  #   levels(bp$pid) = 1:60
  #   xod = order( suppressWarnings(as.numeric(gsub("[A-Z]", "", t[cx]))))
  rtn = bplxOrderFactors(rtn)
  if(length(which(colnames(rtn$bp) == 'Dilution')) > 0)
    rtn$Dilution = factor(rtn$Dilution, ordered=T)
  
  return(processbplxInput(list(Header = header, bp=rtn)))
}
bplxMER = function(bp, bp.mer, alpha = 0.01,cex = 1, adjust = FALSE, doOutlier = T)
{
  par(las =1)
  myfunc = function(x,ci,effect, cex) { 
    
    pv   <- attr(x, "postVar")
    cols <- 1:(dim(pv)[1])
    se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    x1 = x[,1] - 2*se
    x2 = x[,1] + 2*se
    
    ci = abs(ci)
    minv = floor(min(-ci,min(x1)))
    xlim = c(minv, ceiling(max(ci,max(x2))))
    INC = ( par("fin")[2])
    FS = (par("cin")[2])* nrow(x)
    MagicNumber = INC/FS
    
  
    dotchart(x[,1], 
             xlab="Random Effect", 
             xlim=xlim,
             cex = cex,
             col = ifelse(MagicNumber < 1, 'grey40', 'black'),
             labels=rownames(x), 
             pch=16
    )
    for(i in 1:nrow(x)) {
      lines(x=c(x1[i], x2[i]), y = c(i,i), col='blue', lwd=par("cex")*2)
    }
    mtext(effect,3, line=0.3, at = minv ,font=2, cex=cex * 1.5)
    abline(v=0, col='red',lwd=1.2)
    abline(v= ci, col='red',lwd=2,lty=4)
    abline(v= -ci, col='red',lwd=2,lty=4)
    
  }
  
  
  compareNA <- function(v1,v2) {
    same <- (v1 == v2) | (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
  }

  
  
  
  xcov = as.data.frame(VarCorr(bp.mer))
  ranf= ranef(bp.mer, condVar=T)
  nms = names(ranf)
  df = data.frame()
  Action = NULL;
  for(k in nms) {
    kp = ranef(bp.mer, condVar=T)[[k]]
    effect = k
    std = xcov[xcov[,1] == k,][1,5]
    #stderr = xcov[k,5]/sqrt(nrow(kp))
    z =  kp[,1]/ std
    
    zscore = abs(qnorm(p=alpha/2))  # two tailed test
    
    myfunc(ranf[[k]], (zscore) * std, effect,cex)
    if(doOutlier==F) {
      qqnorm(kp[,1], pch=16, col='blue'); qqline(kp[,1], col='red', lwd=2) 
    }
    
    if(doOutlier) {
      pv = 2*pnorm(-abs(z))  # because two tailed test
      
      if(adjust) 
        pv = p.adjust(pv, method='fdr')
      
      
      symp <- symnum(pv, corr = FALSE,
                     cutpoints = c(0,  0.001,.01,.05, .1, 1),
                     symbols = c("***","**","*","."," "))
      
      
      df = rbind(df, data.frame(Effect = effect, level = rownames(kp),
                                Zscore =z,
                                p.value = pv,
                                sig. = symp))
      
      idz = pv  < alpha
      Action =c(Action, ifelse( idz, 'REMOVE', 'ok') )
      
      
      if(length(which(idz))) {
        d = data.frame(idz = idz, x = kp[,1], y=1:nrow(kp), lab = rownames(kp))
        a = droplevels(d[ d$idz == T,])
        text(a$x, a$y, a$lab,col='darkblue', pos=4,font=2, offset=1, cex=cex)
      }
    }
    
  }
  
  if(doOutlier) {
    id = rep(F, nrow(bp))
    df$Action = Action
    y = droplevels(df[ compareNA(df$Action,'REMOVE'), ])
    
    if(nrow(y) > 0) {
      for(i in 1:nrow(y)) {
        
        effect = as.character(y[i,'Effect'])
        level = as.character(y[i,'level'])
        if( effect %in% colnames(bp))
          id = id | bp[, effect] == level
        else {
          
          s = strsplit(effect,':')[[1]]
          f = interaction(bp[,s], sep=':', drop=T)
          
          id = id | f == level
          
        }
      }
    }
    
    df[,3:4] = round(df[,3:4], 3)
    # par(pa)
    return(invisible(list( remIdx = id, df = df, remove = y)))
  }
  
}
bplxPlateNos = function(ns, max.spp=37, rp = 1, sample='X', code='X')
{
  # rp replication level
  nsa = ns * rp; # total number of assays needed
  
  pr = nsa/max.spp  # number of plates needed
  rem = pr%%1 
    
  np = ceiling(pr)  # actual number of plates needed
  app = floor (nsa/np)  # assays per plate
  
  nnp= ceiling(ns/np) # number of samples per plate

  plate = rep(app, np)
  samples = rep(nnp,np)

  rem = nsa - app * np
  rem2 = ns - sum(samples)
  
  if(rem > 1) {
    pp = sample(1:np, rem)
    plate[pp] = plate[pp] + 1
  }
  
  # if(rp > 1 ) {
  #   pp = sample(1:np, rem2)
  #   samples[pp] = samples[pp] + 1
  # }
  
  TYPE = lapply(samples, function(x) noquote(sprintf("%s%d",code,1:x)))

  pid = sprintf("p%03d",1:np)
 
  #names(TYPE) = pid

  fr = max.spp-plate
  plates =data.frame( pid = pid, Na=plate, Ns = samples, Rep = round(plate/samples,3), fr)
  
  lst = list(sample=data.frame(Nsample = ns, max.spp = max.spp, rl = rp, Plates = np, Assays = nsa, FreeWells = sum(fr)),
      plates = plates
       )

  for(i in 1:length(TYPE)) {
    type = TYPE[[i]]
    id = (((0:plates[i, 'Na'])  %% plates[i, 'Ns'])+1)
    id = id[-length(id)]
    
    TYPE[[i]] = type[id]
  }


  lst = c(lst, type =list(TYPE))
  
  lst
}

bplxPlateNosRep = function(ns, max.spp = 76, rp = 2, rndize = F)
{

  ip = as.integer(rp)
  fp = (rp %% 1)/ip
  
  mxp = floor(max.spp/rp)

  
  
  X = sprintf("X%d", 1:mxp)
  Xr = rep(X, each=ip)
  
  if(fp > 0) {
     X = c(Xr, X[1:floor(fp * mxp)])
    
  }
  
  
  print(X)
  print(length(X))
  
  return(list(samples=bplxPlateNos(ns, max.spp), X=X))
  
 
  
}


bplxReplicates = function(ins){
  
  str=as.character(ins)
  od = order(str)
  s = rep(0,length(str))
  s[od] = sequence(rle(str[od])$length)
  
  return (data.frame(Rep=sprintf("R%d",s),ins))
}


bplxICC = function(bp, clevel = "Analyte", str="(1|Description)", 
                   DV = 'log2(FI)', rnd=2)
{
  bplxInstalls('e1071')

  dolog = length(grep('log2', DV))>0
  
  f = as.formula(paste(DV, "~", str, sep='' ))
  
  lst = list()
  
  #fkit = factor(bp$kit:bp$Cytokine)
  
  for(kit in levels(bp$kit)) {
        kit = kit
        ibp = droplevels(bp[bp$kit==kit,])

        df = data.frame()
        for(i in levels(ibp[,clevel])) {
              id = ibp[,clevel] == i
              y = droplevels(ibp[id,])
              mod = lmer(f,y)
              
              
              t = as.data.frame(VarCorr(mod))
              n = dim(t)[1]
              betwn = sum(t[1:(n-1), 4])
              within = t[n,4]
              
              if(dolog) {
                    m = sqrt(2^betwn) * 2^(fixef(mod)[1])
                    v1 = ((2^betwn) - 1) * m^2
                    v2 = ((2^within) - 1) * m^2
                    val = log2(y$FI)
                    
              } else{
                    m = fixef(mod)[1]
                    v1 = betwn
                    v2 = within
                    val = y$FI
              }
              
             
              #icc =  t[1]/sum(t)
             
              icc = betwn/(betwn+within)
              
              # sem = sqrt(sum(t)) * sqrt(1-icc) 
              
              # q = quantile(y$FI, 0.95)
              df = rbind(df, data.frame(Cytokine=i, 
                                        N = nrow(y),
                                        mean = m, 
                                        med = median(y$FI), min = min(y$FI), 
                                        sd.s = sqrt(v1), sd.e = sqrt(v2), ICC = icc, 
                                        # SEM = sem, MD = sem * 1.96 * sqrt(2),
                                        log2.skew = skewness(log2(y$FI)),
                                        SNR = v1/v2, CV = v2/v1))
        }
        rownames(df) = NULL
        df[,-1] = round(df[,-1], rnd)
        lst[[kit]] = df[ order(df$SNR), ]
  }
  return(lst)
}


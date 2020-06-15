############### SYNOPSIS ###################
### Generate GTEx figure

### === STEPS performed in script - behind the scenes === ###
# 1) Read DEPICT tissue enrichment results.
# 2) Auto-detect format type
# 3) Process files and build plot variables
# 4) Construct plots via ggplot2
# 5) Export plots to pdf files

################## LOAD PACKAGES and ARGUMENTS #####################
#suppressPackageStartupMessages()
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid)) # needed for "unit" function

#options(warn=-1) # turn off warnings globally

######################################################
################# GLOBAL variables ###################

n.blank_spacers <- 4 # integer | number of empty bars between "groups"

### =========== INTERNAL PARAMETERS =========== ###
colnames.gtex <- make.names(c("Name","Nominal P value","False discovery rate", "Label"))

group_variable.gtex <- make.names("Name")

######################################################
##################### FUNCTIONS ######################

function.input_procces <- function(df, cols2keep) {
  ### USE: UNIVERSIAL
  ### Subsetting data ###
  df.res <- subset(df, select=cols2keep)
  ### SORTING
  #df.res <- df.res[order(df.res[,sort.colname]),] # IMPORTANT: sorting data frame
  # ^ <--- no longer needed to sort here. The functions function.*.set_variables* does this!
  return(df.res)
}

function.gtex.set_variables <- function(df) {
  ### USE: GTEx
  ### STEPS: adding new variables
  df.res <- df
  ### Set "significance" variable
  df.res$FDR.significant <- with(df.res, ifelse( (False.discovery.rate=="<0.01"|False.discovery.rate=="<0.05"), TRUE, FALSE))
  ### Set "group"/tissue variable
  tmp.strsplit <- strsplit(as.character(df.res[,group_variable.gtex]), " - ") # ALTERNATIVELY: #tmp.strsplit <- with(df.res, strsplit(as.character(Name), " - "))
  df.res$group <- sapply(tmp.strsplit, "[[", 1) # or use unlist() | or "lapply(strsplit(XX," - "), function(x) x[1])"
  
  ### ADDITIONAL VARIABLES | not needed
  #df.res <- df.res %>% group_by(group) %>% mutate(group.width=n()) # calculating the number of observations in each group ("group width")
  
  ### ***SORTING*** 
  # df.res <- df.res[order(df.res[,group_variable.gtex]),] # IMPORTANT: sorting data frame
  df.res <- df.res[order(df.res[,"group"], df.res[,"Nominal.P.value"]),] # IMPORTANT: sorting data frame | NEW 2020
  
  return(df.res)
}

function.add.spacers <- function(df) {
  ### USE: UNIVERSIAL
  ### Function to INSERT spacers (dummy rows) in data frame between "groups"
  ### *IMPORTANT*:we are relying the data frame being (correctly) SORTED!!
  df.tissue_enrichment.spaced <- data.frame()
  group.previously <- df$group[1]
  for (i in 1:nrow(df)) {
    ### *OBS*: df MUST be sorted at this point!
    group.now <- df$group[i]
    if ( group.now != group.previously ) { # add "spacer" (blank) rows between "groups"
      df.dummy <- data.frame(Name="dummy",Nominal.P.value=1,False.discovery.rate=">=0.20")
      for (j in 1:n.blank_spacers) {
        # use dplyr::bind_rows() to fill unmatched columns with NA [rbind() will complain]
        df.tissue_enrichment.spaced <- suppressWarnings(dplyr::bind_rows(df.tissue_enrichment.spaced, df.dummy))
      }
      group.previously <- df$group[i]
    }
    df.tissue_enrichment.spaced <- dplyr::bind_rows(df.tissue_enrichment.spaced, df[i,])
    # GTEx --> Name  Nominal.P.value  False.discovery.rate
  }
  return(df.tissue_enrichment.spaced)
}


function.plot.barplot <- function(df.tissue_enrichment, xlabel="MISSING X-LABEL") {
  ### USE: UNIVERSIAL
  ### Function creates a ggplot object [geom_bar()]
  ### INPUT: 
  # df.tissue_enrichment: a sorted "cleaned" data frame
  # xlabel: a character string used as xlabel
  ### OUTPUT: 
  # p: a ggplot object
  
  ### STEPs
  # 1) create breaks and labels for plot
  # 2) *call function* "function.add.spacers()"
  # 3) make plot and adjust theme
  
  
  ##########################################################
  ################## Preparing for PLOTTING ################
  
  ### Constructing labels and break positions for x-axis
  # *OBS* df.tissue_enrichment is used for input and NOT the "spaced" data frame
  df.breaks_and_labels <- df.tissue_enrichment %>% group_by(group) %>% summarize(group.width=n()) # *OBS*: this code relies on the "group" being in the correct order
  df.breaks_and_labels$order.group.numeric <- seq(0,nrow(df.breaks_and_labels)-1) # 0,1,...,3
  df.breaks_and_labels$group.width.cumsum <- with(df.breaks_and_labels, cumsum(group.width)) # cumsum()
  df.breaks_and_labels$group.width.cumsum.shift <- with(df.breaks_and_labels, c(0,group.width.cumsum[-length(group.width.cumsum)])) # *OBS*: "shifting" position by one (pop array). Inserting zero in first position | try diff(x)?
  df.breaks_and_labels$break_position <- with(df.breaks_and_labels, 0.5+(n.blank_spacers*order.group.numeric)+group.width.cumsum.shift+group.width/2) 
  
  ### Adding spacers to data frame | *OBS*: we are relying on CORRECT SORTING of the data frame
  df.tissue_enrichment.spaced <- function.add.spacers(df.tissue_enrichment)
  df.tissue_enrichment.spaced$order.numeric <- 1:nrow(df.tissue_enrichment.spaced) # Add numeric code | maps to x-axis
  #str(df.tissue_enrichment.spaced)
  ##########################################################
  ########################## PLOT ##########################
  

  # ===================== CORE PLOT ====================== # 
  p <- ggplot(df.tissue_enrichment.spaced)
  p <- p + geom_bar(aes(x=order.numeric, y=-log10(Nominal.P.value), fill=FDR.significant), stat="identity")
  p <- p + geom_hline(yintercept=-log10(0.05/54), linetype="dashed", color="darkgray") # *********HACK****************
  p <- p + geom_text(data=df.tissue_enrichment.spaced %>% filter(False.discovery.rate=="<0.01"), aes(x=order.numeric, y=-log10(Nominal.P.value), label=Label), angle=25, nudge_y=0.05, size=rel(2), hjust=0, show.legend=F)
  p <- p + geom_text(data=df.tissue_enrichment.spaced %>% filter(Name=="Brain - Hypothalamus"), aes(x=order.numeric, y=-log10(Nominal.P.value), label=Label), angle=25, nudge_y=0.05, size=rel(2), hjust=0, show.legend=F)
  p <- p + scale_fill_manual(name="", values=c("TRUE"="#de2d26","FALSE"="#606060",guide='legend')) # set colors for FDR significant bars
  p <- p + scale_x_continuous(breaks=df.breaks_and_labels$break_position, labels=df.breaks_and_labels$group) # add x-axis labels
  p <- p + labs(y=expression(-log[10](P[S-LDSC]))) # add axis titles
  p <- p + guides(fill=FALSE) # remove legend
  
  # =============== FINE TUNING - theme() =============== # 
  ### Do not out-comment individual lines in the below block of code. It should either be out-commented or un-commented all together
  p.theme <- theme_classic() # default is "theme_gray()" "theme_classic()" theme removes background and more [white background, eliminates background, gridlines, and chart border]. Note that it contains black colored axis lines. Therefore it is used as a template for further modifications
  # OR --> p <- p + theme_bw() + theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank()) 
  # See also theme_classic(), theme_minimal() and theme_bw()
  # Try the "theme() %+replace% theme()" operator
  
  ### Setting axis line for y-axis (only) | The following did not work --> theme(axis.line.y=element_line(color="black"))
  p.theme$axis.line.x <- element_blank() # empty list of class "element_blank" and "element" | SAME as theme_bw()$axis.line
  #p.theme$axis.line.x <- p.theme$axis.line # copy list
  #p.theme$axis.line.x$colour <- 'black'
  p.theme$axis.line.y <- p.theme$axis.line # copy list
  p.theme$axis.line.y$colour <- 'black' # OBS: sensitive to spelling: "colour" and NOT "color"
  
  
  ### Setting axis tick marks for y-axis (only)
  p.theme$axis.ticks.x <- element_blank()
  p.theme$axis.ticks.y <- p.theme$axis.ticks
  p.theme$axis.ticks.y$colour <- 'black' # redundant: theme_bw() has this by default
  # traditional approach --> theme(axis.ticks=element_blank())
  
  ### Adjusting length of tick marks and expansion
  p.theme$axis.ticks.length <- unit(0.15, "cm") # adjust distance from axis labels to axis | you may need to play with this a little
  # Note that setting different lengths for tickmarks for x and y does not work | that is, manually setting "axis.ticks.length.x"/"axis.ticks.length.x" does *NOT* work.
  # traditional approach --> theme(axis.ticks.length=unit(0, "cm"))
  p <- p + scale_y_continuous(expand = c(0, 0)) # removes expansion of y-axis. Now the y-axis starts at zero!
  
  p <- p + coord_cartesian(clip="off")
  
  ### Adjust x-axis labels (size and rotation)
  p.theme <- p.theme + theme(axis.text.x=element_text(angle=35, hjust=1, size=rel(0.8))) # size=rel(1.15) consider using %+replace% operator
  
  # p.theme <- p.theme
  
  # ===================== Combining theme + plot ====================== # 
  p <- p + p.theme # saving "p.theme" into plot
  
  # ===================== input_type specific ====================== # 
  p <- p + labs(x=xlabel)
  
  # ===================== return value ====================== # 
  return(p) # return ggplot object
}












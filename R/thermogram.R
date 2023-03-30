#' t.gram
#'
#' @description Reads txt thermogram files, performs data cleanup, vizualization, and helps choose fraction sizes
#'
#' @details This function can also be called from the \href{http://soilradiocarbon.org}{ISRaD website}.
#' @param txt Should be raw 'easgraph' file from soliTOC
#' @param perC Sample carbon concentration as a percent. Default is 2.
#' @param reqmgC The mass of C required for analysis of each. Used to calculate total required sample mass. Default is 0.5 mg C, for 14C analysis.
#' @param temps_in Cutoff temperatures for fractions. Function fills high and low bounds. Length should be equal to cutoff temperatures, so fraction count = n + 1.
#' @param fil The window over which the rolling mean is calculated. A higher number produces a smoother thermogram. This does not affect fraction sizes.
#' @import dplyr
#' @importFrom DescTools AUC
#' @importFrom gridExtra grid.table
#' @importFrom grid grid.text
#' @importFrom stringr str_split str_replace
#' @importFrom ("graphics", "abline", "lines", "par")
#' @importFrom ("stats", "median", "sd")
#' @importFrom ("utils", "head", "read.csv")
#' @export
#' @examples
#' \donttest{
#' # Load example dataset Gaudinski_2001
#' entry <- ISRaD::Gaudinski_2001
#' # Save as .xlsx file
#' ISRaD.save.entry(
#'   entry = entry,
#'   template_file = system.file("extdata", "ISRaD_Master_Template.xlsx", package = "ISRaD"),
#'   outfile = file.path(tempdir(), "Gaudinski_2001.xlsx")
#' )
#' # Run QAQC
#' QAQC(file.path(tempdir(), "Gaudinski_2001.xlsx"))
#' }
#'

t.gram <- function(txt, perC = 2, reqmgC = 0.5, temps_in = c(250, 300, 350, 450), fil = 20){

  library(grid)
  library(gridExtra)
  #library(mclust)
  library(zoo)
  library(stringr)
  library(DescTools)
  library(pROC)
  library(dplyr)
  #library(ggplot2)

  title = str_replace( unlist(str_split(colnames(head(read.csv(txt),1)[1]), '[...]'))[9:length(unlist(str_split(colnames(head(read.csv(txt),1)[1]), '[...]')))],
                       pattern= "\n", replacement = "")

  read.csv(txt,
         sep = '\t', skip = 6,
         col.names = c('Time', "Temp", 'Post.tube', 'Flow.N', "Flow.O","In.mBar",'Out.mBar', 'IR')) %>%
    dplyr::mutate(CO2_scaled = (IR - 1000) / max(IR - 1000)) %>%
    dplyr::mutate(CO2_scaled = ifelse(CO2_scaled < 0, 0, CO2_scaled)) %>%
    dplyr::mutate(CO2_scaled = ifelse(CO2_scaled < 0, 0, CO2_scaled)) %>%
    slice(1:grep(max(Temp), Temp)[1] + 30) %>%
    select(c("Temp", "CO2_scaled")) %>%
    group_by(Temp) %>%
    summarise(Temp_av = mean(CO2_scaled),
              Temp_med = median(CO2_scaled),
              Temp_sd = sd(CO2_scaled)) %>%
    ungroup() %>%
    right_join(data.frame(Temp = seq(50, 1000))) %>%
    arrange(Temp) %>%
    mutate(Temp_av = ifelse(Temp > grep(max(Temp), Temp)[1], 0, Temp_av)) %>%
    mutate(Temp_med = ifelse(Temp > grep(max(Temp), Temp)[1], 0, Temp_med)) %>%
    mutate(Temp_sd = ifelse(Temp > grep(max(Temp), Temp)[1], 0, Temp_sd)) %>%
    mutate(sp.Mean = na.spline(Temp_av, xout = seq(50, 1000))) %>%
    mutate(sp.Med = na.spline(Temp_med, xout = seq(50, 1000))) %>%
    dplyr::mutate(sp.Mean = ifelse(sp.Mean < 0, 0, sp.Mean)) %>%
    dplyr::mutate(sm.Mean = rollmean(sp.Mean, k = fil, fill = 0)) %>%
    dplyr::mutate(sm.Med = rollmean(sp.Med, k = fil, fill = 0)) -> t.gram

  area.tot = AUC(t.gram$Temp, t.gram$sp.Mean)
  area.vec = c(0)
  for(i in seq(min(t.gram$Temp) + 1, max(t.gram$Temp))){
    area.vec = c(area.vec, AUC(t.gram$Temp, t.gram$Temp_av, from = t.gram$Temp[i], to = t.gram$Temp[i + 1]) / area.tot)
  }

  t.gram$Area <- area.vec
  t.gram$AreaCumu <- cumsum(t.gram$Area)

  #### Calculate fraction sizes

  fracs <- c(0, dplyr::filter(t.gram, Temp %in% temps_in)$AreaCumu, 1)

  par(mfrow = c(2,1))
  plot(t.gram$Temp, t.gram$sm.Mean, type = 'l', col = 'white', lwd = 4,
       xlab = 'Temperature (C)', ylab = 'Normalized C', main = title, xlim = c(50, 800))
  abline(,,0)
  lines(t.gram$Temp, t.gram$sm.Mean, type = 'l', col = 'steelblue4', lwd = 4)
  abline(,,,dplyr::filter(t.gram, Temp %in% temps_in)$Temp, lty = 2, lwd = 3, col = 'brown4')

  pvp = viewport(x = .3, y = .3)
  pushViewport(pvp)
  grid.table(as.matrix(data.frame(`Upper Temp` = c(dplyr::filter(t.gram, Temp %in% temps_in)$Temp, 1000),
                                  Prop = round(diff(fracs),4))))
  grid.text(paste("In order to collect", reqmgC, "mg C\n in the smallest fraction,\n",
                  toString(round(reqmgC * (1/ min(fracs[-1])) * (100/perC))), "mg sample required."), x=.9, y = .55)

  return(t.gram)

}

tout = t.gram('/Users/shane/14Constraint Dropbox/Shane Stoner/IMPRS/ThermalAnalysis/Ch3/Thermograms/Feb2022/easgraph007.txt',
              fil = 15)

tout %>%
  select(Temp, sp.Mean, sm.Mean) %>%
  write.csv('/Users/shane/Desktop/ANwf_thermogram.csv')

library(ggplot2)
library(data.table)
library(plotly)
library(htmlwidgets)

result_path <- "C:/Santos/VIH/Prueba_push/data7/output/20250311/EFA_results"
setwd(result_path)


for(i in 1:5) {
  if(!dir.exists(paste0("Prev_", i))){dir.create(paste0("Prev_", i))}
  files <- list.files(pattern = paste0(i,"%.*\\.csv$"))  
  
  for(resultFile in files){
    t <- fread(resultFile)
    setorder(t, V1)
    t[, V1 := factor(x = V1, levels = rev(V1))] # Invert factor level for alphabetical 
    temp <- melt(data = t, id.vars = "V1")
    temp[value == 0, value := NA]
    temp[, text := paste0("y: ", V1, "\n", "x: ", variable, "\n", "Value: ", value)]
    
    p <- ggplot(temp, aes(x = variable, y = V1, fill = value, text = text)) + 
      geom_tile() + theme_bw() + theme(
        axis.title = element_blank(),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        panel.background = element_blank(), # Remove background color
        axis.line = element_line(color = "black") # Add axis lines
      )
    
    pp <- ggplotly(p, tooltip = "text")
    fileName <- gsub("csv$", "html", resultFile)
    fileName <- gsub("efa_loadings_", "", fileName)
    fileName <- gsub("%", "prev", fileName)
    
    saveWidget(pp, file=paste0(getwd(),"/Prev_", i, "/", fileName)) 
    # full path without special characters, % not allowed.

  }
}


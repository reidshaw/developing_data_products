library(data.table); library(reshape2); library(ggplot2); library(plotly); library(shiny)


shinyServer(function(input, output) {
     
     
     output$interactive_plot <- renderPlotly({
     
          x <- list(title = "Log2 (KRAS-mut / KRAS-wt)", zeroline = TRUE,
                    showline = TRUE,
                    mirror = "ticks",
                    gridcolor = toRGB("white"),
                    gridwidth = 0,
                    zerolinecolor = toRGB("black"),
                    zerolinewidth = 1,
                    linecolor = toRGB("black"),
                    linewidth = 3)
          y <- list(title = "-Log10 p-value", zeroline = TRUE,
                    showline = TRUE,
                    mirror = "ticks",
                    gridcolor = toRGB("white"),
                    gridwidth = 0,
                    zerolinecolor = toRGB("black"),
                    zerolinewidth = 1,
                    linecolor = toRGB("black"),
                    linewidth = 3)
          
          p <- plot_ly(df, x = ~log_2_fc, y = ~neg_log_pval, 
                       type = "scatter",
                       text = ~paste('Compound: ', drug_common, 
                                     '</br> Target: ', drug_target),
                       color = ~drug_pathway,
                       size = ~abs_diff * 2) %>%
               layout(title = "Lund Adenocarcinoma",
                      hovermode="closest", xaxis = x, yaxis = y) 
          p
          
     })
     
     output$boxplot <- renderPlot({
          b <- response[response$DRUG.NAME == input$drug,]
          c <- df[df$drug_common == input$drug,]
          
          a <- boxplot(b$AUC ~ b$KRAS_status, las = 1,
                       names=c("KRAS-WT", "KRAS-mut"), outcol = "white",
                       main = input$drug,
                       ylab = "AUC", 
                       ylim = c(0,1.2))
          stripchart(AUC ~ KRAS_status, data = b, vertical = TRUE,
                     method = "jitter", add = TRUE, pch = 19, col = alpha('steelblue2', 1/2))
          text(x = 1.5, y = 0.1, labels = paste("p-value: ", round(c$p_value,3)))
          text(c(1:2) , a$stats[nrow(a$stats), ] +0.1 , paste("n = ", table(b$KRAS_status), sep=""))
     })
     
})

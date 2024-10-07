library(shiny)
library(ggplot2)
library(dplyr)
library(patchwork)

calculate_model_stats <- function(maf, model, N) {
  q <- maf
  p <- 1 - q
  freq_p2 <- p^2
  freq_2pq <- 2 * p * q
  freq_q2 <- q^2
  AA <- round(N * freq_p2)
  AB <- round(N * freq_2pq)
  BB <- round(N * freq_q2)
  
  add <- c(rep(0, AA), rep(1, AB), rep(2, BB))
  y <- ifelse(add == 2, model[3], ifelse(add == 1, model[2], model[1]))
  
  a <- sum(add == 2) / length(add)
  h <- sum(add == 1) / length(add)
  r <- sum(add == 0) / length(add)
  
  dom_hom_ref <- -h * a
  dom_het <- 2 * a * r
  dom_hom_alt <- -h * r
  
  dom <- add
  dom[which(add == 0)] <- dom_hom_ref 
  dom[which(add == 1)] <- dom_het
  dom[which(add == 2)] <- dom_hom_alt
  
  rec <- ifelse(add == 2, 1, 0)
  dom_011 <- ifelse(add > 0, 1, 0)
  
  model_add <- lm(y ~ add)
  model_dom <- lm(y ~ dom)
  
  summary_add <- summary(model_add)
  summary_dom <- summary(model_dom)
  
  return(list(
    genotypes = add,
    phenotypes = y,
    model_add = model_add,
    model_dom = model_dom,
    frequencies = c(AA, AB, BB) / N,
    counts = c(AA, AB, BB),
    r_squared_add = summary_add$r.squared,
    r_squared_dom = summary_dom$r.squared,
    p = p,
    q = q
  ))
}

plot_genetic_architecture <- function(maf, model, N) {
  stats <- calculate_model_stats(maf, model, N)
  
  df <- data.frame(
    Genotype = c(0, 1, 2),
    Real = model,
    Frequency = stats$frequencies,
    Count = stats$counts
  )
  
  genotype_seq <- seq(-0.1, 2.1, length.out = 100)
  additive_values <- predict(stats$model_add, newdata = data.frame(add = genotype_seq))
  additive_df <- data.frame(Genotype = genotype_seq, Additive = additive_values)
  
  # Create the updated title
  title <- sprintf("p=%.2f, q=%.2f | Model: [%.2f, %.2f, %.2f]", 
                   stats$p, stats$q, model[1], model[2], model[3])
  
  main_plot <- ggplot() +
    geom_point(data = df, 
               aes(x = Genotype, y = Real, size = Frequency),
               color = "black", fill = "#7C9DCD", shape = 21, stroke = 1) +
    geom_line(data = df, 
              aes(x = Genotype, y = Real, color = "Real", linetype = "Real"),
              size = 0.6) +
    geom_line(data = additive_df, 
              aes(x = Genotype, y = Additive, color = "Additive", linetype = "Additive"),
              size = 1.2) +
    scale_color_manual(values = c("Real" = "black", "Additive" = "#30589E")) +
    scale_linetype_manual(values = c("Real" = "solid", "Additive" = "dashed")) +
    scale_size_continuous(range = c(2, 20), name = "Genotype\nFrequency") +
    scale_x_continuous(breaks = c(0, 1, 2), limits = c(-0.2, 2.2)) +
    scale_y_continuous(limits = c(min(model) - 0.1, max(model) + 0.2)) +
    labs(
      title = title,
      x = "Genotype",
      y = "Phenotype"
    ) +
    theme_classic() +
    theme(
      legend.position = "none"
    )
  
  variance_df <- data.frame(
    Model = c("Additive", "Non-additive"),
    Variance = c(stats$r_squared_add, stats$r_squared_dom)
  )
  
  variance_plot <- ggplot(variance_df, aes(x = 1, y = Variance, fill = Model)) +
    geom_bar(stat = "identity", width = 0.5) +
    coord_flip() +
    scale_fill_manual(values = c("Additive" = "#30589E", "Non-additive" = "#BA2025")) +
    theme_void() +
    theme(legend.position = "none") +
    geom_text(aes(label = sprintf("%.1f%%", Variance * 100)), 
              position = position_stack(vjust = 0.5), color = "white") +
    labs(title = "Variance Explained")
  
  combined_plot <- variance_plot / main_plot + plot_layout(heights = c(1, 10))
  
  return(list(
    plot = combined_plot,
    r_squared_add = stats$r_squared_add,
    r_squared_dom = stats$r_squared_dom,
    genotype_counts = stats$counts,
    genotype_frequencies = stats$frequencies
  ))
}

ui <- fluidPage(
  titlePanel("Interactive Genetic Architecture Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("maf", "Minor Allele Frequency:", 
                  min = 0.01, max = 0.5, value = 0.3, step = 0.01),
      sliderInput("wildtype", "Wildtype Effect:", 
                  min = -2, max = 2, value = 0, step = 0.01),
      sliderInput("heterozygous", "Heterozygous Effect:", 
                  min = -2, max = 2, value = 0.01, step = 0.01),
      sliderInput("homozygous", "Homozygous Effect:", 
                  min = -2, max = 2, value = 1, step = 0.01),
      numericInput("N", "Population Size:", 
                   min = 1000, max = 1000000, value = 10000)
    ),
    
    mainPanel(
      plotOutput("geneticPlot")
    )
  )
)

server <- function(input, output) {
  result <- reactive({
    model <- c(input$wildtype, input$heterozygous, input$homozygous)
    plot_genetic_architecture(input$maf, model, input$N)
  })
  
  output$geneticPlot <- renderPlot({
    result()$plot
  })
}

shinyApp(ui = ui, server = server)
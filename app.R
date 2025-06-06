library(shiny)
library(ggplot2)
library(dplyr)
library(patchwork)
library(bslib)

# Function to calculate HWE chi-square test
calculate_hwe_test <- function(n_AA, n_AB, n_BB) {
  n_total <- n_AA + n_AB + n_BB
  if (n_total == 0) return(list(p_value = 1, chi_sq = 0))
  
  # Observed allele frequencies
  p <- (2 * n_AA + n_AB) / (2 * n_total)
  q <- 1 - p
  
  # Expected counts under HWE
  exp_AA <- n_total * p^2
  exp_AB <- n_total * 2 * p * q
  exp_BB <- n_total * q^2
  
  # Chi-square test (1 degree of freedom for HWE)
  observed <- c(n_AA, n_AB, n_BB)
  expected <- c(exp_AA, exp_AB, exp_BB)
  
  # Avoid division by zero
  chi_sq <- sum((observed - expected)^2 / pmax(expected, 1e-10))
  p_value <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
  
  return(list(p_value = p_value, chi_sq = chi_sq, allele_freq_p = p, allele_freq_q = q))
}

calculate_model_stats_flex <- function(n_AA, n_AB, n_BB, model) {
  n_total <- n_AA + n_AB + n_BB
  
  # Create genotype and phenotype vectors
  add <- c(rep(0, n_AA), rep(1, n_AB), rep(2, n_BB))
  y <- ifelse(add == 2, model[3], ifelse(add == 1, model[2], model[1]))
  
  # Calculate proportions for dominance encoding
  a <- n_BB / n_total  # homozygous alternate
  h <- n_AB / n_total  # heterozygous
  r <- n_AA / n_total  # homozygous reference
  
  # Dominance encoding (from your paper's equation)
  dom_hom_ref <- -h * a
  dom_het <- 2 * a * r
  dom_hom_alt <- -h * r
  
  dom <- add
  dom[which(add == 0)] <- dom_hom_ref 
  dom[which(add == 1)] <- dom_het
  dom[which(add == 2)] <- dom_hom_alt
  
  # Fit models
  model_add <- lm(y ~ add)
  model_dom <- lm(y ~ dom)
  
  summary_add <- summary(model_add)
  summary_dom <- summary(model_dom)
  
  # HWE test
  hwe_result <- calculate_hwe_test(n_AA, n_AB, n_BB)
  
  return(list(
    genotypes = add,
    phenotypes = y,
    model_add = model_add,
    model_dom = model_dom,
    frequencies = c(r, h, a),
    counts = c(n_AA, n_AB, n_BB),
    r_squared_add = summary_add$r.squared,
    r_squared_dom = summary_dom$r.squared,
    hwe_p_value = hwe_result$p_value,
    hwe_chi_sq = hwe_result$chi_sq,
    allele_freq_p = hwe_result$allele_freq_p,
    allele_freq_q = hwe_result$allele_freq_q,
    total_n = n_total
  ))
}

plot_genetic_architecture_flex <- function(n_AA, n_AB, n_BB, model) {
  stats <- calculate_model_stats_flex(n_AA, n_AB, n_BB, model)
  
  df <- data.frame(
    Genotype = c(0, 1, 2),
    Real = model,
    Frequency = stats$frequencies,
    Count = stats$counts
  )
  
  genotype_seq <- seq(-0.1, 2.1, length.out = 100)
  additive_values <- predict(stats$model_add, newdata = data.frame(add = genotype_seq))
  additive_df <- data.frame(Genotype = genotype_seq, Additive = additive_values)
  
  # Create title with HWE information
  hwe_status <- ifelse(stats$hwe_p_value < 0.05, 
                       sprintf("HWE: p=%.2e (VIOLATED)", stats$hwe_p_value),
                       sprintf("HWE: p=%.3f", stats$hwe_p_value))
  
  title <- sprintf("Counts: [%d, %d, %d] | %s | Model: [%.2f, %.2f, %.2f]", 
                   n_AA, n_AB, n_BB, hwe_status, model[1], model[2], model[3])
  
  # Subtitle with allele frequencies
  subtitle <- sprintf("Allele frequencies: p=%.3f, q=%.3f", 
                      stats$allele_freq_p, stats$allele_freq_q)
  
  main_plot <- ggplot() +
    geom_point(data = df, 
               aes(x = Genotype, y = Real, size = Frequency),
               color = "black", fill = "#3498DB", shape = 21, stroke = 1.5) +
    geom_line(data = df, 
              aes(x = Genotype, y = Real, color = "Real", linetype = "Real"),
              size = 1.2, alpha = 0.8) +
    geom_line(data = additive_df, 
              aes(x = Genotype, y = Additive, color = "Additive", linetype = "Additive"),
              size = 1.5) +
    scale_color_manual(values = c("Real" = "#2C3E50", "Additive" = "#E74C3C")) +
    scale_linetype_manual(values = c("Real" = "solid", "Additive" = "dashed")) +
    scale_size_continuous(range = c(3, 25), name = "Genotype\nFrequency") +
    scale_x_continuous(breaks = c(0, 1, 2), limits = c(-0.2, 2.2)) +
    scale_y_continuous(limits = c(min(model) - 0.15, max(model) + 0.25)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Genotype (Number of Variant Alleles)",
      y = "Phenotypic Effect"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, color = "#2C3E50", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "#7F8C8D"),
      axis.title = element_text(size = 11, color = "#34495E", face = "bold"),
      axis.text = element_text(size = 10, color = "#2C3E50"),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.5),
      panel.grid.minor = element_line(color = "#F8F9FA", size = 0.3),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  variance_df <- data.frame(
    Model = c("Additive", "Non-additive"),
    Variance = c(stats$r_squared_add, stats$r_squared_dom)
  )
  
  variance_plot <- ggplot(variance_df, aes(x = 1, y = Variance, fill = Model)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", size = 1) +
    coord_flip() +
    scale_fill_manual(values = c("Additive" = "#3498DB", "Non-additive" = "#E74C3C")) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, color = "#2C3E50"),
      plot.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(t = 10, b = 0)
    ) +
    geom_text(aes(label = sprintf("%.1f%%", Variance * 100)), 
              position = position_stack(vjust = 0.5), 
              color = "white", size = 4, fontface = "bold") +
    labs(title = "Variance Explained") +
    theme(plot.title = element_text(hjust = 0.5, color = "#2C3E50", face = "bold", size = 11)) +
    guides(fill = guide_legend(direction = "horizontal", 
                               override.aes = list(size = 0.8)))
  
  combined_plot <- variance_plot / main_plot + plot_layout(heights = c(1, 10))
  
  return(list(
    plot = combined_plot,
    r_squared_add = stats$r_squared_add,
    r_squared_dom = stats$r_squared_dom,
    genotype_counts = stats$counts,
    genotype_frequencies = stats$frequencies,
    hwe_p_value = stats$hwe_p_value,
    hwe_chi_sq = stats$hwe_chi_sq
  ))
}

ui <- fluidPage(
  theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
  
  titlePanel(
    div(
      h1("Interactive Genetic Architecture Explorer", 
         style = "color: #2C3E50; font-weight: bold; margin-bottom: 10px;"),
      h4("Visualizing Additive vs. Non-additive Genetic Effects", 
         style = "color: #7F8C8D; font-weight: normal; margin-top: 0px;")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      style = "background-color: #F8F9FA; padding: 20px; border-radius: 10px;",
      
      div(
        style = "background-color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
        h4("🧬 Input Method", style = "color: #2C3E50; margin-top: 0px;"),
        radioButtons("input_method", NULL,
                     choices = list("Assume Hardy-Weinberg Equilibrium" = "hwe",
                                    "Custom Genotype Counts" = "counts"),
                     selected = "hwe")
      ),
      
      conditionalPanel(
        condition = "input.input_method == 'hwe'",
        div(
          style = "background-color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
          h4("📊 HWE Parameters", style = "color: #2C3E50; margin-top: 0px;"),
          sliderInput("maf", "Minor Allele Frequency:", 
                      min = 0.01, max = 0.5, value = 0.3, step = 0.01,
                      ticks = TRUE),
          sliderInput("N_hwe", "Population Size:", 
                      min = 1000, max = 50000, value = 10000, step = 500,
                      ticks = TRUE)
        )
      ),
      
      conditionalPanel(
        condition = "input.input_method == 'counts'",
        div(
          style = "background-color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
          h4("🔢 Custom Genotype Counts", style = "color: #2C3E50; margin-top: 0px;"),
          div(
            style = "margin-bottom: 15px;",
            div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                span("🟦", style = "margin-right: 8px; font-size: 16px;"),
                strong("Wildtype (0/0):", style = "color: #34495E;")
            ),
            sliderInput("n_AA", NULL,
                        min = 0, max = 200000, value = 200000, step = 2,
                        ticks = TRUE)
          ),
          div(
            style = "margin-bottom: 15px;",
            div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                span("🟨", style = "margin-right: 8px; font-size: 16px;"),
                strong("Heterozygous (0/1):", style = "color: #34495E;")
            ),
            sliderInput("n_AB", NULL,
                        min = 0, max = 10000, value = 100, step = 2,
                        ticks = TRUE)
          ),
          div(
            style = "margin-bottom: 15px;",
            div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                span("🟥", style = "margin-right: 8px; font-size: 16px;"),
                strong("Homozygous (1/1):", style = "color: #34495E;")
            ),
            sliderInput("n_BB", NULL,
                        min = 0, max = 100, value = 5, step = 1,
                        ticks = TRUE)
          ),
          div(
            style = "background-color: #ECF0F1; padding: 10px; border-radius: 5px; text-align: center;",
            textOutput("total_n_display", container = span),
            style = "font-weight: bold; color: #2C3E50;"
          )
        )
      ),
      
      div(
        style = "background-color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
        h4("🎯 Phenotypic Effects", style = "color: #2C3E50; margin-top: 0px;"),
        sliderInput("wildtype", "Wildtype Effect:", 
                    min = -2, max = 2, value = 0, step = 0.01,
                    ticks = TRUE),
        sliderInput("heterozygous", "Heterozygous Effect:", 
                    min = -2, max = 2, value = 0.01, step = 0.01,
                    ticks = TRUE),
        sliderInput("homozygous", "Homozygous Effect:", 
                    min = -2, max = 2, value = 1, step = 0.01,
                    ticks = TRUE)
      ),
      
      div(
        style = "background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
        h4("📈 HWE Test Results", style = "color: #2C3E50; margin-top: 0px;"),
        div(
          style = "background-color: #F8F9FA; padding: 12px; border-radius: 6px; font-family: 'Courier New', monospace;",
          verbatimTextOutput("hwe_results")
        )
      )
    ),
    
    mainPanel(
      width = 8,
      div(
        style = "background-color: white; padding: 20px; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
        plotOutput("geneticPlot", height = "650px")
      )
    )
  )
)

server <- function(input, output) {
  
  # Calculate genotype counts based on input method
  genotype_counts <- reactive({
    if (input$input_method == "hwe") {
      # HWE-based calculation
      q <- input$maf
      p <- 1 - q
      n_total <- input$N_hwe
      
      n_AA <- round(n_total * p^2)
      n_AB <- round(n_total * 2 * p * q)
      n_BB <- round(n_total * q^2)
      
      # Adjust for rounding errors
      actual_total <- n_AA + n_AB + n_BB
      if (actual_total != n_total) {
        # Adjust the largest group
        largest_idx <- which.max(c(n_AA, n_AB, n_BB))
        if (largest_idx == 1) n_AA <- n_AA + (n_total - actual_total)
        else if (largest_idx == 2) n_AB <- n_AB + (n_total - actual_total)
        else n_BB <- n_BB + (n_total - actual_total)
      }
      
      return(list(n_AA = n_AA, n_AB = n_AB, n_BB = n_BB))
    } else {
      # Direct count input
      return(list(n_AA = input$n_AA, n_AB = input$n_AB, n_BB = input$n_BB))
    }
  })
  
  result <- reactive({
    counts <- genotype_counts()
    model <- c(input$wildtype, input$heterozygous, input$homozygous)
    plot_genetic_architecture_flex(counts$n_AA, counts$n_AB, counts$n_BB, model)
  })
  
  output$geneticPlot <- renderPlot({
    result()$plot
  })
  
  output$total_n_display <- renderText({
    if (input$input_method == "counts") {
      total <- input$n_AA + input$n_AB + input$n_BB
      paste("Total N =", total)
    }
  })
  
  output$hwe_results <- renderText({
    res <- result()
    paste0(
      "HWE χ² = ", round(res$hwe_chi_sq, 3), "\n",
      "HWE p-value = ", format(res$hwe_p_value, scientific = TRUE, digits = 3), "\n",
      ifelse(res$hwe_p_value < 0.05, "HWE VIOLATED", "HWE not violated")
    )
  })
}

shinyApp(ui = ui, server = server)

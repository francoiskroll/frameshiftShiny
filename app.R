library(shiny)
library(ggplot2)
library(dplyr)


# function pKos to calculate knockout probabilities -----------------------

pKos <- function (pmuts, pfras, ploidy) {
  
  # calculate each Pmiss, i.e. probability that one gRNA does not make a frameshift
  pmiss <- sapply(1:length(pmuts), function(g) { # g is 1 for gRNA 1, 2 for gRNA 2 etc.
    (1 - (pmuts[g] * pfras[g])) 
  })
  
  # then calculate Pko if gRNA 1 injected, gRNA 1 and 2, etc.
  pkos <- sapply(1:length(pmuts), function(ng) { # ng for number of gRNAs injected together
    (1 - prod(pmiss[1:ng])) ^ ploidy
  })
  
  # build a dataframe
  kodf <- data.frame(ngs=1:length(pmuts), pko=pkos)
  kodf$ngs <- as.factor(kodf$ngs)
  
  # this is what we return
  return(kodf)
  
}




# user interface ----------------------------------------------------------

ui <- fluidPage(titlePanel ('Knockout by frameshift model'),
                h5('Francois Kroll'),
                #a(href='https://twitter.com/francois_kroll', img(src='twit.png')),
                
                tags$head(tags$link(rel="icon", type="image/png", href="favicon.png")),
                
                # side bar with fields for user to enter values
                sidebarLayout(sidebarPanel(
                  h4('gRNA 1'),
                  numericInput('mut1', h6('mutation probability'), value=0.7, min=0.0, max=1.0, step=0.1),
                  numericInput('fra1', h6('frameshift probability'), value=0.66, min=0.0, max=1.0, step=0.1),
                  
                  h4('gRNA 2'),
                  numericInput('mut2', h6('mutation probability'), value=0.9, min=0.0, max=1.0, step=0.1),
                  numericInput('fra2', h6('frameshift probability'), value=0.66, min=0.0, max=1.0, step=0.1),
                  
                  h4('gRNA 3'),
                  numericInput('mut3', h6('mutation probability'), value=0.8, min=0.0, max=1.0, step=0.1),
                  numericInput('fra3', h6('frameshift probability'), value=0.66, min=0.0, max=1.0, step=0.1),
                ),
                
                # plot
                mainPanel(plotOutput('ggModel'),
                          h3('Notes'),
                          withMathJax("For gRNA 1 alone: $$P_{KO}=(1-(1-P_{mutation1} \\cdot P_{frameshift1}))^{ploidy}$$"),
                          withMathJax("For gRNA 1+2: $$P_{KO}=(1 - (1-P_{mutation1} \\cdot P_{frameshift1}) \\cdot (1-P_{mutation2} \\cdot P_{frameshift2}) )^{ploidy}$$"),
                          withMathJax("etc."),
                          p(''),
                          p('• Probability of knockout is likely an underestimation as the model assumes that frameshift is the only knockout mechanism.
                            In practice, other mechanisms may cause loss of function of the targeted gene,
                            such as mutations of amino acids at key domains or large deletions that span two target loci.'),
                          p('• Please give frameshift probability as the proportion of indels whose length is not a multiple of 3,
                            not as a proportion of all reads/genomes.'),
                          p('• If you do not know the proportions of frameshift mutations,
                            please leave frameshift probabilities as random, i.e. 2/3 = 0.66.'),
                          p('• Model assumes ploidy = 2.'),
                          h3('References'),
                          p('• I originally read about the frameshift model in', a(href="https://www.cell.com/developmental-cell/fulltext/S1534-5807(18)30457-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS153458071830457X%3Fshowall%3Dtrue", HTML("Wu et al., 2018. <em>Developmental Cell</em>"))),
                          p('• Can read more about well it matches in vivo results in', a(href="https://elifesciences.org/articles/59683", HTML('Kroll et al., 2021. <em>eLife</em>'))),
                          p('• Source code is available at', a(href="https://github.com/francoiskroll/frameshiftShiny", HTML('GitHub/frameshiftShiny'))))))



# server ------------------------------------------------------------------

server <- function (input, output, session) {
  
  # check the inputs one by one
  # req prevents app from crashing if field is emptied by the user
  # then we check that input is not below 0 or above 1
  # in which case we replace by 0 or 1, respectively
  
  # check mut1
  observeEvent(input$mut1, {
    
    req(input$mut1)
    
    if (input$mut1 < 0) {
      updateNumericInput(session, 'mut1', value=0)
    } else if (input$mut1 > 1) {
      updateNumericInput(session, 'mut1', value=1)
    }
  })
  
  # check mut2
  observeEvent(input$mut2, {
    
    req(input$mut2)
    
    if (input$mut2 < 0) {
      updateNumericInput(session, 'mut2', value=0)
    } else if (input$mut2 > 1) {
      updateNumericInput(session, 'mut2', value=1)
    }
  })
  
  # check mut3
  observeEvent(input$mut3, {
    
    req(input$mut3)
    
    if (input$mut3 < 0) {
      updateNumericInput(session, 'mut3', value=0)
    } else if (input$mut3 > 1) {
      updateNumericInput(session, 'mut3', value=1)
    }
  })
  
  # check fra1
  observeEvent(input$fra1, {
    
    req(input$fra1)
    
    if (input$fra1 < 0) {
      updateNumericInput(session, 'fra1', value=0)
    } else if (input$fra1 > 1) {
      updateNumericInput(session, 'fra1', value=1)
    }
  })
  
  # check fra2
  observeEvent(input$fra2, {
    
    req(input$fra2)
    
    if (input$fra2 < 0) {
      updateNumericInput(session, 'fra2', value=0)
    } else if (input$fra2 > 1) {
      updateNumericInput(session, 'fra2', value=1)
    }
  })
  
  # check fra3
  observeEvent(input$fra3, {
    
    req(input$fra3)
    
    if (input$fra3 < 0) {
      updateNumericInput(session, 'fra3', value=0)
    } else if (input$fra3 > 1) {
      updateNumericInput(session, 'fra3', value=1)
    }
  })

  output$ggModel <-
    renderPlot({
      
      pKos(pmuts=c(input$mut1, input$mut2, input$mut3),
           pfras=c(input$fra1, input$fra2, input$fra3),
           ploidy=2) %>%
        
      ggplot(., aes(x=ngs, y=pko, group=1, label=round(pko, 2))) +
        geom_line(size=1) +
        geom_point(size=3.5, shape=21, fill='white', stroke=1.5) +
        geom_text(hjust=1, vjust=-1, size=5) +
        coord_cartesian(ylim=c(0,1.1)) +
        ylab('probability of knockout') +
        scale_x_discrete(breaks=c(1, 2, 3), labels=c('gRNA 1 alone', 'gRNA 1+2', 'gRNA 1+2+3')) +
        scale_y_continuous(breaks=seq(0, 1.0, 0.2)) +
        theme_minimal() +
        theme(panel.grid.minor=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=18, margin=margin(t=0, r=20, b=0, l=0)),
              axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12))
    })
}


# run the app -------------------------------------------------------------

shinyApp(ui=ui, server=server)
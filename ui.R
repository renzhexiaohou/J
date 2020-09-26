#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)


# Define UI for application that draws a histogram
shinyUI(
    fluidPage(
        theme = shinytheme("flatly"),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")),   
        navbarPage(
            
            title = strong("JAB3312 real-time PKPD MS"),
            
            tabPanel(
                strong("PK"),
                fluidRow(
                    column(4, 
                           # h3("Dose regimen"),
                           fluidRow(
                               column(6),
                               column(6,
                                      sliderInput("cycle_pk", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("amt_pk1", label = "dose (mg)", min = 0, max = 30, step = 0.5, value = c(4))
                               ),
                               column(6,
                                      sliderInput("amt_pk2", label = "dose interruption (mg)", min = 0, max = 30, step = 0.5, value = c(0))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("interval_pk1", label = "dose interval (h)", min = 0, max = 168, step = 12, value = c(24))
                               ),
                               column(6,
                                      sliderInput("interval_pk2", label = "interrupted interval (h)", min = 0, max = 168, step = 12, value = c(24))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("n_pk1", label = "dose times", min = 1, max = 28, step = 1, value = c(28))
                               ),
                               column(6,
                                      sliderInput("n_pk2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                               )
                           )
                    ),
                column(8,
                       plotOutput("distPlotPK"))
                ),
                
                fluidRow(
                    column(2,
                           sliderInput("cl_pk", label = "CL (L/h)", min = 1, max = 10, step = 0.01, value = c(3.11))),
                    column(2,
                           sliderInput("v2_pk", label = "Vc (L)", min = 1, max = 200, step = 0.1, value = c(23.6))),
                    column(2,
                           sliderInput("v3_pk", label = "Vp (L)", min = 1, max = 500, step = 1, value = c(249))),
                    column(2,
                           sliderInput("q_pk", label = "Q (L/h)", min = 1, max = 100, step = 0.1, value = c(39.7))),
                    column(2,
                           sliderInput("ka_pk", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.44)))
                ),
                absolutePanel(
                    top = 63, right = 230, width = 170, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>C<sub>", "max_last", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkcmax1", inline = T), "</code> ng/mL</strong>")
                    )
                ),
                absolutePanel(
                    top = 63, right = 60, width = 170, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>C<sub>", "min_last", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkcmin1", inline = T), "</code> ng/mL</strong>")
                    )
                ),
                absolutePanel(
                    top = 63, right = 410, width = 220, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>AUC<sub>", "tau_last", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkauc1", inline = T), "</code> hr*ng/mL</strong>")
                    )
                ),
                
                fluidRow(
                        column(2,
                               numericInput("t_start", label = "NCA_lower_t1 (h)", value = c(0))),
                        column(2,
                               numericInput("t_end", label = "NCA_upper_t2 (h)", value = c(24)))
                ),
                
                absolutePanel(
                    top = 630, right = 270, width = 170, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>C<sub>", "max_t1t2", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkcmax2", inline = T), "</code> ng/mL</strong>")
                    )
                ),
                absolutePanel(
                    top = 630, right = 90, width = 170, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>C<sub>", "min_t1t2", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkcmin2", inline = T), "</code> ng/mL</strong>")
                    )
                ),
                
                absolutePanel(
                    top = 630, right = 450, width = 220, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>AUC<sub>", "tau_t1t2", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkauc2", inline = T), "</code> hr*ng/mL</strong>")
                    )
                ),
                
                absolutePanel(
                    top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                    img(src="LOGOdMed.png", height = 50)
                )
            ),
            
            
            tabPanel(
                strong("TGI"),
                fluidRow(
                    fluidRow(
                        column(4, 
                               # h3("Dose regimen"),
                               fluidRow(
                                   column(6),
                                   column(6,
                                          sliderInput("cycle_tumor", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("amt_tumor1", label = "dose (mg)", min = 0, max = 30, step = 0.5, value = c(4))
                                   ),
                                   column(6,
                                          sliderInput("amt_tumor2", label = "dose interruption (mg)", min = 0, max = 30, step = 0.5, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("interval_tumor1", label = "dose interval (h)", min = 0, max = 168, step = 1, value = c(24))
                                   ),
                                   column(6,
                                          sliderInput("interval_tumor2", label = "interrupted interval (h)", min = 0, max = 168, step = 1, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("n_tumor1", label = "dose times", min = 1, max = 28, step = 1, value = c(26))
                                   ),
                                   column(6,
                                          sliderInput("n_tumor2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                                   )
                               )
                        ),
                        column(4, plotOutput("distPlotTumorPK")),
                        column(4, plotOutput("distPlotTumorPD"))
                    ),    
                    
                    fluidRow(
                        column(2,
                               sliderInput("cl_tumor", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(3.11))),
                        column(2,
                               sliderInput("v2_tumor", label = "V2 (L)", min = 1, max = 200, step = 0.1, value = c(23.6))),
                        column(2,
                               sliderInput("v3_tumor", label = "V3 (L)", min = 10, max = 500, step = 1, value = c(249))),
                        column(2,
                               sliderInput("q_tumor", label = "Q (L/h)", min = 1, max = 100, step = 0.1, value = c(39.7))),
                        column(2,
                               sliderInput("ka_tumor", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.44)))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("emax_tumor", label = "EMAX", min = 0.0001, max = 0.01, step = 0.0001, value = c(0.0069))),
                        column(3,
                               sliderInput("ec50_tumor", label = "EC50 (ng/mL)", min = 1, max = 200, step = 5, value = c(56.4*4.99/1.8))),
                        column(2,
                               sliderInput("lamda0_tumor", label = "Lambda0", min = 0.0001, max = 0.01, step = 0.0001, value = c(0.0041)))
                        # column(2,
                        #        sliderInput("w0_tumor", label = "W0 (mm^3)", min = 100, max = 100, step = 10, value = c(100))),
                        # column(2,
                        #        sliderInput("frac_tumor", label = "frac", min = 0.05, max = 0.05, step = 0.05, value = c(0.098)))
                    ),
                    absolutePanel(
                        top = 150, right = 120, width = 80, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>TGI<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "tgi", inline = T), "</code> %</strong>")
                        )
                    ),
                    # absolutePanel(
                    #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                    #     h3("ing")
                    # ),
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            
            
            tabPanel(
                strong("pERK"),
                fluidRow(
                    fluidRow(
                        column(4, 
                               # h3("Dose regimen"),
                               fluidRow(
                                   column(6),
                                   column(6,
                                          sliderInput("cycle_biomarker", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("amt_biomarker1", label = "dose (mg)", min = 0, max = 30, step = 0.5, value = c(4))
                                   ),
                                   column(6,
                                          sliderInput("amt_biomarker2", label = "dose interruption (mg)", min = 0, max = 30, step = 0.5, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("interval_biomarker1", label = "dose interval (h)", min = 0, max = 168, step = 12, value = c(24))
                                   ),
                                   column(6,
                                          sliderInput("interval_biomarker2", label = "interrupted interval (h)", min = 0, max = 168, step = 12, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("n_biomarker1", label = "dose times", min = 1, max = 28, step = 1, value = c(28))
                                   ),
                                   column(6,
                                          sliderInput("n_biomarker2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                                   )
                               )
                        ),
                        column(4, plotOutput("distPlotEfficacyPKPDtime")),
                        column(4, plotOutput("distPlotEfficacyPKPD"))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("cl_biomarker", label = "CL (L/h)", min = 1, max = 10, step = 0.01, value = c(3.11))),
                        column(3,
                               sliderInput("v2_biomarker", label = "V2 (L)", min = 1, max = 200, step = 0.1, value = c(23.6))),
                        column(3,
                               sliderInput("v3_biomarker", label = "V3 (L)", min = 1, max = 500, step = 1, value = c(249))),
                        column(3,
                               sliderInput("q_biomarker", label = "Q (L/h)", min = 1, max = 200, step = 0.1, value = c(39.7))),
                        column(3,
                               sliderInput("ka_biomarker", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.44)))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("gamma_biomarker", label = "Gamma", min = 0.1, max = 10, step = 0.1, value = c(1))),
                        column(3,
                               sliderInput("emax_biomarker", label = "EMAX", min = 0.1, max = 1, step = 0.05, value = c(1))),
                        column(3,
                               sliderInput("ec50_biomarker", label = "EC50 (ng/mL)", min = 1, max = 50, step = 1, value = c(11.1)))
                    ),
                    absolutePanel(
                        top = 200, right = 180, width = 130, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>Inhibition<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "inhi", inline = T), "</code> %</strong>")
                        )
                    ),
                    # absolutePanel(
                    #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                    #     h3("ing")
                    # ),
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            
            
            
            tabPanel(
                strong("PLT"),
                fluidRow(
                    fluidRow(
                        column(4, 
                               # h3("Dose regimen"),
                               fluidRow(
                                   column(6),
                                   column(6,
                                          sliderInput("cycle_plt", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("amt_plt1", label = "dose (mg)", min = 0, max = 30, step = 0.5, value = c(4))
                                   ),
                                   column(6,
                                          sliderInput("amt_plt2", label = "dose interruption (mg)", min = 0, max = 30, step = 0.5, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("interval_plt1", label = "dose interval (h)", min = 0, max = 168, step = 12, value = c(24))
                                   ),
                                   column(6,
                                          sliderInput("interval_plt2", label = "interrupted interval (h)", min = 0, max = 168, step = 12, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("n_plt1", label = "dose times", min = 1, max = 28, step = 1, value = c(28))
                                   ),
                                   column(6,
                                          sliderInput("n_plt2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                                   )
                               )
                        ),
                        
                        column(8, 
                               plotOutput("distPlotSafetyPKPD"))
                    ),
                    fluidRow(
                        column(2,
                               sliderInput("cl_plt", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(3.11))),
                        column(2,
                               sliderInput("v2_plt", label = "V2 (L)", min = 1, max = 200, step = 0.1, value = c(23.6))),
                        column(2,
                               sliderInput("v3_plt", label = "V3 (L)", min = 1, max = 500, step = 1, value = c(249))),
                        column(2,
                               sliderInput("q_plt", label = "Q (L/h)", min = 10, max = 100, step = 0.1, value = c(39.7))),
                        column(2,
                               sliderInput("ka_plt", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.44)))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("slope_plt", label = "Slope", min = 0.01, max = 10, step = 0.01, value = c(2.2))),
                        column(3,
                               sliderInput("circ0_plt", label = "Circ0 (*10^9)", min = 1, max = 500, step = 1, value = c(230))),
                        column(3,
                               sliderInput("gamma_plt", label = "Gamma", min = 0.001, max = 1, step = 0.001, value = c(0.15))),
                        column(3,
                               sliderInput("mtt_plt", label = "MTT (h)", min = 50, max = 200, step = 1, value = c(110)))
                    ),
                    absolutePanel(
                        top = 290, right = 30, width = 170, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>PLT nadir<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "nadir", inline = T), "</code> *10^9</strong>")
                        )
                    ),
                    
                    absolutePanel(
                        top = 180, right = 10, width = 210, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong> Gr3 AE: <sub>", "</sub>  <code style='color:#ec4c3c;background-color:#F8F9F9'>","25~<50", 
                                   "</code> *10^9</strong>")
                        )
                    ),
                    absolutePanel(
                        top = 210, right = 0, width = 190, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong> Gr4 AE: <sub>", "</sub>  <code style='color:#ec4c3c;background-color:#F8F9F9'>","<25", 
                                   "</code> *10^9</strong>")
                        )
                    ),
                    
                    # absolutePanel(
                    #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                    #     h3("ing")
                    # ),
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            
            tabPanel(
                strong("Population Based Simulation"),
                fluidRow(
                    fluidRow(
                        # h3("Dose regimen"),
                        column(3, 
                               br(),br(),br(),br(),br(),
                               sliderInput("cycle_pop", label = "cycles", min = 1, max = 12, step = 1, value = c(4))
                        ),
                        column(3,
                               sliderInput("amt1_pop", label = "dose1 (mg)", min = 0, max = 50, step = 0.5, value = c(14)),
                               sliderInput("amt2_pop", label = "dose2 (mg)", min = 0, max = 50, step = 0.5, value = c(0))
                        ),
                        column(3,
                               sliderInput("interval1_pop", label = "dose1 interval (h)", min = 0, max = 168, step = 12, value = c(24)),
                               sliderInput("interval2_pop", label = "dose2 interval (h)", min = 0, max = 168, step = 12, value = c(24))
                        ),
                        column(3,
                               sliderInput("n1_pop", label = "dose1 times", min = 1, max = 28*6, step = 1, value = c(2)),
                               sliderInput("n2_pop", label = "dose2 times", min = 1, max = 28*6, step = 1, value = c(5))
                        )
                    ),
                    
                    br(),
                    hr(),
                    
                    fluidRow(
                        column(4,
                               
                               br(),
                               fluidRow(
                                   column(6,
                                          sliderInput("cl_pkpop", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(3.11))),
                                   column(6,
                                          sliderInput("v2_pkpop", label = "V2 (L)", min = 1, max = 200, step = 0.1, value = c(23.6)))),
                               fluidRow(
                                   column(6,
                                          sliderInput("v3_pkpop", label = "V3 (L)", min = 1, max = 500, step = 1, value = c(249))),
                                   column(6,
                                          sliderInput("q_pkpop", label = "Q (L/h)", min = 10, max = 100, step = 0.1, value = c(39.7)))),
                               fluidRow(
                                   column(6,
                                          sliderInput("ka_pkpop", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.44))),
                                   column(6)),
                               br(),
                               
                               
                               hr(), 
                               
                               br(),  
                                fluidRow(
                                   column(6,
                                          sliderInput("slope_pltpop", label = "Slope (mL/ug)", min = 0.01, max = 10, step = 0.01, value = c(2.2))),
                                   column(6,
                                          sliderInput("gamma_pltpop", label = "Gamma", min = 0.001, max = 1, step = 0.001, value = c(0.187)))),
                               fluidRow(
                                   column(6,
                                          sliderInput("circ0_1_pltpop", label = "Circ0 1 (10^9)", min = 1, max = 500, step = 1, value = c(226))),
                                   column(6,
                                          sliderInput("circ0_2_pltpop", label = "Circ0 2 (10^9)", min = 1, max = 500, step = 1, value = c(148)))),
                               fluidRow(
                                   column(6,
                                          sliderInput("mtt_pltpop", label = "MTT (h)", min = 50, max = 200, step = 1, value = c(109))),
                                   column(6)),
                               br(),
                               
                               hr(),
                               
                               br(),
                               fluidRow(
                                   column(6,
                                          sliderInput("emax_tumorpop", label = "EMAX", min = 0.0001, max = 0.01, step = 0.0001, value = c(0.0069))),
                                   column(6,
                                          sliderInput("ec50_tumorpop", label = "EC50 (ng/mL)", min = 1, max = 200, step = 5, value = c(56.4*4.99/1.8)))),
                               fluidRow(
                                   column(6,
                                          sliderInput("lamda0_tumorpop", label = "Lambda0", min = 0.0001, max = 0.01, step = 0.0001, value = c(0.0041))),
                                   column(6))
                        ),
                        
                        column(8,
                               
                               plotOutput("distPlotPKpop"),
                               plotOutput("distPlotPLTpop"),
                               plotOutput("distPlotTumorpop"))
                    ),
                    
                    
                    br(),
                    br(),
                    hr(),
                    
                    plotOutput("distPlotPKPDpop", 
                               # width = 16*200, 
                               height = 12*100),
                    
                    # absolutePanel(
                    #     top = 400, right = 30, width = 170, height = 10, draggable = TRUE,
                    #     HTML(
                    #         paste0("<strong>PLTpop nadir<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                    #                textOutput(outputId = "nadirpop", inline = T), "</code> *10^9</strong>")
                    #     )
                    # ),
                    # 
                    # absolutePanel(
                    #     top = 150, right = 40, width = 210, height = 10, draggable = TRUE,
                    #     HTML(
                    #         paste0("<strong> Grade 3 PLT ↓: <sub>", "</sub>  <code style='color:#ec4c3c;background-color:#F8F9F9'>","25~<50", 
                    #                "</code> *10^9</strong>")
                    #     )
                    # ),
                    # absolutePanel(
                    #     top = 180, right = 40, width = 190, height = 10, draggable = TRUE,
                    #     HTML(
                    #         paste0("<strong> Grade 4 PLT ↓: <sub>", "</sub>  <code style='color:#ec4c3c;background-color:#F8F9F9'>","<25", 
                    #                "</code> *10^9</strong>")
                    #     )
                    # ),
                    
                    
                    absolutePanel(
                        top = 780, right = 30, width = 170, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>PLT nadir1<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "nadirpop1", inline = T), "</code> *10^9</strong>")
                        )
                    ),
                    absolutePanel(
                        top = 810, right = 30, width = 170, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>PLT nadir2<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "nadirpop2", inline = T), "</code> *10^9</strong>")
                        )
                    ),
                    
                    
                    absolutePanel(
                        top = 1200, right = 130, width = 80, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>TGI<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "tgipop1", inline = T), "</code> %</strong>")
                        )
                    ),
                    
                    # absolutePanel(
                    #     top = 2470, right = 130, width = 80, height = 10, draggable = TRUE,
                    #     HTML(
                    #         paste0("<strong>TGI<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                    #                textOutput(outputId = "tgipop2", inline = T), "</code> %</strong>")
                    #     )
                    # ),
                    
                    
                    absolutePanel(
                        top = 290, right = 400, width = 160, height = 10, draggable = FALSE,
                        titlePanel(h5(strong("Subject")))
                    ),
                    absolutePanel(
                        top = 290, right = 325, width = 160, height = 10, draggable = FALSE,
                        sliderInput("idv_pop", label = NULL, min = 0, max = 1000, step = 10, value = 300)
                    ),
                    
                    absolutePanel(
                        top = 290, right = 115, width = 160, height = 10, draggable = FALSE,
                        titlePanel(h5(strong("Proportion")))
                    ),
                    absolutePanel(
                        top = 290, right = 35, width = 160, height = 10, draggable = FALSE,
                        sliderInput("p_pop", label = NULL, min = 0, max = 1, step = 0.05, value = c(0.1,0.9))
                    ),
                    
                    
                    
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            
            footer = h5(HTML("dMed Copyright 2020 : 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> E </strong>arly 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> C </strong>linical 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> D </strong>evelopment"), align = "right")
        )
    )
)

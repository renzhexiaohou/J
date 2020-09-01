#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(mrgsolve)
library(deSolve)
library(PKNCA)

code_pk <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double ke = CLi/V2i;

        $INIT GUT = 0, CENT = 0, PERH = 0
        
        $ODE 
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Q/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2*CENT - Qi/V3i*PERH;

        $OMEGA @name IIV @labels ETA_CL ETA_V2 ETA_V3 ETA_Q
        0 0 0 0
        
        $SIGMA 0
        
        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
    '
mod_pk <- mcode("pk_model", code_pk)


code_tumor <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, LAMDA0 = 0.005, EMAX = 0.005, EC50 = 50, W0 = 100, FRAC = 0

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double ke = CLi/V2i;
        double EMAXi = EMAX*exp(ETA_EMAX);
        double EC50i = EC50*exp(ETA_EC50);
        
        
        $INIT GUT = 0, CENT = 0, PERH = 0, TUMOR1 = 100, TUMOR2 = 0
        
        $ODE 
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Q/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2*CENT - Qi/V3i*PERH;
        
        double EFFECT = EMAXi*PKCONC/(EC50i + PKCONC);
        
        dxdt_TUMOR1 = LAMDA0*TUMOR1 - EFFECT*TUMOR1;
        dxdt_TUMOR2 = LAMDA0*TUMOR2;
        
        
        $OMEGA @name IIV @labels ETA_CL ETA_V2 ETA_V3 ETA_Q ETA_EMAX ETA_EC50
        0 0 0 0 0 0
        
        $SIGMA 0
        
        $TABLE
        capture TUMOR = TUMOR1 + TUMOR2;

        $CAPTURE PKCONC TUMOR
    '
mod_tumor <- mcode("tumor_growth_inhibition", code_tumor)



code_biomarker <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, EMAX = 1, EC50 = 11.1, GAMMA = 1, RECL = 0, REC50 = 0
        
        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double ke = CLi/V2i;
        double EMAXi = EMAX*exp(ETA_EMAX);
        double EC50i = EC50*exp(ETA_EC50);
        double GAMMAi = GAMMA;
        
        
        $INIT GUT = 0, CENT = 0, PERH = 0
        
        $ODE 
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Q/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2*CENT - Qi/V3i*PERH;

        $OMEGA @name IIV @labels ETA_CL ETA_V2 ETA_V3 ETA_Q ETA_EMAX ETA_EC50
        0 0 0 0 0 0
        
        $SIGMA 0
        
        $TABLE
        capture EFFECT = EMAXi*pow(PKCONC, GAMMAi) / (EC50i + pow(PKCONC, GAMMAi));

        $CAPTURE PKCONC EFFECT
    '
mod_biomarker <- mcode("efficacy_biomarker", code_biomarker)



code_plt <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, SLOPE = 0.0007, CIRC0 = 226, GAMMA = 0.0566, MTT = 59.8

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double ke = CLi/V2i;
        double KTR = 4/MTT;
        double KPROL = KTR;
        double KCIRC = KTR;
        double CIRC0i = CIRC0;


        $INIT GUT = 0, CENT = 0, PERH = 0, PROL = 226, TRANSIT1 = 226, TRANSIT2 = 226, TRANSIT3 = 226, CIRC = 226

        $ODE
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Q/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2*CENT - Qi/V3i*PERH;

        double EFFECT = SLOPE*PKCONC;

        dxdt_PROL = KPROL*PROL*(1-EFFECT)*pow(CIRC0/CIRC, GAMMA) - KTR*PROL;
        
        dxdt_TRANSIT1 = KTR*PROL - KTR*TRANSIT1;
        dxdt_TRANSIT2 = KTR*TRANSIT1 - KTR*TRANSIT2;
        dxdt_TRANSIT3 = KTR*TRANSIT2 - KTR*TRANSIT3;
        
        dxdt_CIRC = KTR*TRANSIT3 - KCIRC*CIRC;


        $OMEGA @name IIV @labels ETA_CL ETA_V2 ETA_V3 ETA_Q
        0 0 0 0

        $SIGMA 0

        $TABLE
        capture PLT = CIRC;

        $CAPTURE PKCONC PLT
    '
mod_plt <- mcode("safety_biomarker_plt", code_plt)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    dat_pk <- reactive({
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        cl <-input$cl_pk
        v2 <- input$v2_pk
        v3 <- input$v3_pk
        q <- input$q_pk
        ka <- input$ka_pk
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V2 = v2, 
                          V3 = v3, 
                          Q = q, 
                          KA = ka,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V2 = v2, 
                          V3 = v3, 
                          Q = q, 
                          KA = ka,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V2 = v2, 
                         V3 = v3, 
                         Q = q, 
                         KA = ka,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V2 = rnorm(1, 0, 0),
                         ETA_V3 = rnorm(1, 0, 0),
                         ETA_Q = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        # union_all(dos1, dos2) %>% arrange(ID, time, addl, ii)
        union_all(dos1, dos2) %>%
            # union_all(dos) %>% 
            arrange(ID, time, addl, ii)
    })
    
    output$distPlotPK <- renderPlot({
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        # mod_pk %>%
        #     data_set(dat_pk()) %>% 
        #     mrgsim(end = cyc * (tau1*n1 + tau2*n2),
        #            delta = 0.05) %>% 
        #     plot(PKCONC~time)
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        # browser()
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        ggplot(data = pk,
               aes(x = time/24+1, y = PKCONC)) +
            geom_line(size = 0.8) +
            # geom_hline(aes(yintercept = 11.14), colour = "grey", linetype = "dashed", size = 1, alpha = 0.8) +
            # geom_hline(aes(yintercept = 44.55), colour = "grey", linetype = "dashed", size = 1, alpha = 0.8) +
            # geom_hline(aes(yintercept = 100.24), colour = "grey", linetype = "dashed", size = 1, alpha = 0.8) +
            # geom_hline(aes(yintercept = 24), colour = "red", linetype = "dashed", size = 1, alpha = 0.8) +
            # geom_hline(aes(yintercept = 337), colour = "red", linetype = "dashed", size = 1, alpha = 0.8) +
            # geom_text(data = concrange.txt, aes(label = name, fontface = 2), colour = c("red", "#0080ff"), size = 5) +
            # scale_x_continuous(breaks = seq(0, 168*100, 168)) +
            scale_x_continuous(breaks = seq(1, 7*100, 7)) +
            scale_y_continuous(breaks = seq(0, 500, 20), limits = c(0, 160)) +
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 20),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
            xlab("Nominal Day") + ylab("Concentration (ng/mL)")
        
        
    })
    
    
    output$halflife <- renderText({
        # browser()
        round(0.693/(input$cl_pk/input$v3_pk), digits = 1)
    })
    
    
    reacnca <- reactive({
        # PKNCA packages code as followed
        # (since ncappc packages can't be used in linux sys)
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        
        
        nca <- PKNCAconc(pk, conc~time|subject) %>%
            PKNCAdata(.,
                      intervals=data.frame(
                          start= cyc*(tau1*(n1) + tau2*(n2)) - tau1 - tau2*n2,
                          end= cyc*(tau1*(n1) + tau2*(n2)) - tau2*n2,
                          cmax=TRUE,
                          cmin=TRUE,
                          tmax=TRUE,
                          # aucinf.obs=TRUE,
                          auclast=TRUE
                      )
            ) %>%
            pk.nca() %>%
            as.data.frame(.$result)
    })
    
    output$pkcmax1 <- renderText({
        # browser()
        reacnca() %>% filter(PPTESTCD == "cmax") %>% .$PPORRES %>% 
            round(digits = 1)
        
    })
    output$pkcmin1 <- renderText({
        # browser()
        reacnca() %>% filter(PPTESTCD == "cmin") %>% .$PPORRES %>% 
            round(digits = 1)
    })
    output$pkauc1 <- renderText({
        # browser()
        reacnca() %>% filter(PPTESTCD == "auclast") %>% .$PPORRES %>% 
            round(digits = 1)
    })
    
    
    reacnca2 <- reactive({
        # PKNCA packages code as followed
        # (since ncappc packages can't be used in linux sys)
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        
        t_start <- input$t_start
        t_end <- input$t_end
        
        nca2 <- PKNCAconc(pk, conc~time|subject) %>%
            PKNCAdata(.,
                      intervals=data.frame(
                          start= t_start,
                          end= t_end,
                          cmax=TRUE,
                          cmin=TRUE,
                          tmax=TRUE,
                          # aucinf.obs=TRUE,
                          auclast=TRUE
                      )
            ) %>%
            pk.nca() %>%
            as.data.frame(.$result)
    })
    
    output$pkcmax2 <- renderText({
        # browser()
        reacnca2() %>% filter(PPTESTCD == "cmax") %>% .$PPORRES %>% 
            round(digits = 1)
        
    })
    output$pkcmin2 <- renderText({
        # browser()
        reacnca2() %>% filter(PPTESTCD == "cmin") %>% .$PPORRES %>% 
            round(digits = 1)
    })
    output$pkauc2 <- renderText({
        # browser()
        reacnca2() %>% filter(PPTESTCD == "auclast") %>% .$PPORRES %>% 
            round(digits = 1)
    })
    
    
    
    
    output$downloadPKData <- downloadHandler(
        filename = function() {
            paste("PKdata-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(data, file)
        }
    )
    
    dat_tumor <- reactive({
        
        cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        cl <-input$cl_tumor
        v2 <- input$v2_tumor
        v3 <- input$v3_tumor
        q <- input$q_tumor
        ka <- input$ka_tumor
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V2 = v2,
                          V3 = v3,
                          Q = q,
                          KA = ka,
                          LAMDA0 = input$lamda0_tumor,
                          EMAX = input$emax_tumor,
                          EC50 = input$ec50_tumor,
                          # W0 = input$w0_tumor,
                          # FRAC = input$frac_tumor,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V2 = v2,
                          V3 = v3,
                          Q = q,
                          KA = ka,
                          LAMDA0 = input$lamda0_tumor,
                          EMAX = input$emax_tumor,
                          EC50 = input$ec50_tumor,
                          # W0 = input$w0_tumor,
                          # FRAC = input$frac_tumor,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V2 = v2,
                         V3 = v3,
                         Q = q,
                         KA = ka,
                         LAMDA0 = input$lamda0_tumor,
                         EMAX = input$emax_tumor,
                         EC50 = input$ec50_tumor,
                         # W0 = input$w0_tumor,
                         # FRAC = input$frac_tumor,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V2 = rnorm(1, 0, 0),
                         ETA_V3 = rnorm(1, 0, 0),
                         ETA_Q = rnorm(1, 0, 0),
                         ETA_EMAX = rnorm(1, 0, 0),
                         ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        union_all(dos1, dos2) %>% 
            union_all(dos) %>%
            arrange(ID, time, addl, ii)
        
    })
    
    output$distPlotTumorPK <- renderPlot({
        
        cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC~time)
    })
    
    output$distPlotTumorPD <- renderPlot({
        cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(TUMOR~time)
    })
    
    output$tgi <- renderText({
        # browser()
        
        cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        # browser()
        
        out <- mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        t0 <- tmp %>% filter(ID == 2) %>% .$TUMOR %>% max()
        t <- tmp %>% filter(ID == 1) %>% .$TUMOR %>% max()
        
        # browser()
        
        round(abs(t - t0)/t0 * 100, digits = 0)
        
        # reacauc()
    })
    
    
    
    dat_biomarker <- reactive({
        
        cyc <- input$cycle_biomarker
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        
        cl <-input$cl_biomarker
        v2 <- input$v2_biomarker
        v3 <- input$v3_biomarker
        q <- input$q_biomarker
        ka <- input$ka_biomarker
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V2 = v2, 
                          V3 = v3, 
                          Q = q, 
                          KA = ka,
                          EMAX = input$emax_biomarker,
                          EC50 = input$ec50_biomarker,
                          GAMMA = input$gamma_biomarker,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V2 = v2, 
                          V3 = v3, 
                          Q = q, 
                          KA = ka,
                          EMAX = input$emax_biomarker,
                          EC50 = input$ec50_biomarker,
                          GAMMA = input$gamma_biomarker,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V2 = v2, 
                         V3 = v3, 
                         Q = q, 
                         KA = ka,
                         EMAX = input$emax_biomarker,
                         EC50 = input$ec50_biomarker,
                         GAMMA = input$gamma_biomarker,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V2 = rnorm(1, 0, 0),
                         ETA_V3 = rnorm(1, 0, 0),
                         ETA_Q = rnorm(1, 0, 0),
                         ETA_EMAX = rnorm(1, 0, 0),
                         ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        union_all(dos1, dos2) %>% 
            union_all(dos) %>%
            arrange(ID, time, addl, ii)
    })
    
    output$distPlotEfficacyPKPDtime <- renderPlot({
        cyc <- input$cycle_biomarker
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        
        mod_biomarker %>%
            data_set(dat_biomarker()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC + EFFECT~time)
    })
    
    output$distPlotEfficacyPKPD <- renderPlot({
        
        cyc <- input$cycle_biomarker
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        
        mod_biomarker %>%
            data_set(dat_biomarker()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(EFFECT~PKCONC)
    })
    
    output$inhi <- renderText({
        
        cyc <- input$cycle_biomarker
        cyc <- 1
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        out <- mod_biomarker %>%
            data_set(dat_biomarker()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        # browser()
        
        eff <- tmp %>% filter(ID == 1, 
                              time < cyc * (tau1*n1 + tau2*(n2-1)),
                              time > cyc * (tau1*(n1-1) + tau2*(n2-1))) %>% .$EFFECT %>% min()
        
        round(abs(eff) * 100, digits = 0)
    })
    
    
    
    
    dat_plt <- reactive({
        
        cyc <- input$cycle_plt
        
        amt1 <- input$amt_plt1
        tau1 <- input$interval_plt1
        n1 <- input$n_plt1
        
        amt2 <- input$amt_plt2
        tau2 <- input$interval_plt2
        n2 <- input$n_plt2
        
        cl <-input$cl_plt
        v2 <- input$v2_plt
        v3 <- input$v3_plt
        q <- input$q_plt
        ka <- input$ka_plt
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V2 = v2, 
                          V3 = v3, 
                          Q = q, 
                          KA = ka,
                          SLOPE = input$slope_plt,
                          CIRC0 = input$circ0_plt,
                          GAMMA = input$gamma_plt,
                          MTT = input$mtt_plt,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V2 = v2, 
                          V3 = v3, 
                          Q = q, 
                          KA = ka,
                          SLOPE = input$slope_plt,
                          CIRC0 = input$circ0_plt,
                          GAMMA = input$gamma_plt,
                          MTT = input$mtt_plt,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V2 = rnorm(1, 0, 0),
                          ETA_V3 = rnorm(1, 0, 0),
                          ETA_Q = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V2 = v2, 
                         V3 = v3, 
                         Q = q, 
                         KA = ka,
                         SLOPE = input$slope_plt,
                         CIRC0 = input$circ0_plt,
                         GAMMA = input$gamma_plt,
                         MTT = input$mtt_plt,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V2 = rnorm(1, 0, 0),
                         ETA_V3 = rnorm(1, 0, 0),
                         ETA_Q = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        
        union_all(dos1, dos2) %>% 
            union_all(dos) %>%
            arrange(ID, time, addl, ii)
        
    })
    
    output$distPlotSafetyPKPD <- renderPlot({
        cyc <- input$cycle_plt
        
        amt1 <- input$amt_plt1
        tau1 <- input$interval_plt1
        n1 <- input$n_plt1
        
        amt2 <- input$amt_plt2
        tau2 <- input$interval_plt2
        n2 <- input$n_plt2
        
        mod_plt %>%
            data_set(dat_plt()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC + PLT~time)
    })
    
    output$nadir <- renderText({
        # browser()
        cyc <- input$cycle_plt
        
        amt1 <- input$amt_plt1
        tau1 <- input$interval_plt1
        n1 <- input$n_plt1
        
        amt2 <- input$amt_plt2
        tau2 <- input$interval_plt2
        n2 <- input$n_plt2
        
        # browser()
        
        out <- mod_plt %>%
            data_set(dat_plt()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        nadir <- tmp %>% filter(ID == 1) %>% .$PLT %>% min()
        
        # browser()
        
        round(nadir, digits = 2)
        
        # reacauc()
    })
    
    
})


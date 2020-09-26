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
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;

        $OMEGA @name IIV @labels ETA_CL ETA_V2 ETA_V3 ETA_Q
        0 0 0 0
        
        $SIGMA 0
        
        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
    '
mod_pk <- mcode("pk_model", code_pk)


code_pkpop <- '
        $PARAM KA = 0.44, CL = 3.11, V2 = 23.6, V3 = 249, Q = 39.7

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double KAi = KA*exp(ETA_KA);
        double ke = CLi/V2i;

        $INIT GUT = 0, CENT = 0, PERH = 0
        
        $ODE 
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT  = -KAi*GUT;
        dxdt_CENT = KAi*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;

        
        $OMEGA @annotated @block 
        ETA_CL:     0.09 : Clearance 0.097
        ETA_V2:     0   0.09 : central volume 0.44
        ETA_V3:     0   0   0.09 : peripheral volume 0.31
        ETA_Q:      0   0   0   0 : inter clearance
        ETA_KA:     0   0   0   0   0   : absorb rate 0.141
        
        $SIGMA 0
        
        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
    '
mod_pkpop <- mcode("pkpop_model", code_pkpop)


code_tumor <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, LAMDA0 = 0.005, EMAX = 0.005, EC50 = 50, W0 = 100, FRAC = 0.1

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double KAi = KA*exp(ETA_KA);
        double ke = CLi/V2i;
        double EMAXi = EMAX*exp(ETA_EMAX);
        double EC50i = EC50*exp(ETA_EC50);
        double LAMDA0i = LAMDA0*exp(ETA_LAMDA0);
        
        
        $INIT GUT = 0, CENT = 0, PERH = 0, TUMOR1 = 100, TUMOR2 = 9.8
        
        $ODE 
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KAi*GUT;
        dxdt_CENT = KAi*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;
        
        double EFFECT = EMAXi*PKCONC/(EC50i + PKCONC);
        
        dxdt_TUMOR1 = LAMDA0i*TUMOR1 - EFFECT*TUMOR1;
        dxdt_TUMOR2 = LAMDA0i*TUMOR2;
        
        
        $OMEGA @annotated @block 
        ETA_CL:     0    : Clearance 0.097
        ETA_V2:     0   0    : central volume 0.44
        ETA_V3:     0   0   0    : peripheral volume 0.31
        ETA_Q:      0   0   0   0       : inter clearance
        ETA_KA:     0   0   0   0   0       : absorb rate 0.141
        ETA_EMAX:   0   0   0   0   0   0       : EMAX
        ETA_EC50:   0   0   0   0   0   0   0   : EC50
        ETA_LAMDA0: 0   0   0   0   0   0   0   0      : LAMDA0
        
        $SIGMA 0
        
        $TABLE
        capture TUMOR = TUMOR1 + TUMOR2;

        $CAPTURE PKCONC TUMOR
    '
mod_tumor <- mcode("tumor_growth_inhibition", code_tumor)


code_tumorpop <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, LAMDA0 = 0.005, EMAX = 0.005, EC50 = 50, W0 = 100, FRAC = 0.1

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double KAi = KA*exp(ETA_KA);
        double ke = CLi/V2i;
        double EMAXi = EMAX*exp(ETA_EMAX);
        double EC50i = EC50*exp(ETA_EC50);
        double LAMDA0i = LAMDA0*exp(ETA_LAMDA0);
        
        
        $INIT GUT = 0, CENT = 0, PERH = 0, TUMOR1 = 100, TUMOR2 = 9.8
        
        $ODE 
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KAi*GUT;
        dxdt_CENT = KAi*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;
        
        double EFFECT = EMAXi*PKCONC/(EC50i + PKCONC);
        
        dxdt_TUMOR1 = LAMDA0i*TUMOR1 - EFFECT*TUMOR1;
        dxdt_TUMOR2 = LAMDA0i*TUMOR2;
        
        $OMEGA @annotated @block 
        ETA_CL:     0.09    : Clearance 0.097
        ETA_V2:     0   0.09    : central volume 0.44
        ETA_V3:     0   0   0.09    : peripheral volume 0.31
        ETA_Q:      0   0   0   0       : inter clearance
        ETA_KA:     0   0   0   0   0       : absorb rate 0.141
        ETA_EMAX:   0   0   0   0   0   0       : EMAX
        ETA_EC50:   0   0   0   0   0   0   0.145   : EC50
        ETA_LAMDA0: 0   0   0   0   0   0   0   0      : LAMDA0
        
        $SIGMA 0
        
        $TABLE
        capture TUMOR = TUMOR1 + TUMOR2;

        $CAPTURE PKCONC TUMOR
    '
mod_tumorpop <- mcode("tumor_growth_inhibitionpop", code_tumorpop)


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
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;

        $OMEGA @name IIV @labels ETA_CL ETA_V2 ETA_V3 ETA_Q ETA_EMAX ETA_EC50
        0 0 0 0 0 0
        
        $SIGMA 0
        
        $TABLE
        capture EFFECT = EMAXi*pow(PKCONC, GAMMAi) / (EC50i + pow(PKCONC, GAMMAi));

        $CAPTURE PKCONC EFFECT
    '
mod_biomarker <- mcode("efficacy_biomarker", code_biomarker)



code_plt <- '
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, SLOPE = 2.4, CIRC0 = 230, GAMMA = 0.15, MTT = 100

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


        $INIT GUT = 0, CENT = 0, PERH = 0, PROL = 230, TRANSIT1 = 230, TRANSIT2 = 230, TRANSIT3 = 230, CIRC = 230

        $ODE
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = KA*GUT - CLi/V2i*CENT - Q/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2*CENT - Qi/V3i*PERH;

        double EFFECT = SLOPE*PKCONC/1000;

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


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}



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
        
        union_all(dos1, dos2) %>%
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
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        
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
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        
        mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC~time)
    })
    
    output$distPlotTumorPD <- renderPlot({
        cyc <- input$cycle_tumor
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        
        mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(TUMOR~time)
    })
    
    output$tgi <- renderText({
        # browser()
        
        cyc <- input$cycle_tumor
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        amt2 <- input$amt_tumor2
        tau2 <- input$interval_tumor2
        n2 <- input$n_tumor2
        
        # browser()
        
        out <- mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        t0 <- tmp %>% filter(ID == 2) %>% .$TUMOR %>% max()
        t <- tmp %>% filter(ID == 1) %>% .$TUMOR %>% max()
        
        # browser()
        
        round(abs(t - t0)/(t0-100) * 100, digits = 0)
        
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
    
    
    
    
    PK <- reactive({
        
        cyc <- input$cycle_pop
        
        amt1 <- input$amt1_pop
        tau1 <- input$interval1_pop
        n1 <- input$n1_pop
        
        amt2 <- input$amt2_pop
        tau2 <- input$interval2_pop
        n2 <- input$n2_pop
        
        
        idv <- input$idv_pop
        
        
        cl <-input$cl_pkpop
        v2 <- input$v2_pkpop
        v3 <- input$v3_pkpop
        q <- input$q_pkpop
        ka <- input$ka_pkpop
        
        
        dos1 <- expand.grid(cmt = 1,
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
                            ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt)
        
        dos2 <- expand.grid(cmt = 1,
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
                            ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt)
        
        dos <- union_all(dos1, dos2) %>%
            mutate(evid = 1) %>% 
            arrange(ID, time, addl, ii)
        
        out <- mod_pkpop %>%
            data_set(dos) %>%
            mrgsim(set.seed(2020),
                   end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.5)
        
        # browser()
        PK <- out@data %>%
            select(ID, time, PKCONC) %>% 
            mutate(subject = ID, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        PK
    })
    
    
    PKSUM <- reactive({
        tmp <- PK() %>% group_by(time) %>% mutate(
            p_upper = quantile(conc, probs = c(input$p_pop[2]), na.rm = TRUE),
            p_lower = quantile(conc, probs = c(input$p_pop[1]), na.rm = TRUE),
            p_50 = quantile(conc, probs = c(0.5), na.rm = TRUE)) %>% 
            select(., subject, time, conc, p_upper, p_lower, p_50) %>% 
            ungroup() %>% 
            filter(subject == 1)
    })
    
    output$distPlotPKpop <- renderPlot({
        
        PK <- PK()
        PKSUM <- PKSUM()
        # browser()
        
        ggplot(PK, aes(x = time, y = conc, group=subject)) +
            scale_x_continuous(breaks = seq(0, 168*4*12, 168)) +
            scale_y_continuous(breaks = c(0, 25, 50, 100, 200, 300, 400, 500, 600), limits = c(0, 350)) +
            # geom_line(alpha=0.2) +
            geom_ribbon(data=PKSUM, aes(x=time, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.3, fill="#566573") +
            geom_line(data=PKSUM, aes(x=time, y=p_50), linetype="solid", size = 0.8, colour="black") +
            
            geom_hline(aes(yintercept = 25), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            geom_hline(aes(yintercept = 50), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            geom_hline(aes(yintercept = 100), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            geom_hline(aes(yintercept = 200), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 21),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                  axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            xlab("Time (h)") + ylab("Plasma concentration (ng/mL)")
    })
    
    
    
    PLT <- reactive({
        
        cyc <- input$cycle_pop
        
        amt1 <- input$amt1_pop
        tau1 <- input$interval1_pop
        n1 <- input$n1_pop
        
        amt2 <- input$amt2_pop
        tau2 <- input$interval2_pop
        n2 <- input$n2_pop
        
        idv <- input$idv_pop
        
        
        cl <-input$cl_pkpop
        v2 <- input$v2_pkpop
        v3 <- input$v3_pkpop
        q <- input$q_pkpop
        ka <- input$ka_pkpop
        
        slope <- input$slope_pltpop
        circ0_1 <- input$circ0_1_pltpop
        circ0_2 <- input$circ0_2_pltpop
        gamma <- input$gamma_pltpop
        mtt <- input$mtt_pltpop
        
        
        dos1_1 <- expand.grid(cmt = 1,
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
                          SLOPE = slope,
                          CIRC0 = circ0_1,
                          GAMMA = gamma,
                          MTT = mtt,
                          ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt, evid = 1) %>% 
            arrange(ID, time, addl, ii)
        
        dos2_1 <- expand.grid(cmt = 1,
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
                          SLOPE = slope,
                          CIRC0 = circ0_1,
                          GAMMA = gamma,
                          MTT = mtt,
                          ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt, evid = 1) %>% 
            arrange(ID, time, addl, ii)
        
        dos_1 <- union_all(dos1_1, dos2_1) %>% 
            arrange(ID, time, addl, ii)
        
        
        
        code_pltpop1 <- paste0('
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, SLOPE = 2.4, CIRC0 = 230, GAMMA = 0.15, MTT = 100

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double KAi = KA*exp(ETA_KA);
        double ke = CLi/V2i;
        double SLOPEi = SLOPE*exp(ETA_SLOPE);
        double CIRC0i = CIRC0*exp(ETA_CIRC0);
        double GAMMAi = GAMMA*exp(ETA_GAMMA);
        double MTTi = MTT*exp(ETA_MTT);
        double KTRi = 4/MTTi;
        double KPROLi = KTRi;
        double KCIRCi = KTRi;


        $INIT GUT = 0, CENT = 0, PERH = 0,','PROL =', circ0_1, ',TRANSIT1 = ', circ0_1, ',TRANSIT2 = ', circ0_1, ',TRANSIT3 = ', circ0_1, ',CIRC = ', circ0_1,
        '
        $ODE
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KAi*GUT;
        dxdt_CENT = KAi*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;

        double EFFECT = SLOPEi*PKCONC/1000;

        dxdt_PROL = KPROLi*PROL*(1-EFFECT)*pow(CIRC0i/CIRC, GAMMAi) - KTRi*PROL;
        
        dxdt_TRANSIT1 = KTRi*PROL - KTRi*TRANSIT1;
        dxdt_TRANSIT2 = KTRi*TRANSIT1 - KTRi*TRANSIT2;
        dxdt_TRANSIT3 = KTRi*TRANSIT2 - KTRi*TRANSIT3;
        
        dxdt_CIRC = KTRi*TRANSIT3 - KCIRCi*CIRC;

        $OMEGA @annotated @block 
        ETA_CL:     0.09 : Clearance 0.097
        ETA_V2:     0   0.09 : central volume 0.44
        ETA_V3:     0   0   0.09 : peripheral volume 0.31
        ETA_Q:      0   0   0   0 : inter clearance
        ETA_KA:     0   0   0   0   0   : absorb rate 0.141
        ETA_SLOPE:  0   0   0   0   0   0.25   : effec slope
        ETA_CIRC0:  0   0   0   0   0   0   0 : circulation plt
        ETA_GAMMA:  0   0   0   0   0   0   0   0.18 : negative feedback
        ETA_MTT:    0   0   0   0   0   0   0   0   0 : mean time

        $SIGMA 0

        $TABLE
        capture PLT = CIRC;

        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
        PLT: PLT (10^9)
        
    ')
        mod_pltpop1 <- mcode("safety_biomarker_pltpop1", code_pltpop1)
        
        
        
        dos1_2 <- expand.grid(cmt = 1,
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
                            SLOPE = slope,
                            CIRC0 = circ0_2,
                            GAMMA = gamma,
                            MTT = mtt,
                            ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt, evid = 1) %>% 
            arrange(ID, time, addl, ii)
        
        dos2_2 <- expand.grid(cmt = 1,
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
                            SLOPE = slope,
                            CIRC0 = circ0_2,
                            GAMMA = gamma,
                            MTT = mtt,
                            ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt, evid = 1) %>% 
            arrange(ID, time, addl, ii)
        
        dos_2 <- union_all(dos1_2, dos2_2) %>% 
            arrange(ID, time, addl, ii)
        
        
        code_pltpop2 <- paste0('
        $PARAM KA = 0.32, CL = 2.91, V2 = 19.1, V3 = 259, Q = 29, SLOPE = 2.4, CIRC0 = 230, GAMMA = 0.15, MTT = 100

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double V2i = V2*exp(ETA_V2);
        double V3i = V3*exp(ETA_V3);
        double Qi = Q*exp(ETA_Q);
        double KAi = KA*exp(ETA_KA);
        double ke = CLi/V2i;
        double SLOPEi = SLOPE*exp(ETA_SLOPE);
        double CIRC0i = CIRC0*exp(ETA_CIRC0);
        double GAMMAi = GAMMA*exp(ETA_GAMMA);
        double MTTi = MTT*exp(ETA_MTT);
        double KTRi = 4/MTTi;
        double KPROLi = KTRi;
        double KCIRCi = KTRi;


        $INIT GUT = 0, CENT = 0, PERH = 0,','PROL =', circ0_2, ',TRANSIT1 = ', circ0_2, ',TRANSIT2 = ', circ0_2, ',TRANSIT3 = ', circ0_2, ',CIRC = ', circ0_2,
                               '
        $ODE
        double PKCONC = CENT/V2i*1000;
        dxdt_GUT = -KAi*GUT;
        dxdt_CENT = KAi*GUT - CLi/V2i*CENT - Qi/V2i*CENT + Qi/V3i*PERH;
        dxdt_PERH = Qi/V2i*CENT - Qi/V3i*PERH;

        double EFFECT = SLOPEi*PKCONC/1000;

        dxdt_PROL = KPROLi*PROL*(1-EFFECT)*pow(CIRC0i/CIRC, GAMMAi) - KTRi*PROL;
        
        dxdt_TRANSIT1 = KTRi*PROL - KTRi*TRANSIT1;
        dxdt_TRANSIT2 = KTRi*TRANSIT1 - KTRi*TRANSIT2;
        dxdt_TRANSIT3 = KTRi*TRANSIT2 - KTRi*TRANSIT3;
        
        dxdt_CIRC = KTRi*TRANSIT3 - KCIRCi*CIRC;

        $OMEGA @annotated @block 
        ETA_CL:     0.09 : Clearance 0.097
        ETA_V2:     0   0.09 : central volume 0.44
        ETA_V3:     0   0   0.09 : peripheral volume 0.31
        ETA_Q:      0   0   0   0 : inter clearance
        ETA_KA:     0   0   0   0   0   : absorb rate 0.141
        ETA_SLOPE:  0   0   0   0   0   0.25   : effec slope
        ETA_CIRC0:  0   0   0   0   0   0   0 : circulation plt
        ETA_GAMMA:  0   0   0   0   0   0   0   0.18 : negative feedback
        ETA_MTT:    0   0   0   0   0   0   0   0   0 : mean time

        $SIGMA 0

        $TABLE
        capture PLT = CIRC;

        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
        PLT: PLT (10^9)
        
    ')
        mod_pltpop2 <- mcode("safety_biomarker_pltpop2", code_pltpop2)
        
        
        
        out1 <- mod_pltpop1 %>%
            data_set(dos_1) %>%
            mrgsim(set.seed(2020),
                   end = cyc * (tau1*n1 + tau2*n2),
                   delta = 24)
        
        out2 <- mod_pltpop2 %>%
            data_set(dos_2) %>%
            mrgsim(set.seed(2020),
                   end = cyc * (tau1*n1 + tau2*n2),
                   delta = 24)
        
        # browser()
        PLT1 <- out1@data %>%
            select(ID, time, PKCONC, PLT) %>% 
            mutate(subject = ID, time = time, conc = PKCONC, plt = PLT, circ0 = circ0_1) %>%
            distinct()
        
        PLT2 <- out2@data %>%
            select(ID, time, PKCONC, PLT) %>% 
            mutate(subject = ID, time = time, conc = PKCONC, plt = PLT, circ0 = circ0_2) %>%
            distinct()
        
        
        # browser()
        PLT <- union_all(PLT1, PLT2) %>% as.data.frame()
        PLT
    })
    
    
    PLTSUM <- reactive({
        # browser()
        
        tmp <- PLT() %>% group_by(time, circ0) %>%
            mutate(
                p_upper = quantile(plt, probs = c(input$p_pop[2]), na.rm = TRUE),
                p_lower = quantile(plt, probs = c(input$p_pop[1]), na.rm = TRUE),
                p_50 = quantile(plt, probs = c(0.5), na.rm = TRUE)
            ) %>% 
            select(., subject, time, plt, p_upper, p_lower, p_50, circ0) %>% 
            ungroup() %>% 
            filter(subject == 1)
    })
    
    # output$nadirpop <- renderText({
    #     # browser()
    #     
    #     PLTSUM <- PLTSUM()
    #     nadirpop <- PLTSUM$p_lower %>% min()
    #     round(nadirpop, digits = 0)
    # })
    
    output$nadirpop1 <- renderText({
        # browser()
        PLTSUM <- PLTSUM()
        nadiravg <- filter(PLTSUM, circ0 == input$circ0_1_pltpop)$p_50 %>% min()
        round(nadiravg, digits = 0)
    })
    output$nadirpop2 <- renderText({
        # browser()
        PLTSUM <- PLTSUM()
        nadiravg <- filter(PLTSUM, circ0 == input$circ0_2_pltpop)$p_50 %>% min()
        round(nadiravg, digits = 0)
    })
    
    output$distPlotPLTpop <- renderPlot({
        PLT <- PLT()
        PLTSUM <- PLTSUM()
        
        circ0_1 <- input$circ0_1_pltpop
        circ0_2 <- input$circ0_2_pltpop
        
        # browser()
        
        ggplot(PLTSUM, aes(x = time/24 + 1, y = plt, group = subject)) + 
            # geom_line(alpha=0.2) +
            scale_x_continuous(breaks = seq(1, 28*12+1, 7)) +
            scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 150 , 200, 250, 300, 400, 500)) +
            # scale_x_continuous(limits = c(input$obs[1], input$obs[2]), breaks = seq(0,100*24,24)) +
            # scale_y_continuous(limits = c(0, 500)) +
            
            geom_ribbon(data=filter(PLTSUM
                                    , circ0 == circ0_1), 
                        aes(x=time/24 + 1, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.3, fill="#002288") +
            geom_line(data=filter(PLTSUM
                                  , circ0 == circ0_1), 
                      aes(x=time/24 + 1, y=p_50), linetype="solid", size = 1, colour="#002288") +
            
            geom_ribbon(data=filter(PLTSUM, circ0 == circ0_2), 
                        aes(x=time/24 + 1, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.2, fill="#880022") +
            geom_line(data=filter(PLTSUM, circ0 == circ0_2), 
                      aes(x=time/24 + 1, y=p_50), linetype="solid", size = 1, colour="#880022") +
            
            # geom_hline(yintercept = 75, linetype="dashed", size = 0.5, colour = "black", alpha = 0.8) +
            geom_hline(yintercept = 50, linetype="dashed", size = 0.5, colour = "black", alpha = 0.8) +
            geom_hline(yintercept = 25, linetype="dashed", size = 0.8, colour = "black", alpha = 0.8) +
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 21),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                  axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            
            xlab("Days after treatment") + ylab("Platelets (10^9)")
    })
    
    
    TUMOR <- reactive({
        
        cyc <- input$cycle_pop
        
        amt1 <- input$amt1_pop
        tau1 <- input$interval1_pop
        n1 <- input$n1_pop
        
        amt2 <- input$amt2_pop
        tau2 <- input$interval2_pop
        n2 <- input$n2_pop
        
        idv <- input$idv_pop
        
        
        
        dur <- cyc*(n1*tau1 + n2*tau2)
        dur_tgi <- 26*24
        amt <- cyc*(n1*amt1 + n2*amt2)
        amt_tgi <- cyc*amt/dur*dur_tgi
        
        
        amt1_tgi <- amt1/dur*dur_tgi
        tau1_tgi <- tau1/dur*dur_tgi
        # n1_tgi <- n1
        
        amt2_tgi <- amt2/dur*dur_tgi
        tau2_tgi <- tau2/dur*dur_tgi
        # n2_tgi <- n2
        
        cl <- input$cl_pkpop
        v2 <- input$v2_pkpop
        v3 <- input$v3_pkpop
        q <- input$q_pkpop
        ka <- input$ka_pkpop
        
        lamda0 <- input$lamda0_tumorpop
        emax <- input$emax_tumorpop
        ec50 <- input$ec50_tumorpop
        
        
        dos1 <- expand.grid(cmt = 1,
                            time = c(0,
                                     c(1:cyc) * (tau1_tgi*n1 + tau2_tgi*n2)),
                            amt = amt1_tgi,
                            ii = tau1_tgi,
                            addl = n1 - 1,
                            CL = cl,
                            V2 = v2,
                            V3 = v3,
                            Q = q,
                            KA = ka,
                            LAMDA0 = lamda0,
                            EMAX = emax,
                            EC50 = ec50,
                            ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt)
        
        dos2 <- expand.grid(cmt = 1,
                            time = c(c(tau1_tgi*n1 + tau2_tgi*n2, 
                                       c(1:cyc) * (tau1_tgi*n1 + tau2_tgi*n2)) - tau2_tgi*n2)[-1],
                            amt = amt2_tgi,
                            ii = tau2_tgi,
                            addl = n2 - 1,
                            CL = cl,
                            V2 = v2,
                            V3 = v3,
                            Q = q,
                            KA = ka,
                            LAMDA0 = lamda0,
                            EMAX = emax,
                            EC50 = ec50,
                            ID = seq(1, idv, 1)) %>% 
            mutate(dose = amt)
        
        dos0 <- expand.grid(cmt = 1,
                            time = 0,
                            amt = 0,
                            ii = 24,
                            addl = 100,
                            CL = cl,
                            V2 = v2,
                            V3 = v3,
                            Q = q,
                            KA = ka,
                            LAMDA0 = lamda0,
                            EMAX = emax,
                            EC50 = ec50,
                            ID = 9999) %>%
            mutate(dose = amt)
        
        dos <- union_all(dos1, dos2) %>% 
            union_all(dos0) %>%
            mutate(evid = 1) %>% 
            arrange(ID, time, addl, ii)
        
        out <- mod_tumorpop %>%
            data_set(dos) %>%
            mrgsim(set.seed(2020),
                   end = dur_tgi,
                   delta = 12)
        
        # browser()
        TUMOR <- out@data %>%
            select(ID, time, TUMOR) %>% 
            mutate(subject = ID, time = time, tumor = TUMOR) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        TUMOR
    })
    
    
    TUMORSUM <- reactive({
        # browser()
        tmp <- TUMOR() %>% filter(ID != 9999) %>% 
            group_by(time) %>% mutate(
                p_upper = quantile(tumor, probs = c(input$p_pop[2]), na.rm = TRUE),
                p_lower = quantile(tumor, probs = c(input$p_pop[1]), na.rm = TRUE),
                p_50 = quantile(tumor, probs = c(0.5), na.rm = TRUE)) %>% 
            select(., subject, time, tumor, p_upper, p_lower, p_50) %>% 
            ungroup() %>% 
            filter(subject == 1)
    })
    
    output$tgipop1 <- renderText({
        
        # browser()
        TUMOR <- TUMOR()
        TUMORSUM <- TUMORSUM()
        TUMOR0 <- filter(TUMOR, ID == 9999)
        
        
        t0 <- TUMOR0 %>% filter(ID == 9999) %>% .$tumor %>% max()
        t <- TUMORSUM %>% .$p_50 %>% max()
        
        # browser()
        round(abs(t - t0)/(t0-100) * 100, digits = 0)
    })
    
    output$tgipop2 <- renderText({
        
        # browser()
        TUMOR <- TUMOR()
        TUMORSUM <- TUMORSUM()
        TUMOR0 <- filter(TUMOR, ID == 9999)
        
        
        t0 <- TUMOR0 %>% filter(ID == 9999) %>% .$tumor %>% max()
        t <- TUMORSUM %>% .$p_50 %>% max()
        
        # browser()
        round(abs(t - t0)/(t0-100) * 100, digits = 0)
    })
    
    
    output$distPlotTumorpop <- renderPlot({
        
        TUMOR <- TUMOR()
        TUMORSUM <- TUMORSUM()
        TUMOR0 <- filter(TUMOR, ID == 9999)
        # browser()
        
        ggplot(TUMOR, aes(x = time, y = tumor, group = subject)) +
            
            scale_y_continuous(breaks = seq(0, 10000, 200)) +
            # geom_line(alpha=0.2) +
            
            geom_ribbon(data=TUMORSUM, 
                        aes(x=time, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.3, fill="#1D8348") +
            geom_line(data=TUMORSUM, 
                      aes(x=time, y=p_50), linetype="solid", size = 0.8, colour="black") +
            geom_line(data=TUMOR0, 
                      aes(x=time, y=tumor), linetype="solid", size = 0.8, colour="red") +
            
            # geom_hline(aes(yintercept = 100.24), colour = "grey", linetype = "dashed", size = 1, alpha = 0.8) +
            
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 21),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                  axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            xlab("Time (h)") + ylab("Tumor size (mmm^3)")
    })
    
    
    
    
    
    output$distPlotPKPDpop <- renderPlot({
        
        
        # PK
        PK <- PK()
        PKSUM <- PKSUM()
        # browser()
        
        plotpk <- ggplot(PK, aes(x = time, y = conc, group=subject)) +
            scale_x_continuous(breaks = seq(0, 168*4*12, 168)) +
            scale_y_continuous(breaks = c(0, 25, 50, 100, 200, 300, 400, 500, 600), limits = c(0, 350)) +
            # geom_line(alpha=0.2) +
            geom_ribbon(data=PKSUM, aes(x=time, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.3, fill="#566573") +
            geom_line(data=PKSUM, aes(x=time, y=p_50), linetype="solid", size = 0.8, colour="black") +
            
            geom_hline(aes(yintercept = 25), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            geom_hline(aes(yintercept = 50), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            geom_hline(aes(yintercept = 100), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            geom_hline(aes(yintercept = 200), colour = "#CB4335", linetype = "dashed", size = 0.5, alpha = 0.8) +
            
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 21),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                  axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            xlab("Time (h)") + ylab("Plasma concentration (ng/mL)")   
        
        
        # PK-PLT
        PLT <- PLT()
        PLTSUM <- PLTSUM()
        
        circ0_1 <- input$circ0_1_pltpop
        circ0_2 <- input$circ0_2_pltpop
        
        # browser()
        
        plotplt <- ggplot(PLT, aes(x = time/24 + 1, y = plt, group = subject)) + 
            # geom_line(alpha=0.2) +
            scale_x_continuous(breaks = seq(1, 28*12+1, 7)) +
            scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 150 , 200, 250, 300, 400, 500)) +
            # scale_x_continuous(limits = c(input$obs[1], input$obs[2]), breaks = seq(0,100*24,24)) +
            # scale_y_continuous(limits = c(0, 500)) +
            
            geom_ribbon(data=filter(PLTSUM
                                    , circ0 == circ0_1), 
                        aes(x=time/24 + 1, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.3, fill="#002288") +
            geom_line(data=filter(PLTSUM
                                  , circ0 == circ0_1), 
                      aes(x=time/24 + 1, y=p_50), linetype="solid", size = 1, colour="#002288") +
            
            geom_ribbon(data=filter(PLTSUM, circ0 == circ0_2), 
                        aes(x=time/24 + 1, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.2, fill="#880022") +
            geom_line(data=filter(PLTSUM, circ0 == circ0_2), 
                      aes(x=time/24 + 1, y=p_50), linetype="solid", size = 1, colour="#880022") +
            
            # geom_hline(yintercept = 75, linetype="dashed", size = 0.5, colour = "black", alpha = 0.8) +
            geom_hline(yintercept = 50, linetype="dashed", size = 0.5, colour = "black", alpha = 0.8) +
            geom_hline(yintercept = 25, linetype="dashed", size = 0.8, colour = "black", alpha = 0.8) +
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 21),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                  axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            
            xlab("Days after treatment") + ylab("Platelets (10^9)")
        
        # PK-TUMOR
        TUMOR <- TUMOR()
        TUMORSUM <- TUMORSUM()
        TUMOR0 <- filter(TUMOR, ID == 9999)
        # browser()
        
        plottumor <- ggplot(TUMOR, aes(x = time, y = tumor, group = subject)) +
            
            scale_y_continuous(breaks = seq(0, 10000, 200)) +
            # geom_line(alpha=0.2) +
            
            geom_ribbon(data=TUMORSUM, 
                        aes(x=time, ymin=p_lower, ymax=p_upper), linetype=0, alpha=0.3, fill="#1D8348") +
            geom_line(data=TUMORSUM, 
                      aes(x=time, y=p_50), linetype="solid", size = 0.8, colour="black") +
            geom_line(data=TUMOR0, 
                      aes(x=time, y=tumor), linetype="solid", size = 0.8, colour="red") +
            
            # geom_hline(aes(yintercept = 100.24), colour = "grey", linetype = "dashed", size = 1, alpha = 0.8) +
            
            theme_bw(base_rect_size = 1) +
            theme(axis.text = element_text(size = 20),
                  axis.title = element_text(size = 21),
                  axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                  axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            xlab("Time (h)") + ylab("Tumor size (mmm^3)")
        
        
        
        gridExtra::grid.arrange(plotpk, plotplt, plottumor)
        
    })
    
})


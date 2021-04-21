## R shiny app to Simulate Disease Seasonality Using Temporally Forced Models
## and Test how seasonality changes with R0 and amplitude of seasonality
## Model: SEIR; Parameters based on childhood infections
## 4/15/18, by Wan Yang

library(shiny)
library(deSolve)


HIVTB = function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = caplambda - (beta*c*S*(I+J3)/N) - lambda*sigma*S*(Jas/R) - mu*S
    dL = beta*c * (S + Tb)*(I+J3)/N - lambda*sigma*L*(Jas/R) - ( mu+k+r1)*L
    dI = k*L - (mu+d+r2)*I
    dTb = r1*L + r2*I - beta*c*Tb*((I+J3)/N) - lambda*sigma*Tb*(Jas/R) - mu*Tb
    dJ1 = lambda*sigma*(S+Tb)*(Jas/R) - beta*c*J1*(I+J3)/N - (alpha1 + mu)*J1 + ras*J2
    dJ2 = lambda*sigma*L*(Jas/R) + beta*c*J1*(I+J3)/N - (alpha1 + mu+kas+ras)*J2
    dJ3 = kas*J2 - (alpha3 + mu + das)*J3
    dA = alpha1*J1 + alpha2*J2 + alpha3*J3 - (mu+f)*A
    list(c(dS, dL, dI, dTb, dJ1, dJ2, dJ3, dA))
  })
}
N=1; # 100%

ui <- fluidPage(h4('Simulating TB-HIV Co-infections with changes in R0'),
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'R01',h6('beta for TB'),value=1.7,min=0,max=20,step=.1),
                               sliderInput(inputId = 'R02',h6('sigma for HIV '),value=2.5,min=0,max=20,step=.01),
                               sliderInput(inputId = 'N',h6('Population size'),value=6800000,min=5000000,max=20000000,step=100000),
                               #sliderInput(inputId = '',h6('R0 for HIV '),value=0,min=0,max=20,step=.01),
                               h6('beta: TB transmission rate'),
                               h6('lambda: HIV transmission rate')
                               
                  ),
                  
                  mainPanel(
                    plotOutput(outputId = 'plots',width = '100%', height = "550px")
                  )
                )
                
)

server <- function(input, output){
  
  output$plots=renderPlot({
    
    #R01=input$R01         # R0 for all three runs
    #R02=input$R02
    
    N = input$N
    
    # other parameters and initial conditions
    TN = 8336817; N = N * TN; L = 0.16 * N; I = 6.9/100000 * N;
    Tb = 0.82 * I; J1 = 1458/100000 *N; J2 = .14*J1; J3 = 0.6 * I;
    A = 0.166 * J1; mu = 555.1/100000; gamma1 = .95; 
    R = N - I - J3 - A;
    R01 = 1.7;
    beta = R01 * gamma1; gamma2 = 0.67;
    #lambda = R02*(gamma2 + mu);
    
    lambda = 4.025/100000
    
    
    state = c(
      S = N - L - J1 - J2 - J3 , L = L, I = I, Tb = Tb, J1 = J1,
      J2 = J2, J3 = J3, A = A
    )
    
    
    params = c(
      N = N,
      R = R,
      Jas = J1 + J2 + J3,
      caplambda = ((.175*TN)/15),
      beta = beta,
      lambda = lambda,
      c = (10 * I) / R,
      sigma = 0.009/100000 ,
      mu = mu ,
      k = 0.08,
      kas = 0.1,
      d = 0.9,
      das = 4.7/1000 ,
      f = 21.35/100000 ,
      r1 = 0.78 ,
      r2 = 0.82 ,
      ras = 1,
      alpha1 = 10/100000 ,
      alpha2 = 0.25 ,
      alpha3 = 0.25
    )
    
    #tm_step=1; # run it by week
    times=seq(0,3650,by=1) # need to run for ~10 yrs to have the 3rd sim stable!
    
    sim=ode(y=state,times=times,func=HIVTB,parms=params)
    
    
    # plot results
    par(mfrow=c(3,1),cex=1,mar=c(3,3,1,1),oma=c(0,0,1,0),mgp=c(1.8,.5,0))
    
    plot(sim[,'time'],(sim[,'I']),ylab='TB infectious',xlab='Time (year)',
         type='l',lwd=2, ylim=c(-5,20000))
    mtext(bquote(lambda[1]==.(lambda)),side=3,line=-2.5,outer=F,adj=.95,cex=1.2)
    mtext(bquote(beta[1]==.(beta)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    #mtext(bquote(R[0]==.(R01)),side=3,line=0,outer=F,adj=.5,cex=1.5)
    #abline(v=990:1000,col='grey',lty=2)
    
    plot(sim[,'time'],(sim[,'J1']),ylab='HIV infectious',xlab='Time (year)',
         type='l',lwd=2, ylim=c(0,1500000))
    mtext(bquote(lambda[1]==.(lambda)),side=3,line=-2.5,outer=F,adj=.95,cex=1.2)
    mtext(bquote(beta[1]==.(beta)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
   
    #  abline(v=990:1000,col='grey',lty=2)
    # plot(sim[,'time'],(sim[,'J3']),ylab='TB-HIV infectious',xlab='Time (year)',
    #      type='l',lwd=2, ylim=c(-5,2000))
    # mtext(bquote(lambda[1]==.(lambda)),side=3,line=-2.5,outer=F,adj=.95,cex=1.2)
    # mtext(bquote(beta[1]==.(beta)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    # abline(v=990:1000,col='grey',lty=2)
   
     matplot(simR[,6:8], type = 'l', lty = 6:8, col = c('red', 'green', 'blue'), xlab = 'Time', ylab = 'TB')
    legend('topright', c('J1', 'J2', 'J3'), lty = 6:8, col = c('red', 'green', 'blue'))
    
  })
  
}

shinyApp(ui=ui, server = server)


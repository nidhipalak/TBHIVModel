## R shiny app to Simulate Disease Seasonality Using Temporally Forced Models
## and Test how seasonality changes with R0 and amplitude of seasonality
## Model: SEIR; Parameters based on childhood infections
## 4/15/18, by Wan Yang

library(shiny)
library(deSolve)


HIVTB = function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = caplambda - (betac*S*(I+J3)/N) - lambdasigma*S*(Jas/R) - mu*S
    dL = betac * (S + Tb)*(I+J3)/N - lambdasigma*L*(Jas/R) - ( mu+k+r1)*L
    dI = k*L - (mu+d+r2)*I
    dTb = r1*L + r2*I - betac*Tb*((I+J3)/N) - lambdasigma*Tb*(Jas/R) - mu*Tb
    dJ1 = lambdasigma*(S+Tb)*(Jas/R) - betac*J1*(I+J3)/N - (alpha1 + mu)*J1 + ras*J2
    dJ2 = lambdasigma*L*(Jas/R) + betac*J1*(I+J3)/N - (alpha1 + mu+kas+ras)*J2
    dJ3 = kas*J2 - (alpha3 + mu + das)*J3
    dA = alpha1*J1 + alpha2*J2 + alpha3*J3 - (mu+f)*A
    list(c(dS, dL, dI, dTb, dJ1, dJ2, dJ3, dA))
  })
}
N=1; # 100%

ui <- fluidPage(h4('Simulating TB-HIV Co-infections with changes in transmission-contact rate'),
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'betac',h6('transmission-contact rate TB'),value=10,min=0,max=20,step=1),
                               sliderInput(inputId = 'lambdasigma',h6('transmission-contact rate HIV '),value=10,min=0,max=20,step=1)
                               
                  ),
                  
                  mainPanel(
                    plotOutput(outputId = 'plots',width = '100%', height = "550px")
                  )
                )
                
)

server <- function(input, output){
  
  output$plots=renderPlot({
    
    betac=input$betac         # R0 for all three runs
    lambdasigma=input$lambdasigma
    
    
    
    # other parameters and initial conditions
    TN = 8336817; N = .825 * TN; L = 0.16 * N; I = 6.9/100000 * N;
    Tb = 0.82 * I; J1 = 1458/100000 *N; J2 = .14*J1; J3 = 0.6 * I;
    A = 0.166 * J1; mu = 555.1/100000; gamma1 = .95; 
    R = N - I - J3 - A;
    #R01 = 1.7;
    #beta = R01 * gamma1; gamma2 = 0.67;
    #lambda = R02*(gamma2 + mu);lambda = 4.025/100000
    
    
    state = c(
      S = N - L - J1 - J2 - J3 , L = L, I = I, Tb = Tb, J1 = J1,
      J2 = J2, J3 = J3, A = A
    )
    
    
    params = c(
      N = N,
      R = R,
      Jas = J1 + J2 + J3,
      caplambda = ((.175*TN)/15)/N,
      betac = betac,
      lambdasigma = lambdasigma,
      #c = (10 * I) / R,
      #sigma = 0.009/100000 ,
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
    times=seq(0,500,by=1) # need to run for ~10 yrs to have the 3rd sim stable!
    
    sim=ode(y=state,times=times,func=HIVTB,parms=params)
    
    # plot tb infectious
    par(mfrow=c(3,1),cex=1,mar=c(3,3,1,1),oma=c(0,0,1,0),mgp=c(1.8,.5,0))
    matplot(sim[,6:8]/N, type = 'l', lty = 1, lwd = c(2,2,3), 
            col = c('red', 'green', 'blue'), 
            xlab = 'Time (years)', ylab = 'proportion infectious',
            main=paste('Proportion infectious HIV'),cex.main=1)
    lines(((sim[,'I'] + sim[,'J3'])/N), lty= 2, col = "purple", type='l',lwd=2)
    lines(sim[,'I']/N, lty= 1, col = "black", type='l',lwd=2)
    #lines(sim[,'time'], (sim[,'I']/N, lty= 2, col = "black", type='l',lwd=2)
     legend('topright', c('infectious HIV', 'infectious HIV, latent TB', 
                          'infectious HIV, infectious TB', 'infectious TB + infectious HIV-TB',
                          'infectious TB'), lty = c(1,1,1,2,1), lwd = 2,
            col = c('red', 'green', 'blue', 'purple', 'black'))    
    
    # plot tb infectious proportion
    plot(sim[,'time'],sim[,"I"]/N, 
                       ylab= 'proportion infectious TB', xlab='Time (year)',
         type='l',lwd=2, col = c('black'),
         main=paste('Proportion infectious'),cex.main=1)
    lines((sim[,'I'] + sim[,'J3'])/N, type = 'l', lwd = 2, col = c('purple'))
    legend('topright', c("infectious TB", 'infectious TB + infectious HIV-TB'), 
           lwd=2, col = c('black', 'purple'))
    
     # plot AIDS
     plot(sim[,'time'],(sim[,'A']/N),ylab='AIDS',xlab='Time (year)',
          type='l',lwd=2, col="orange",
          main=paste('Proportion with AIDS'),cex.main=1)
     legend('topright', c('AIDS'), 
            lwd=2, col = c('orange'))

  })
  
}

shinyApp(ui=ui, server = server)


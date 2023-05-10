library(shiny)
library(ggplot2)
# themes
########### The app might not call these right; check what folder they are in.
# mine had problems when I moved these files into a package.
# source("RLS_predicted_performance.R")
# source("InfluenceFunctions.R")
# source("Time-Invariant_Functions.R")

library(devtools)
load_all()

#need to do
#   optionally plot performance
#   automatically create day column
#   add 0 to day column if there is no zero day



######!!!!!comments that look like this are things that might need to be fixed.


#### cute functional programing stuff
params_tmp <- c("p_0","k_1","k_2","tau_1","tau_2")

inputPar <- function(id){
  textInput(id, label = paste("Choose",id,sep=" "), value = NULL, placeholder ="initial guess")
}

inputParamsInvariant <- purrr::map(params_tmp,inputPar)


vars <- tibble::tribble(
  ~ id,   ~ label,
  "days",   "Select Day Variable",
  "trainingLoad", "Select Training Load",
  "performance",    "Select Performance"
)

myInput <- function(id,label){
  selectInput(id,label,choices=NULL)
}
inputsData <- purrr::pmap(vars, myInput)


#### Same thing for Time-Varying
num_input_variant <- c("p_0", "alpha", "by_2", "by_1")

inputsTimeVarying <- c(purrr::map(num_input_variant,inputPar),NA,NA)


#### make by_1 and by_2 interactive with the slider ####
inputsTimeVarying[[6]] = inputsTimeVarying[[3]]
inputsTimeVarying[[3]] = sliderInput("bounds_1", "Choose Bounds_1",
                                     1, 100, c(1,50), 1)
inputsTimeVarying[[5]] = sliderInput("bounds_2", "Choose Bounds_1",
                                     1, 100, c(1,50), 1)



#### For switching between Time-Inv and Time-Var
inputParams <- tabsetPanel(
  id="model_type",
  type="hidden",
  tabPanel("Time-Invariant",
           inputParamsInvariant),
  tabPanel("Time-Variant",
           inputsTimeVarying)
)


#### appUI
ui=fluidPage(
  theme = bslib::bs_theme(
    bootswatch = "yeti"
  ),
   # alxdfj;als r
  tabsetPanel(
    id="wizard",
    #can be hidden
    # type="hidden",

    #choosing dataset...
    tabPanel("page_1",
             fluidRow(
               fileInput("file", "Choose CSV file:",
                         multiple=TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),

               #...which lets you go to the next page when one is selected
               uiOutput("nextPage")
             ),

             fluidRow(
               tableOutput("head")
             )
    ),

    #Choosing inputs and graphing
    tabPanel("page_2",
             fluidRow(
               column(3,
                      actionButton("page_21", "Change Data"),
                      inputsData
               ),

               column(2,
                      actionButton("optim", "Optimize parameters"),
                      radioButtons("model_control", label = "Type of Model",
                                   choices=list("Time-Invariant" = "Time-Invariant",
                                                "Time-Variant" = "Time-Variant"),
                                   selected = "Time-Invariant"
                      ),
                      inputParams
               ),

               column(7,
                      actionButton("graph", "Graph!!"),
                      plotOutput("plotMain"),
                      plotOutput("plotTrainingLoad"),
                      plotOutput("plotInfluence"),
                      verbatimTextOutput("out")
               )
             )
    ),

    #to see what's up
    #tableOutput("Preview"),

  )
)


#appServer
server=function(input,output,session){

  #### page wizard
  output$nextPage <- renderUI({
    req(input$file)
    actionButton("page_12", "Change Inputs")
  })

  switch_page <- function(i) {
    updateTabsetPanel(inputId = "wizard", selected = paste0("page_", i))
  }

  observeEvent(input$page_12, switch_page(2))
  observeEvent(input$page_21, switch_page(1))


  #### reading dataset
  df <- eventReactive(input$file,{
    req(input$file)
    out=read.csv(input$file$datapath)
    out
  })

  training_load_on_graph <- eventReactive(input$graph,{
    req(input$trainingLoad)
    df()[[input$trainingLoad]]
  })

  performance_on_graph <- eventReactive(input$graph,{
    req(input$performance)
    df()[[input$performance]]
  })

  #### interactive choices for the columns of the dataset
  name=reactive(names(df()))

  observe({
    updateSelectInput(session,"performance", "Select Performance", choices=name())
    updateSelectInput(session,"trainingLoad","Select Training Load",choices=name())
    updateSelectInput(session,"days", "Select Day Variable", choices=name())
  })


  #### outputs the head of the chosen dataset
  output$head=renderTable(head(df(), c(6,4)))


  #### getting the parameter inputs from the user
  input_Params <- eventReactive(input$graph, {
    tmp=c(input$p_0,input$k_1,input$k_2,
          input$tau_1,input$tau_2)
    tmp=purrr::map_dbl(tmp,as.double)
    tmp
  })


  #### Time-Invariant predicted performance
  time_inv_pred_perf <- reactive({
    invariant_perf(input_Params(),training_load_on_graph())
  })


  #### days variable
  days <- eventReactive(input$graph, {
    df()[[input$days]]
  })


  #for testing
  #output$out=renderPrint(time_inv_pred_perf())


  #### allows the user to optimize the inputed parameters with R's
  #### built-in optim function. Returns the optimized values into the
  #### text boxes.
  updateInputPar <- function(id,value){
    updateTextInput(session=session, inputId=id,
                    label=paste("Choose",id,sep=" "),
                    value=value,
                    placeholder="number"
    )
  }

  optimParams <- observeEvent(input$optim,{
    tmp=optim_par(input_Params(),df()[[input$trainingLoad]],
                  df()[[input$performance]])
    purrr::map2(params_tmp,tmp,updateInputPar)
  })


  #### Time-Varying predicted Performance
  time_var_perf <- eventReactive(input$graph, {
    tmp = RLS_predicted_performance(training_load_on_graph(),
                                    performance_on_graph(),
                                    as.double(input$p_0),
                                    as.double(input$alpha),
                                    delta = 1000,
                                    as.double(input$bounds_1),
                                    as.double(input$by_1),
                                    as.double(input$bounds_2),
                                    as.double(input$by_2)
    )
    unlist(tmp[[1]])

  })


  #### Switching between Time-Inv and Time-Var
  observeEvent(input$model_control, {
    updateTabsetPanel(inputId = "model_type", selected = input$model_control)
  })


  #### Main Plot
  plotMain <- eventReactive(input$graph, {
    ggplot(df(),aes(x=days(),y=time_inv_pred_perf()))+
      geom_line(aes(y=time_inv_pred_perf(),color="black"), size=1)+
      geom_line(aes(y=time_var_perf(),color="blue"), size = 1)+
      geom_point(aes(y = df()[[input$performance]], color="red"), shape = 1)+
      scale_color_manual("", values = c("black","blue", "red"),
                         labels = c("Time-Invariant",
                                    "Time_Variant",
                                    "Actual Performance"))+
      labs(x="Days", y="Performance")
  })

  output$plotMain <- renderPlot({
    plotMain()
  })

  #### Influence plot
  ####!!!! day variable was acting wierd !!!!####
  df_influ <- reactive({
    out=data.frame(rev(-1*days()),
                   Influence(input_Params(),-max(days()), min(days()))
    )
  })

  plotInflu <- eventReactive(input$graph, {
    ggplot(df_influ(),aes(df_influ()[[1]],df_influ()[[2]]))+
      geom_line()+
      labs(x="Days", y="Influence")
  })

  output$plotInfluence <- renderPlot({plotInflu()})



  #### Training load plot
  plotTrain <- eventReactive(input$graph,{
    ggplot(df(),aes(days(), df()[[input$trainingLoad]]))+
      geom_bar(stat = "identity")+
      labs(x = "Days", y = "Training Load")
  })

  output$out <- renderPrint({time_inv_pred_perf})

  output$plotTrainingLoad <- renderPlot(plotTrain())

}


shinyApp(ui,server)




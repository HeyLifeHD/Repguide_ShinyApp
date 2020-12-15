#Joschka Hey
#Shiny web application to run the R library Repguide

# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

#libraries
if (!"shinyWidgets" %in%  installed.packages())
  install.packages("shinyWidgets")
library(shinyWidgets)
if (!"shinyBS" %in%  installed.packages())
  install.packages("shinyBS")
library(shinyBS)
if (!"shiny" %in%  installed.packages())
  install.packages("shiny")
library(shiny)
if (!"Repguide" %in%  installed.packages()) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  remotes::install_github(repo = 'tanaylab/repguide')
}
library(Repguide)
if (!"shinycssloaders" %in%  installed.packages()) {
  if (!requireNamespace("remotes", quietFly = TRUE))
    install.packages("remotes")
  remotes::install_github(repo = "andrewsali/shinycssloaders")
}
library(shinycssloaders)
if (!"data.table" %in%  installed.packages())
  install.packages("data.table")
library(data.table)
if (!"devtools" %in%  installed.packages())
  install.packages("devtools")
library(devtools)
#if (!"dplyr" %in%  installed.packages())
#  install_version("dplyr", version = "0.8.0.1", repos = "http://cran.us.r-project.org")
library(dplyr)
if (!"shinythemes" %in%  installed.packages())
  install.packages("shinythemes")
library(shinythemes)
if (!"GenomicFeatures" %in%  installed.packages())
  install.packages("GenomicFeatures")
library(GenomicFeatures)
# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("flatly"),
  navbarPage(
    strong(
      "Repguide |  Design of guideRNAs for CRISPR/dCas9 targeting of repetitive DNA sequences"
    ),
    id = "page",
    windowTitle = "Repguide",
    tabsetPanel(
      type = "tabs",
      tabPanel(
        "Home",
        titlePanel(title = div(
          img(
            src = "logo.png",
            height = 120,
            width = 120
          ),
          "About Repguide",
          align = "left"
        )),
        mainPanel(
          fluidRow(column(7, includeMarkdown("./home.Rmd")),
                   column(1),
                   column(3, includeMarkdown(
                     "./contributions.RMD"
                   ))),
          
          offset = 2,
          br(),
          #HTML("Repguide Visitor Counter:"),
          #br(),
          #HTML(
          #  "<script type='text/javascript' src='//counter.websiteout.net/js/7/8/0/0'></script>"
          #),
          #br(),
          width = 11
        )
      ),
      
      tabPanel(
        "Repeat exploration",
        titlePanel(title = div(
          img(
            src = "logo.png",
            height = 100,
            width = 100
          ),
          "Repguide - Repeat exploration",
          align = "left"
        )),
        sidebarLayout(
          sidebarPanel(
            align = "left",
            selectInput(
              "repeats_genome",
              label = h5("Select reference genome:"),
              choices = list(
                "Human GRCh38/hg38" = "hg38",
                "Human GRCh37/hg19" = "hg19",
                "Mouse GRCh38/mm10" = "mm10"
              ),
              selected = 2
            )
          ),
          #main panel
          mainPanel(dataTableOutput("repeats_table"))
        )
      ),
      
      tabPanel(
        "Guide design",
        titlePanel(title = div(
          img(
            src = "logo.png",
            height = 100,
            width = 100
          ),
          "Repguide - Workflow",
          align = "left"
        )),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "1.Target exploration",
            # Sidebar input selection
            sidebarLayout(
              sidebarPanel(
                align = "left",
                #h4("Specify Repguide options for Target exploration"),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  selectInput(
                    "ref_genome",
                    label = h5("Select reference genome:"),
                    choices = list(
                      "Human GRCh38/hg38" = "hg38",
                      "Human GRCh37/hg19" = "hg19",
                      "Mouse GRCh38/mm10" = "mm10"
                    ),
                    selected = 2
                  ),
                  circleButton(
                    "ques_ref_genome",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  textInput(
                    "target_repeats",
                    label = h5("Select target repeats:"),
                    value = NULL,
                    placeholder = "Enter repeat class..."
                  ),
                  circleButton(
                    "ques_target_repeats",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  textInput(
                    "whitelist_repeats",
                    label = h5("Select whitelist repeats:"),
                    value = NULL,
                    placeholder = "Enter repeat class..."
                  ),
                  circleButton(
                    "ques_whitelist_repeats",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  h5("Blacklist promoter region of essential genes:"),
                  circleButton(
                    "ques_whitelist_blacklist",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                checkboxInput("blacklist",
                              label = h5(""),
                              value = FALSE),
                
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "gap_frequencies",
                    label = h5("Maximum gap frequency for target alignment:"),
                    min = 0,
                    max = 1,
                    value = 0.8
                  ),
                  circleButton(
                    "ques_gap_frequencies",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                
                actionButton("action_target", "Submit target specifications")
              ),
              #main panel
              mainPanel(
                #plot first output: target repeats
                #h3("Target exploration:"),
                plotOutput("targets_repeats", height = "600px"),
                br(),
                conditionalPanel(
                  condition = "input.action_target ==true",
                  includeMarkdown("./figure_legends/target_repeats.RMD")
                )
              )
            ),
          ),
          tabPanel(
            "2.Guide generation",
            sidebarLayout(
              sidebarPanel(
                align = "left",
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  h5("Restrict guideRNA design to parts on consensus sequence"),
                  circleButton(
                    "ques_start_position",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                textInput(
                  "start_position",
                  label = h5("Relative start position:"),
                  value = "",
                  placeholder = "Enter start position, e.g. 50"
                ),
                textInput(
                  "end_position",
                  label = h5("Relative end position:"),
                  value = NULL,
                  placeholder = "Enter end position, e.g. 600"
                ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "guide_length",
                    label = h5("Length of guideRNAs:"),
                    min = 12,
                    max = 26,
                    value = 16
                  ),
                  circleButton(
                    "ques_guide_length",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "n_mismatches",
                    label = h5("Number of tolerated mismatches:"),
                    min = 0,
                    max = 3,
                    value = 0
                  ),
                  circleButton(
                    "ques_n_mismatches",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  textInput(
                    "five_prime_seq",
                    label = h5("Sequence requirement for 5' end:"),
                    value = NULL,
                    placeholder = "Enter nucleotide, e.g. G for transcription from U6 promoter"
                  ),
                  circleButton(
                    "ques_five_prime_seq",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "gc_content",
                    label = h5(" GC content:"),
                    min = 0,
                    max = 1,
                    value = c(0.4, 0.8)
                  ),
                  circleButton(
                    "ques_gc_content",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "min_Son",
                    label = h5("Minimal on target score of guides:"),
                    min = 0,
                    max = 100,
                    value = 50
                  ),
                  circleButton(
                    "ques_min_Son",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "alpha",
                    label = h5("Off-target score coefficient:"),
                    min = 0,
                    max = 100,
                    value = 0
                  ),
                  circleButton(
                    "ques_alpha",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                actionButton("action_guide", "Submit guide specifications")
              ),
              #main panel
              mainPanel(
                #plot first output: target repeats
                #h3("guideRNA design:"),
                plotOutput("guides", height = "800px"),
                br(),
                conditionalPanel(condition = "input.action_guide ==true",
                                 includeMarkdown("./figure_legends/guides.RMD"))
              )
            ),
          ),
          tabPanel(
            "3.Guide combination",
            
            
            sidebarLayout(
              sidebarPanel(
                align = "left",
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "max_guides",
                    label = h5("Maximum number of guides:"),
                    min = 0,
                    max = 30,
                    value = 5
                  ),
                  circleButton(
                    "ques_max_guides",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                # splitLayout(
                #   cellWidths = c("95%", "5%"),
                #   sliderInput(
                #     "iterations_gs",
                #     label = h5("Number of greedy search iterations:"),
                #     min = 0,
                #     max = 50,
                #     value = 10
                #   ),
                #   circleButton(
                #     "ques_iterations",
                #     icon = icon("question-circle"),
                #     size = "xs"
                #   )
                # ),
                splitLayout(
                  cellWidths = c("95%", "5%"),
                  sliderInput(
                    "alpha_combinations",
                    label = h5("Off-target score coefficient:"),
                    min = 0,
                    max = 100,
                    value = 10
                  ),
                  circleButton(
                    "ques_alpha_combinations",
                    icon = icon("question-circle"),
                    size = "xs"
                  )
                ),
                actionButton(
                  "action_combination",
                  "Submit guide combination specifications"
                ),
                br(),
                br(),
                downloadButton("downloadData", "Download Repguide output ")
              ),
              
              
              
              #main panel
              mainPanel(
                #plot first output: target repeats
                #h3("Combinatorial optimization:"),
                plotOutput("combinations", height = "1000px"),
                br(),
                conditionalPanel(
                  condition = "input.action_combination ==true",
                  includeMarkdown("./figure_legends/combinations.RMD")
                )
              )
            )
          )
        )
      ),
      tabPanel(
        "Terms and conditions",
        #main panel
        mainPanel(includeMarkdown("./terms_conditions.RMD")))
      
    )
  )
)




# Define server logic required to draw a histogram
server <- function(input, output) {
  #observer events for target exploration
  observeEvent(input$ques_ref_genome, {
    showModal(
      modalDialog(
        title = "Help",
        "Select reference organism genome for selection of target repeats.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_target_repeats, {
    showModal(
      modalDialog(
        title = "Help",
        "Choose your target repeat family you want to create guides for. You can retrieve the full list of available repeat families in the tab \"Repeat exploratin\".
      For example, LTRs originating from the human endogenous retrovirus 9 all conveniently share a common ‘LTR12’ prefix.fwhit",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_target_repeats, {
    showModal(
      modalDialog(
        title = "Help",
        "Choose your target repeat family you want to create guides for. You can retrieve the full list of available repeat families in the tab \"Repeat exploratin\".
      For example, LTRs originating from the human endogenous retrovirus 9 all conveniently share a common ‘LTR12’ prefix.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_whitelist_repeats, {
    showModal(
      modalDialog(
        title = "Help",
        "Choose one or several repeat families you want to whitlist from your analysis. Seperate multiple repeat families by \",\" Guides binding in whitelisted repeats will not account as off-targets.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_blacklist, {
    showModal(
      modalDialog(
        title = "Help",
        "If you select \"blacklist\", promoters of essential genes will be blacklisted.
      By this guides binding in these regions will be removed from the set. Currently only available for hg19 and hg38",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_gap_frequencies, {
    showModal(
      modalDialog(
        title = "Help",
        "Select the maximum gap frequency of the repeat family multiple sequence alignment. Positions on alignment with a higher gap will be removed.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  #repeats
  repeats_path_invest <- reactive({
    file.path(
      getwd(),
      "reference_data",
      "repeats",
      paste0("repeats_", input$repeats_genome, ".txt.gz")
    )
  })
  repeats_invest <- reactive({
    withProgress(message = "Load repeat data...", {
      if (!"DT" %in%  installed.packages())
        install.packages("DT")
      library(DT)
      req(repeats_path_invest)
      temp_repeats <- fread(repeats_path_invest())
      temp_repeats <- temp_repeats[, c(
        "genoName",
        "genoStart",
        "genoEnd",
        "strand",
        "repName",
        "repClass",
        "repFamily"
      )]
      
      colnames(temp_repeats) <- c(
        "chromosome",
        "start",
        "end",
        "strand",
        "Repeat Name",
        "Repeat Class",
        "Repeat Family"
      )
      temp_repeats1 <- as.data.frame(temp_repeats)
      return(temp_repeats1)
    })
  })
  
  
  output$repeats_table <- renderDataTable(repeats_invest(),
                                          options = list(pageLength = 30))
  #   output$repeats_table <- renderDataTable(
  #     req(repeats_inspect()),
  #     withProgress(message = "Preparae repeats table...", {
  #     as.data.frame(repeats_inspect()),
  #                                           options = list(pageLength = 5))
  #
  # })
  #
  
  #select organism
  org <- reactive({
    withProgress(message = "Load reference genome...", {
      if (input$ref_genome %in% c("hg38", "hg19")) {
        org <- "Hs"
      } else if (input$ref_genome %in% c("mm10", "mm9")) {
        org <- "Mm"
      } else{
        print("No appropriate reference genome chosen")
      }
      #load packages
      if (!paste0("org.", org, ".eg.db") %in%  installed.packages())
        install.packages(paste0("org.", org, ".eg.db"))
      library(paste0("org.", org, ".eg.db"), character.only = TRUE)
      org <- get(paste0("org.", org, ".eg.db"))
      return(org)
    })
  })
  
  #get appropriate reference genome
  BS_genome <- reactive({
    withProgress(message = "Load BS genome...", {
      if (!"BSgenome" %in%  installed.packages())
        install.packages("BSgenome")
      if (input$ref_genome %in% c("hg38")) {
        if (!"BSgenome.Hsapiens.UCSC.hg38" %in%  installed.packages())
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
      } else if (input$ref_genome %in% c("hg19")) {
        if (!"BSgenome.Hsapiens.UCSC.hg19" %in%  installed.packages())
          BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
      } else if (input$ref_genome %in% c("mm10")) {
        if (!"BSgenome.Mmusculus.UCSC.mm10" %in%  installed.packages())
          BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
      } else{
        print("No appropriate BSgenome object found")
      }
      BSgenome::getBSgenome(genome = input$ref_genome)
    })
  })
  
  #get appropriate txdb file
  txdb <- reactive({
    withProgress(message = "Load TXDB...", {
      txdb_objects <- grep(
        pattern = "^TxDb",
        x = rownames(installed.packages()),
        value = TRUE
      )
      txdb_objects <-
        as.data.frame(sapply(txdb_objects, function(x) {
          x <- unlist(strsplit(x, ".", fixed = TRUE))
          x
        }))
      txdb <-
        colnames(txdb_objects[, txdb_objects[4,] == input$ref_genome, drop = FALSE])
      if (!txdb %in%  installed.packages())
        BiocManager::install(txdb)
      library(txdb, character.only = TRUE)
      txdb <- get(txdb)
      return(txdb)
    })
  })
  
  #load data
  #fantom
  fantom_path <- reactive({
    file.path(file.path(
      getwd(),
      "reference_data",
      "fantom",
      paste0("fantom_", input$ref_genome, ".bed.gz")
    ))
  })
  #repeats
  repeats_path <- reactive({
    file.path(
      getwd(),
      "reference_data",
      "repeats",
      paste0("repeats_", input$ref_genome, ".txt.gz")
    )
  })
  
  
  
  #select whitelist repeats
  whitelist_regions <- reactive({
    withProgress(message = "Load Whitelist regions...", {
      whitelist_repeats <-
        unlist(ifelse(
          input$whitelist_repeats == "",
          list(NULL),
          input$whitelist_repeats
        ))
      whitelist_repeats <- if (is.null(whitelist_repeats)) {
        return(NULL)
      } else {
        repeats_dt <-
          fread(repeats_path())
        repeats_dt <-
          makeGRangesFromDataFrame(
            as.data.frame(repeats_dt),
            seqnames.field = "genoName",
            start.field = "genoStart",
            end.field = "genoEnd",
            keep.extra.columns = TRUE
          )
        whitelist_list <-
          as.vector(unlist(strsplit(whitelist_repeats, ",", fixed = TRUE)))
        whitelist_list <- gsub(" ", "", whitelist_list)
        if (is.element(whitelist_list, repeats_dt$repName)) {
          return(repeats_dt[repeats_dt$repName %in% whitelist_list,])
        } else {
          showNotification(
            paste0(
              whitelist_repeats,
              " not found in guideSet annotation."
            ),
            type = "warning",
            duration = 10
          )
          return(NULL)
        }
      }
      return(whitelist_repeats)
    })
  })
  
  #select blacklist regions
  blacklist_regions <- reactive({
    withProgress(message = "Load blacklist regions...", {
      if (isFALSE(input$blacklist)) {
        return(NULL)
      } else{
        GENE <-  genes(txdb())
        essentials <- fread(
          file.path(
            getwd(),
            "reference_data",
            "blacklist",
            "core-essential-genes-sym_HGNCID.txt"
          )
        )
        Promoter <-
          promoters(GENE, upstream = 1500, downstream = 1500)
        promoter_ids <-
          mapIds(
            org(),
            keys = essentials$AAMP,
            keytype = "SYMBOL",
            column = "ENTREZID",
            multiVals = "first"
          )
        # select(
        #     org(),
        #     keys = essentials$AAMP,
        #     columns = c("ENTREZID", "SYMBOL"),
        #     keytype = "SYMBOL"
        # )
        essentials_promoter <-
          Promoter[Promoter$gene_id %in% promoter_ids,]
        return(essentials_promoter)
      }
    })
  })
  
  #create guideset
  gs <- eventReactive(input$action_target, {
    withProgress(message = "Create guide set...", {
      message("Check: start createGuideSet")
      gs_temp <- createGuideSet(
        genome = BS_genome(),
        alt_chromosomes = FALSE,
        #input$alt_chromosomes,
        tes = repeats_path(),
        temp = tempdir(),
        # Directory for temporary files
        cis = fantom_path(),
        blacklist = blacklist_regions(),
        #input$blacklist,
        whitelist = whitelist_regions(),
        n_cores = 12,
        refdir = "",
        seed = 19
      )
      return(gs_temp)
      message("Check: end createGuideSet")
    })
    
  })
  gs1 <- reactive({
    req(gs())
    withProgress(message = "Add targets...", {
      message("Check: start addTargets")
      gs1_temp <- addTargets(gs(),
                             targets = input$target_repeats,
                             force = TRUE)
      return(gs1_temp)
      message("Check: end addTargets")
    })
  })
  
  gs2 <- reactive({
    req(gs1())
    withProgress(message = "Add alignments...", {
      message("Check: start addAlignments")
      gs2_temp <- addAlignments(
        guideSet = gs1(),
        max_gap_freq = input$gap_frequencies,
        iterations = 2,
        refinements = 1,
        force = TRUE
      )
      return(gs2_temp)
      message("Check: end addAlignments")
    })
    
  })
  
  output$targets_repeats <- renderPlot({
    req(gs2())
    #plot target repeat class
    withProgress(message = "Plot target sequences...", {
      plotTargets(gs2())
    })
  })
  
  #observer events for create guides
  observeEvent(input$ques_start_position, {
    showModal(
      modalDialog(
        title = "Help",
        "Select the start and end position on the consensus model of your target repeats. Guides binding outside the consensus range will be scored neutrally.
      This allows you to select guides that have a preferential binding in certain areas of the repeat, e.g. upstream of potential promoters, indicated by enrichment
      of Cis regulatory elements.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_n_mismatches, {
    showModal(
      modalDialog(
        title = "Help",
        "Maximal number of tolerated mismatches when assessing guideRNA binding targets.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_min_Son, {
    showModal(
      modalDialog(
        title = "Help",
        "Minimal on target score of guides.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_guide_length, {
    showModal(
      modalDialog(
        title = "Help",
        "Basepair size of the guideRNAs. guideRNA lenght between 12 and 26 are allowed",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_gc_content, {
    showModal(
      modalDialog(
        title = "Help",
        "Select range of GC content for guideRNAs. Will blacklist guides with GC content outside the range.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_alpha, {
    showModal(
      modalDialog(
        title = "Help",
        "Select the start and end position on the consensus model of your target repeats. Guides binding outside the consensus range will be scored neutrally.
      This allows you to select guides that have a preferential binding in certain areas of the repeat, e.g. upstream of potential promoters, indicated by enrichment
      of Cis regulatory elements.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_five_prime_seq, {
    showModal(
      modalDialog(
        title = "Help",
        "Sequence requirement for 5' start of guideRNAs, e.g. G nucleotide for transcription from U6 promoter.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  #get consensus sequence for guide design
  consensus_df <- eventReactive(input$action_guide, {
    message("Check: start consensus")
    consensus_df_temp <- data.frame(
      repname = c(input$target_repeats),
      start = c(input$start_position),
      end = c(input$end_position)
    ) #untill tss
    return(consensus_df_temp)
    message("Check: end consensus")
    
  })
  
  #Create Guides
  gs3 <-  reactive({
    req(gs2())
    message("Check: start addGuides")
    withProgress(message = "Add guides...", {
      gs3_temp <- addGuides(
        gs2(),
        # our guideSet
        n_mismatches = input$n_mismatches,
        # max allowed mismatch of reported binding sites
        min_Son = input$min_Son,
        # minimal score requirement for valid guides
        guide_length = input$guide_length,
        # length of the guideRNAs to design
        consensus_range = consensus_df(),
        # restrict guideRNA design to parts on consensus
        gc_content = input$gc_content,
        # allowed guide GC content (between 40 and 80 percent)
        alpha = input$alpha,
        five_prime_seq = input$five_prime_seq,
        n_clust = 12,
        force = TRUE
      )
      return(gs3_temp)
      message("Check: end addGuides")
    })
    
  })
  
  output$guides <- renderPlot({
    req(gs3())
    #plot guides
    withProgress(message = "Plot guides...", {
      plotGuides(gs3())
    })
    
  })
  
  
  #observer events for add combinations
  observeEvent(input$ques_alpha_combinations, {
    showModal(
      modalDialog(
        title = "Help",
        "Select the start and end position on the consensus model of your target repeats. Guides binding outside the consensus range will be scored neutrally.
      This allows you to select guides that have a preferential binding in certain areas of the repeat, e.g. upstream of potential promoters, indicated by enrichment
      of Cis regulatory elements.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_alpha_combinations, {
    showModal(
      modalDialog(
        title = "Help",
        "Numeric. Off-target score coefficient. Large alpha penalizes combinations with high off-target score while alpha = 0
        ignores off-targets and picks combinations with highest on-target binding.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_max_guides, {
    showModal(
      modalDialog(
        title = "Help",
        "Maximum number of distinct guides to consider when calculating combinations. Do not use a higher number than experimentally feasible.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  observeEvent(input$ques_iterations, {
    showModal(
      modalDialog(
        title = "Help",
        "Selecting the optimal set of guideRNAs for maximal on- and minimal off-targeting is not trivial.
        Repguide first computes the binding profiles of all possible combinations of previously selected guide representatives,
        i.e. the best guides per cluster. The best combination then serves as seed to initialize a greedy algorithm to further optimize the set.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  
  #calculate guide combinations
  gs4 <-  eventReactive(input$action_combination, {
    req(gs3())
    #browser()
    withProgress(message = "Add guide combinations...", {
      message("Check: start addCombinations")
      gs4_temp <- addCombinations(
        gs3(),
        # our guideSet
        #iterations = 10,
        # number of greedy search iterations
        greedy = TRUE,
        # run greedy algorithm
        alpha = input$alpha_combinations,
        # off-target score penalty coefficient
        max_guides = input$max_guides,
        force = TRUE
      )  # maximal number of guides to consider
      message("Check: end addCombinations")
      return(gs4_temp)
    })
  })
  
  #plot guide combination
  output$combinations <- renderPlot({
    req(gs4())
    withProgress(message = "Plot guide combinations...", {
      #plot guides
      plotCombinations(gs4())
    })
  })
  
  #prepare download data
  output$downloadData <- downloadHandler(
    filename = "Repguide_Results.zip",
    content = function(file) {
      withProgress(message = "Writing Files to Disk. Please wait...", {
        temp <- setwd(tempdir())
        on.exit(setwd(temp))
        files <- c("Repguide_output")
        dir.create("Repguide_output")
        dir.create(file.path("Repguide_output", "full_stats"))
        if (length(gs4()@targets) > 0) {
          data.table::fwrite(
            as_tibble(gs4()@targets),
            file = paste0('Repguide_output/full_stats/',
                          'targets.txt'),
            sep = '\t'
          )
        }
        # export kmers
        library(tidyr)
        library(ggplot2)
        kmers <- as_tibble(gs4()@kmers)
        if (length(kmers) > 0)
        {
          data.table::fwrite(
            kmers,
            file = paste0('Repguide_output/full_stats/', 'kmers.txt'),
            sep = '\t'
          )
        }
        
        # export combinations
        if (length(gs4()@combinations) > 0)
        {
          data.table::fwrite(
            unnest(gs4()@combinations),
            file = paste0('Repguide_output/full_stats/',
                          'combinations.txt'),
            sep = '\t'
          )
        }
        
        #export guides
        combinations_subs <-
          gs4()@combinations %>% filter(best) %>% unnest
        kmers <- as_tibble(gs4()@kmers)
        kmers_subs <-
          inner_join(kmers, combinations_subs, by = 'kmer_id')
        
        kmers_by_nguides <-
          kmers_subs %>%
          group_by(n_guides) %>%
          do(data.frame = as_tibble(.))
        
        kmer_seqs_by_nguides <-
          kmers_subs %>%
          group_by(n_guides) %>%
          do(guide_seq = unique(DNAStringSet(
            structure(.$guide_seq, names = as.character(.$kmer_id)),
            use.names = TRUE
          )))
        
        for (i in kmers_by_nguides$n_guides)
        {
          data.table::fwrite(
            kmers_by_nguides$data.frame[[i]],
            file = paste0("Repguide_output/", i, '_guides_binding.txt'),
            sep = '\t'
          )
          Biostrings::writeXStringSet(
            kmer_seqs_by_nguides$guide_seq[[i]],
            filepath = paste0("Repguide_output/",
                              i,
                              '_guides_sequence.fasta'),
            format = 'fasta'
          )
        }
        
        # export target plots
        i <- 1
        dpi <- 300
        gs_plot <- plotTargets(gs4())
        gs_plot <- plotGuides(gs_plot)
        gs_plot <- plotCombinations(gs_plot)
        
        for (p in gs_plot@plots$targets)
        {
          fn <- paste0("Repguide_output/", 'targets', '_', i, '.png')
          ggsave(fn, p, device = 'png', dpi = dpi)
          i <- i + 1
        }
        
        # export guide plots
        i <- 1
        for (p in gs_plot@plots$guides)
        {
          fn <- paste0("Repguide_output/", 'guides', '_', i, '.png')
          ggsave(fn, p, device = 'png', dpi = dpi)
          i <- i + 1
        }
        
        # export guide plots
        i <- 1
        for (p in gs_plot@plots$combinations)
        {
          fn <- paste0("Repguide_output/",
                       'combinations',
                       '_',
                       i,
                       '.png')
          ggsave(fn, p, device = 'png', dpi = dpi)
          i <- i + 1
        }
        saveRDS(gs_plot,
                paste0(
                  "Repguide_output/",
                  'combinations',
                  "Repguide_results.RDS"
                ))
        zip(zipfile = file, files = files)
        
      })
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)

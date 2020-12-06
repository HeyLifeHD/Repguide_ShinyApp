#Joschka Hey
#Shiny web application to run the R library Repguide

# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

#libraries
if(!"shiny"%in%  installed.packages()) install.packages("shiny")
library(shiny)
if(!"Repguide"%in%  installed.packages()){
    if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
    remotes::install_github(repo = 'tanaylab/repguide')
}
library(Repguide)
if(!"data.table"%in%  installed.packages()) install.packages("data.table")
library(data.table)

# Define UI for application that draws a histogram
ui <- fluidPage(# Application title
    #Logo of Repguide
    titlePanel(
        title = div(
            img(
                src = "logo.png",
                height = 200,
                width = 200
            ),
            "Repguide - design of guideRNAs for CRISPR/dCas9 targeting of repetitive DNA sequences"
        )
    ),
    
    # Sidebar input selection
    sidebarLayout(
        sidebarPanel(
            h4("Specify Repguide options for Target exploration"),
            selectInput(
                "ref_genome",
                label = h6("Select reference genome:"),
                choices = list(
                    "Human GRCh38/hg38" = "hg38",
                    "Human GRCh37/hg19" = "hg19",
                    "Mouse GRCh38/mm10" = "mm10"
                ),
                selected = 2
            ),
            textInput(
                "target_repeats",
                label = h6("Select target class of repeats:"),
                value = NULL,
                placeholder = "Enter repeat class..."
            ),
            textInput(
                "whitelist_repeats",
                label = h6("Select repeat class or classes to be whitelisted:"),
                value = NULL,
                placeholder = "Enter repeat class ..."
            ),
            sliderInput(
                "gap_frequencies",
                label = h6("Maximum gap frequency for target alignment:"),
                min = 0,
                max = 1,
                value = 0.8
            ),
            actionButton("action_target", "Submit target specifications"),
            
            h4("Specify Repguide options for guide generation"),
            h5("Restrict guideRNA design to parts on consensus sequence"),
            textInput(
                "start_position",
                label = h6("Relative start position:"),
                value = NULL,
                placeholder = "Enter start position, e.g. 50"
            ),
            textInput(
                "end_position",
                label = h6("Relative end position:"),
                value = NULL,
                placeholder = "Enter end position, e.g. 600"
            ),
            sliderInput(
                "guide_length",
                label = h6("Basepair size of the guideRNAs:"),
                min = 12,
                max = 26,
                value = 16
            ),
            sliderInput(
                "n_mismatches",
                label = h6(
                    "Maximal number of tolerated mismatches when assessing guideRNA binding targets:"
                ),
                min = 0,
                max = 3,
                value = 0
            ),
            textInput(
                "five_prime_seq",
                label = h6("Sequence requirement for 5' start of guideRNAs:"),
                value = NULL,
                placeholder = "Enter nucleotide, e.g. G for transcription from U6 promoter"
            ),
            sliderInput(
                "gc_content",
                label = h6("Allowed GC content for guidesRNAs:"),
                min = 0,
                max = 1,
                value = c(0.4, 0.8)
            ),
            sliderInput(
                "min_Son",
                label = h6("Minimal on target score of guides:"),
                min = 0,
                max = 100,
                value = 50
            ),
            sliderInput(
                "alpha",
                label = h6(
                    "Off-target score coefficient (large alpha penalizes guides with high off-target score):"
                ),
                min = 0,
                max = 100,
                value = 0
            ),
            h4("Specify Repguide options for guide combination"),
            sliderInput(
                "max_guides",
                label = h6(
                    "Maximum number of distinct guides to consider when calculating combinations:"
                ),
                min = 0,
                max = 30,
                value = 5
            ),
            sliderInput(
                "iterations",
                label = h6("Number of greedy search iterations:"),
                min = 0,
                max = 50,
                value = 10
            ),
            sliderInput(
                "alpha_combinations",
                label = h6(
                    "Off-target score coefficient (large alpha penalizes combinations with high off-target score):"
                ),
                min = 0,
                max = 100,
                value = 10
            )
            
        ),
        #main panel
        mainPanel(
            h1("Repguide workflow:"),
            #Image ofRepguide workflow
            img(
                src = "schematic.png" ,
                height = 280,
                width = 800
            ),
            #plot first output: target repeats
            h1("Target exploration:"),
            plotOutput("targets_repeats", height="600px"),
            h1("guideRNA design:"),
            plotOutput("guides", height="800px"),
            h1("Combinatorial optimization:"),
            plotOutput("combinations", height="1000px")
        )
    ))



# Define server logic required to draw a histogram
server <- function(input, output) {
    #select organism
    org <- reactive({
        if (input$ref_genome %in% c("hg38", "hg19")) {
            return("Hs")
        } else if (input$ref_genome %in% c("mm10", "mm9")) {
            return("Mm")
        } else{
            print("No appropriate reference genome chosen")
        }
        #load packages
        if(!paste0("org.", org, ".eg.db")%in%  installed.packages()) install.packages(paste0("org.", org, ".eg.db"))
        library(paste0("org.", org, ".eg.db"), character.only = TRUE)
    })
    
    #get appropriate reference genome
    BS_genome <- reactive({
        if(!"BSgenome"%in%  installed.packages()) install.packages(paste0("BSgenome"))
        if (input$ref_genome %in% c("hg38")) {
            if(!"BSgenome.Hsapiens.UCSC.hg38" %in%  installed.packages()) install.packages("BSgenome.Hsapiens.UCSC.hg38")
        } else if (input$ref_genome %in% c( "hg19")) {
            if(!"BSgenome.Hsapiens.UCSC.hg19" %in%  installed.packages()) install.packages("BSgenome.Hsapiens.UCSC.hg19")
        } else if (input$ref_genome %in% c( "mm10")) {
            if(!"BSgenome.Mmusculus.UCSC.mm10" %in%  installed.packages()) install.packages("BSgenome.Mmusculus.UCSC.mm10")
        } else{
            print("No appropriate BSgenome object found")
        }
        BSgenome::getBSgenome(genome = input$ref_genome)
    })
    
    #get appropriate txdb file
    txdb <- reactive({
        grep(
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
            colnames(txdb_objects[, txdb_objects[4, ] == input$ref_genome, drop = FALSE])
        if(!txdb %in%  installed.packages()) install.packages(txdb)
        library(txdb, character.only = TRUE)
        return(txdb)
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
        if (!is.null(input$whitelist_repeats)) {
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
            whitelist_regions <-
                repeats_dt[repeats_dt$repName %in% input$whitelist_repeats,]
        } else{
            whitelist_regions <- NULL
        }
        return(whitelist_regions)
    })
    
    #create guideset
    gs <- reactive({
        message("Check: start createGuideSet")
        gs_temp <- createGuideSet(
            genome = BS_genome(),
            alt_chromosomes = FALSE,
            #input$alt_chromosomes,
            tes = repeats_path(),
            temp = tempdir(),
            # Directory for temporary files
            cis = fantom_path(),
            blacklist = NULL,
            #input$blacklist,
            whitelist = whitelist_regions(),
            n_cores = 12,
            refdir = "",
            seed = 19
        )
        return(gs_temp)
        message("Check: end createGuideSet")
        
    })
    gs1 <- reactive({
        req(gs())
        message("Check: start addTargets")
        gs1_temp <- addTargets(gs(),
                               targets = input$target_repeats,
                               force = TRUE)
        return(gs1_temp)
        message("Check: end addTargets")
    })
    
    gs2 <- reactive({
        req(gs1())
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
    
    output$targets_repeats <- renderPlot({
        req(gs2())
        #plot target repeat class
        plotTargets(gs2())
    })
    
    #get consensus sequence for guide design
    consensus_df <- reactive({
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
    
    output$guides <- renderPlot({
        req(gs3())
        #plot guides
        plotGuides(gs3())
        
    })
    
    #calculate guide combinations
    gs4 <-  reactive({
        req(gs3())
        message("Check: start addCombinations")
        #browser()
        gs4_temp <- addCombinations(
            gs3(),
            # our guideSet
            #iterations = input$iterations,
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
    
    #plot guide combination
    output$combinations <- renderPlot({
        req(gs4())
        #plot guides
        plotCombinations(gs4())
    })
    
    
}


# Run the application
shinyApp(ui = ui, server = server)

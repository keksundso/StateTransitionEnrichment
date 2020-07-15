report: "report/workflow.rst"
IDS, = glob_wildcards("data/GOTerms/{GOList}.go")


rule report:
    input:
        "output/HKGList_BinomialModle_summaryPlot.pdf",
        expand("output/{GOList}_BinomialModle_summaryPlot.pdf", GOList=IDS),
        #
        "output/HKGList_BinomialModle_PCC.pdf",
        expand("output/{GOList}_BinomialModle_PCC.pdf", GOList=IDS),
        #
        "output/Sim_BinomialModle_PCC.pdf",
        "output/Sim_BinomialModle_summaryPlot.pdf"   ,   
        "output/recoverdSim.pdf",
        #
        expand("output/{GOList}_posteriorSummaryExplain.pdf", GOList=IDS),
        #
        expand("output/GOAnalysis/{GOList}/{GOList}_GOEnrichment.pdf", GOList=IDS)

    #output:
    #    report("dag.svg", category="0. A Workflow Overview")
        
    shell:
        "snakemake --dag -np | dot -Tsvg > dag.svg"

rule creatingGOList:
    input:
        GOListName="data/GOTerms/{GOList}.go"
    output:
        GOVector="data/ListOfGenes/{GOList}.rda"
    script:
        "scripts/creatGOList.R"


rule creatingHKGList:
    input:
        HGcsv="data/ListOfGenes/journal.pone.0000898.s001.csv"
    output:
        HGKVector="data/ListOfGenes/HKGList.rda"
    script:
        "scripts/creatingHKGList.R"


rule creatingSimHitTable:
    output:
        plotCreat="ZS/SimCreation.pdf",
        SimCompVar="ZS/Sim_completeDF.rda",
        propDFVar="ZS/SimPropDF.rda"
    script:
        "scripts/creatSimData.R"



rule creatingHitTable:
    input:
        mutHitDF="data/Jan2020gcsc3.mutations.RData",
        GeneSelection="data/ListOfGenes/{GeneListWK}.rda"
    output:
        compDFName="ZS/{GeneListWK}_completeDF.rda"
    script:
        "scripts/ReadingInMutData.R"

rule recoveryOfSim:
    input:
        SimModlePath="ZS/Sim_BinomialModle.rda",
        propDFVar="ZS/SimPropDF.rda"
    output:
        simRecPath=report("output/recoverdSim.pdf", caption="report/figSimRec.rst", category="1.Model-Accuracy_Recovery")
    script:
        "scripts/recoveryOfSimDF.R"


rule BRMSModle:
    input:
        compDFName="ZS/{GeneListWK}_completeDF.rda"
    output:
        brmsModleVar="ZS/{GeneListWK}_BinomialModle.rda",
        priorTable="ZS/{GeneListWK}_PriorOfBinomialModle.csv"
    script:
        "scripts/brmsModle.R"


rule PPC:
    input:  
        compDFName="ZS/{GeneListWK}_completeDF.rda",
        brmsModleVar="ZS/{GeneListWK}_BinomialModle.rda"
    output:
        PPCVar=report("output/{GeneListWK}_BinomialModle_PCC.pdf", caption="report/figPPC.rst", category="2.Model-Accuracy_PPC")
    script:
        "scripts/PPC.R"



rule summaryPlot:
    input:
        compDFName="ZS/{GeneListWK}_completeDF.rda",
        brmsModleVar="ZS/{GeneListWK}_BinomialModle.rda"
    output:
        plotVar=report("output/{GeneListWK}_BinomialModle_summaryPlot.pdf", caption="report/figResults.rst", category="3.Result")
    script:
        "scripts/pVSOdds_plot.R"


rule explanationPlot:
    input:
        brmsModleVar="ZS/{GOList}_BinomialModle.rda"
    output:
        plotExp=report("output/{GOList}_posteriorSummaryExplain.pdf", caption="report/figposteriorSummaryExp.rst", category="4.Clarifications")
    script:
        "scripts/ExplainingPosteriorSummary.R"

rule GOEnrichment:
    input:
        mutHitDF="data/Jan2020gcsc3.mutations.RData",
        GOVector="data/ListOfGenes/{GOList}.rda",
        brmsModleVar="ZS/{GOList}_BinomialModle.rda"
    output:
        pdfEnrichment=report("output/GOAnalysis/{GOList}/{GOList}_GOEnrichment.pdf", caption="report/figGOEnrichment.rst", category="5.Interpretation")
    script:
        "scripts/GOEnrichment.R"

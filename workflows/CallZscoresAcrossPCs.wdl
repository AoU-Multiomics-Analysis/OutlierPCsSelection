version 1.0 

task ComputeZscoresScript {
    input {
        File GCTFile
        Int Memory
    }
    command <<<
    Rscript /tmp/ComputeFitAcrossPCs.R --GCTFile ~{GCTFile}  
    >>>
    
    output {
        File Zscores = "ZscoresAcrossPCs.rds" 
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/OutlierPCsSelection:main"
        preemptible: "1"
        cpu: "4"
        memory: "~{Memory} GB"
        disks: "local-disk 100 HDD"
    }

}


workflow ComputeZscores {
    input {
        File GCTFile
        Int Memory
    }
    
    call ComputeZscoresScript {
        input:
            GCTFile = GCTFile,
            Memory = Memory
    }

    output {
        File Zscores = ComputeZscoresScript.Zscores
    }
} 

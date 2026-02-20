version 1.0 

task ComputeZscoresScript {
    input {
        File GCTFile 
    }
    command <<<
    Rscript /tmp/ComputeFitAcrossPCs.R --GCTFile ~{GCTFile}  
    >>>
    
    output {
        File Zscores = "ZscoresAcrossPCs.rds" 
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        preemptible: "~{NumPrempt}"
        cpu: "4"
        memory: "~{Memory} GB"
        disks: "local-disk 100 HDD"
    }

}


workflow ComputeZscores {
    input {
        File GCTFile
    }
    
    call ComputeZscoresScript {
        input:
            GCTFile = GCTFile
    }

    output {
        File Zscores = ComputeZscoresScript.Zscores
    }
} 

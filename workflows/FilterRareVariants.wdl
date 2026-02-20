version 1.0 

task SubsetVariants {
    input {
        File VCF 
        File VCFIndex
        File CisWindows
        Int Memory
    }
    
    command <<<
    bcftools view -R ~{CisWindows} -Oz -o CisWindowRareVariants.gz ~{VCF} 
    >>>


    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/OutlierPCsSelection:main"
        preemptible: "1"
        cpu: "4"
        memory: "~{Memory} GB"
        disks: "local-disk 100 HDD"
    }
}


workflow FilterRareVariants {
    input {
        File VCF 
        File VCFIndex
        File CisWindows
        File GnomadAFs
        Int Memory
    }
    
    call SubsetVariants {
        input:
            VCF = VCF,
            VCFIndex = VCFIndex,
            CisWindows = CisWindows
            Memory = Memory
    }


} 

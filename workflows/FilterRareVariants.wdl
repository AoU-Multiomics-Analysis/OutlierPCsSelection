version 1.0 

task SubsetVariants {
    input {
        File VCF 
        File VCFIndex
        File CisWindows
    }
    
    command <<<
    bcftools view -R ~{CisWindows} -Oz -o CisWindowRareVariants.gz ~{VCF} 
    >>>
}


workflow FilterRareVariants {
    input {
        File VCF 
        File VCFIndex
        File CisWindows
        File GnomadAFs
    }
    
    call SubsetVariants {
        input:
            VCF = VCF,
            VCFIndex = VCFIndex,
            CisWindows = CisWindows
    }


} 
